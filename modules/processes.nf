

process Concatenate {
	label "small_cpus"
	label "small_mem"

	input:
	file "*"

	output:
	file "allfq.fastq"

	shell:
	"""
	cat * > allfq.fastq
	"""
}


process Demultiplex {
	label "big_cpus"
	label "small_mem"

	input:
	file fq

	output:
	file "porechop_results/*"

	shell:
	"""
	porechop -t ${task.cpus} -i ${fq} -b porechop_results
	rm -f porechop_results/none.fastq
	"""
}

process Filter {
	tag "$barcode_id"
	label "big_cpus"
	label "small_mem"

	input:
	tuple val(barcode_id), file("${barcode_id}.fastq")

	output:
	tuple val(barcode_id), file("Filt_${barcode_id}.fastq")
	// file "${x.baseName}" into nanofilted2

	shell:
	"""
	cat "${barcode_id}.fastq" | NanoFilt --quality ${params.nanofilt_quality} --maxlength ${params.nanofilt_maxlength} > Filt_${barcode_id}.fastq
	"""

}


process NanoPlotNoFilt {
	tag "$barcode_id"
	label "small_cpus"
	label "small_mem"

	publishDir "$params.outdir/NanoPlots", pattern: "NanoPlot*", mode: "copy"

	input:
	tuple val(barcode_id), file("${barcode_id}.fastq")

	output:
	tuple val(barcode_id), path("Nanoplot.NoFilt.${barcode_id}"), emit: nanoplot, optional: true
	tuple val(barcode_id), path("count_NoFilt_${barcode_id}.txt"), emit: counts

	script:
	"""
	COUNT=\$(echo \$(cat ${barcode_id}.fastq | wc -l)/4 | bc)
	TWO=2
	echo \$COUNT > count_NoFilt_${barcode_id}.txt
	if [ "\$COUNT" -gt "\$TWO" ]
	then
		NanoPlot -t ${task.cpus} --fastq ${barcode_id}.fastq -o Nanoplot.NoFilt.${barcode_id}
	fi
	"""
}


process NanoPlotFilt {
	tag "$barcode_id"
	label "small_cpus"
	label "small_mem"

	publishDir "$params.outdir/NanoPlots", pattern: "NanoPlot*", mode: "copy"

	input:
	tuple val(barcode_id), file("Filt_${barcode_id}.fastq")

	output:
	tuple val(barcode_id), path("Nanoplot.${barcode_id}"), emit: nanoplot, optional: true
	tuple val(barcode_id), path("count_Filt_${barcode_id}.txt"), emit: counts

	script:
	"""
	COUNT=\$(echo \$(cat Filt_${barcode_id}.fastq | wc -l)/4 | bc)
	TWO=2
	echo \$COUNT > count_Filt_${barcode_id}.txt
	if [ "\$COUNT" -gt "\$TWO" ]
	then
		NanoPlot -t ${task.cpus} --fastq Filt_${barcode_id}.fastq -o Nanoplot.${barcode_id}
	fi
	"""
}

process SummaryTable{
	label "small_cpus"
	label "small_mem"

	publishDir "$params.outdir/NanoPlots", mode: "copy"

	input:
	file "*"
	//file "*"

	output:
	file "summary.tsv"

	shell:
	"""
	#!/usr/bin/env Rscript
	filt_files <- list.files(pattern = "_Filt_")
	nofilt_files <- list.files(pattern = "_NoFilt_")
	filt_files <- setNames(filt_files, gsub("^count_Filt_|[.]txt\$", "", filt_files))
	nofilt_files <- setNames(nofilt_files, gsub("^count_NoFilt_|[.]txt\$", "", nofilt_files))

	lp <- lapply(list(RawReads = nofilt_files, Filtered = filt_files), function(x) {
		vapply(x, readLines, FUN.VALUE=NA_character_)
	})

	x <- do.call(cbind, lapply(lp, function(x) { x[match(names(lp[[1]]), names(x))] }) )
	x <- as.data.frame(x)
	x\$BarCode <- rownames(x)
	x <- x[, c("BarCode", "RawReads", "Filtered")]

	write.table(x, file = "summary.tsv", sep="\t", quote=F, row.names=FALSE, col.names=TRUE)
	"""
}

process Fastq2Fasta {
	tag "$barcode_id"
	label "big_cpus"
	label "big_mem"

	input:
	tuple val(barcode_id), path("${barcode_id}.fastq")

	output:
	tuple val(barcode_id), path("${barcode_id}.fasta")

	when:
	params.stoptocheckparams == false

	shell:
	"""
	paste - - - - < ${barcode_id}.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > ${barcode_id}.fasta
	"""
}


process LastAL {
	label "big_cpus"
	label "big_mem"
	tag "$barcode_id"
	publishDir "$params.outdir/LastAL", enabled: params.keepmaf, pattern: "*.maf", mode: "copy"

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta")

	output:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.maf")

	when:
	!params.stoptocheckparams

	shell:
	"""
	lastal -P${task.cpus} /opt/silva/silva ${barcode_id}.fasta > ${barcode_id}.maf

	"""
}

process DAAConverter{
	label "big_cpus"
	label "big_mem"
	tag "$barcode_id"
	cpus params.cpus

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.maf")

	output:
	tuple val(barcode_id), file("${barcode_id}.daa")

	when:
	!params.stoptocheckparams

	shell:
	"""
	java -jar /opt/DAA_Converter_v0.9.0.jar -top 20 -p ${task.cpus} -i ${barcode_id}.maf -r ${barcode_id}.fasta -o ${barcode_id}.daa
	"""
}

process DAAMeganizer{
	tag "$barcode_id"
	label "big_cpus"
	label "big_mem"

	input:
	tuple val(barcode_id), file("${barcode_id}.daa")

	output:
	tuple val(barcode_id), file("${barcode_id}.daa")

	when:
	params.stoptocheckparams == false

	shell:
	"""
	daa-meganizer -i ${barcode_id}.daa -p ${task.cpus} -s2t /opt/silva/SSURef_Nr99_138_tax_silva_to_NCBI_synonyms.map --lcaAlgorithm ${params.megan_lcaAlgorithm} --lcaCoveragePercent ${params.megan_lcaCoveragePercent}
	"""
}




process ComputeComparison{
	label "small_cpus"
	label "small_mem"

	publishDir "$params.outdir/Megan_Comparison", mode: "copy"

	input:
	file "*"

	output:
	file "comparison.megan"

	when:
	params.stoptocheckparams == false

	shell:
	"""
	xvfb-run --auto-servernum --server-num=1 /usr/local/bin/compute-comparison -i ./* -o comparison.megan
	"""
}



process ExtractOtuTable {
	label "small_cpus"
	label "small_mem"

	publishDir "$params.outdir/Megan_Comparison", mode: "copy"

	input:
	file "comparison.megan"

	output:
	file "OTU_Table.tsv"

	shell:
	"""
	#!/usr/bin/env Rscript
	rl <- readLines("comparison.megan")
	snam <- strsplit(grep("^@Names", rl, value = TRUE), "\t")[[1]][-1]

	taxs <- grep("^TAX", rl)
	taxs <- strsplit(rl[taxs], "\t")
	taxs <- lapply(taxs, "[", -1)
	names(taxs) <- lapply(taxs, "[", 1)
	taxs <- lapply(taxs, "[", -1)
	taxs <- lapply(taxs, as.integer)

	nsam <- length(snam)
	ln <- sapply(taxs, length)
	addz <- nsam - ln

	mp <- t(mapply(function(x, add){c(x, rep(0, add))}, x=taxs, add=addz))
	colnames(mp) <- snam

	write.table(mp, file = "OTU_Table.tsv", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
	"""
}

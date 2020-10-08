
process gunzip {
	label 'small_cpus'
	label 'small_mem'

	input:
	file "*"

	output:
	file "*"

	script:
	"""
	gunzip -f *.gz
	"""
}


process downloadSilvaFasta {
	label 'internet'
	label 'small_cpus'
	label 'small_mem'

	output:
	path("*.fasta")

	script:
	"""
	wget ${params.silvaFastaURL}
	gunzip *gz
	"""
}

process downloadSilvaTaxNcbiSp {
	label 'internet'
	label 'small_cpus'
	label 'small_mem'

	output:
	path("*.txt")

	script:
	"""
	wget ${params.silvaTaxNcbiSpURL}
	gunzip *gz
	"""
}

process downloadSilvaTaxmap {
	label 'internet'
	label 'small_cpus'
	label 'small_mem'

	output:
	path("*.txt")

	script:
	"""
	wget ${params.silvaTaxmapURL}
	gunzip *gz
	"""
}

process generateSynonyms {
	label 'small_cpus'
	label 'small_mem'

	input:
	file "tax_ncbi-species_ssu_ref_nr99_VERSION.txt"
	file "taxmap_slv_ssu_ref_nr_VERSION.txt"

	output:
	file "SSURef_Nr99_tax_silva_to_NCBI_synonyms.map"

	shell:
	"""
	#!/usr/bin/env Rscript

	taxmap <- "taxmap_slv_ssu_ref_nr_VERSION.txt"
	tms <- read.csv(taxmap, sep = "\\t")

	ncbisp <- "tax_ncbi-species_ssu_ref_nr99_VERSION.txt"
	ns <- read.csv(ncbisp, sep = "\\t", header = FALSE)

	# Clean taxmap organism_name
	cln <- grep("[U|u]ncultured|[E|e]nvironmental|[U|u]nidentified|[M|m]etagenome", tms\$organism_name)
	tms\$organism_name[cln] <- NA_character_

	# Clean taxmap path
	tms\$path <- gsub("uncultured;\$", "", tms\$path)

	# Create full path
	rps <- apply(tms, 1, function(x) grepl(x[["organism_name"]], x[["path"]], fixed = TRUE))
	tms\$organism_name[rps] <- ""
	tms\$organism_name[cln] <- ""
	tms\$path <- paste(tms\$path, tms\$organism_name, sep = "")

	# Identify silva indexes for occurrences
	idx <- split(seq_len(dim(tms)[1]), tms\$path)

	# Clean ncbi taxonomy
	ns\$V1 <- gsub(" <\\\\w+>", "", ns\$V1)

	# Create hashmap for ncbi taxonomy
	ns\$tip <- vapply(strsplit(ns\$V1, ";"), function(x) rev(x)[1], FUN.VALUE = NA_character_)
	ee <- list2env(split(ns\$V2, ns\$tip), hash = TRUE)

	# Find taxids 
	pths <- strsplit(setNames(names(idx), names(idx)), ";")
	pths <- lapply(pths, rev)
	txid <- lapply(pths, function(x){
	r <- NULL
	y <- x
	while(is.null(r)){
		a <- y[1]
		y <- y[-1]
		r <- ee[[a]]
		r <- r[1]
	}
	r
	})

	# Repeat taxids as many times as it was assigned
	rp <- mapply(rep, txid, lengths(idx))

	# Sort taxids by index
	ultx <- unlist(rp, use.names = FALSE)
	tms\$map <- ultx[order(unlist(idx, use.names = FALSE))]

	df <- data.frame(Accs = paste(tms\$primaryAccession, tms\$start, tms\$stop, sep = "."),
					Syno = tms\$map)

	write.table(df, 
				file = "SSURef_Nr99_tax_silva_to_NCBI_synonyms.map", 
				sep = "\\t", 
				quote = FALSE, 
				row.names = FALSE, 
				col.names = FALSE)
	"""
}

/*
process downloadMeganSynMap {
	label 'internet'
	label 'small_cpus'
	label 'small_mem'

	output:
	path("*")

	script:
	"""
	wget ${params.meganSynMapURL}
	gunzip *gz
	"""
}


process trimAccTaxID {
	label 'small_cpus'
	label 'small_mem'

	input:
	file "tax_slv_ssu_138.acc_taxid"


	output:
	path("SSURef_Nr99_tax_silva_to_NCBI_synonyms.map")


	script:
	"""
	cat tax_slv_ssu_138.acc_taxid | awk -F '[.\\t]' '{print \$1 "\\t" \$4}' > SSURef_Nr99_tax_silva_to_NCBI_synonyms.map
	"""
}
*/

process Concatenate {
	label "small_cpus"
	label "small_mem"

	input:
	file "*fastq"

	output:
	file "allfq.fastq"

	script:
	"""
	cat *fastq > allfq.fastq
	"""
}


process Demultiplex {
	label "big_cpus"
	label "big_mem"

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
	label "small_cpus"

	input:
	tuple val(barcode_id), file("${barcode_id}.fastq")

	output:
	tuple val(barcode_id), file("Filt_${barcode_id}.fastq")

	shell:
	"""
	cat "${barcode_id}.fastq" | NanoFilt --quality ${params.nanofilt_quality} --maxlength ${params.nanofilt_maxlength} --length ${params.nanofilt_length} > Filt_${barcode_id}.fastq
	"""

}


process NanoPlotRaw {
	tag "$barcode_id"
	label "small_cpus"
	label "small_mem"

	publishDir "$params.outdir/NanoPlots/${barcode_id}/Raw", mode: "copy"

	input:
	tuple val(barcode_id), file("${barcode_id}.fastq")

	output:
	tuple val(barcode_id), path("${barcode_id}_Raw_*"), emit: nanoplot, optional: true
	tuple val(barcode_id), path("count_Raw_${barcode_id}.txt"), emit: counts

	script:
	"""
	COUNT=\$(echo \$(cat ${barcode_id}.fastq | wc -l)/4 | bc)
	TWO=2
	echo \$COUNT > count_Raw_${barcode_id}.txt
	if [ "\$COUNT" -gt "\$TWO" ]
	then
		NanoPlot -t ${task.cpus} --fastq ${barcode_id}.fastq -p ${barcode_id}_Raw_
	fi
	"""
}


process NanoPlotFilt {
	tag "$barcode_id"
	label "small_cpus"
	label "small_mem"

	publishDir "$params.outdir/NanoPlots/${barcode_id}/Filtered", mode: "copy"

	input:
	tuple val(barcode_id), file("Filt_${barcode_id}.fastq")

	output:
	tuple val(barcode_id), path("${barcode_id}_Filtered_*"), emit: nanoplot, optional: true
	tuple val(barcode_id), path("count_Filt_${barcode_id}.txt"), emit: counts

	script:
	"""
	COUNT=\$(echo \$(cat Filt_${barcode_id}.fastq | wc -l)/4 | bc)
	TWO=2
	echo \$COUNT > count_Filt_${barcode_id}.txt
	if [ "\$COUNT" -gt "\$TWO" ]
	then
		NanoPlot -t ${task.cpus} --fastq Filt_${barcode_id}.fastq -p ${barcode_id}_Filtered_
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
	raw_files <- list.files(pattern = "_Raw_")
	filt_files <- setNames(filt_files, gsub("^count_Filt_|[.]txt\$", "", filt_files))
	raw_files <- setNames(raw_files, gsub("^count_Raw_|[.]txt\$", "", raw_files))

	lp <- lapply(list(RawReads = raw_files, Filtered = filt_files), function(x) {
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
	label "small_cpus"
	label "small_mem"

	input:
	tuple val(barcode_id), path("${barcode_id}.fastq")

	output:
	tuple val(barcode_id), path("${barcode_id}.fasta")

	when:
	params.stoptocheckparams == false

	shell:
	"""
	seqtk seq -A ${barcode_id}.fastq > ${barcode_id}.fasta
	"""
}

process MakeLastDB {
	label "big_cpus"

	input:
	path("silva_SSU_tax.fasta")

	output:
	file "silva.*"

	shell:
	"""
	lastdb -cR01 -P${task.cpus} silva silva_SSU_tax.fasta
	"""
}

process LastTrain {
	tag "$barcode_id"
	label "big_cpus"

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta")
	file "*"
	
	
	output:
	tuple val(barcode_id), path("${barcode_id}.par")
		

	shell:
	"""
	last-train -P${task.cpus} -Q0 silva ${barcode_id}.fasta > ${barcode_id}.par 
	"""
}

process LastAL {
	label "big_cpus"
	label "big_mem"
	tag "$barcode_id"
	publishDir "$params.outdir/LastAL", enabled: params.keepmaf, pattern: "*.maf", mode: "copy"

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta")
	path "*"

	output:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.maf")


	shell:
	"""
	lastal -P${task.cpus} silva ${barcode_id}.fasta > ${barcode_id}.maf
	"""
}

process LastALPar {
	label "big_cpus"
	label "big_mem"
	tag "$barcode_id"
	publishDir "$params.outdir/LastAL", enabled: params.keepmaf, pattern: "*.maf", mode: "copy"

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.par")
	path "*"

	output:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.maf")


	shell:
	"""
	lastal -P${task.cpus} -p ${barcode_id}.par silva ${barcode_id}.fasta > ${barcode_id}.maf
	"""
}

process Maf2Sam {
	label "big_cpus"
	label "small_mem"
	tag "$barcode_id"

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.maf")
	path("silva_SSU_tax.fasta")

	output:
	tuple val(barcode_id), path("${barcode_id}.sam"), path("${barcode_id}.fasta")

	shell:
	"""
	maf-convert sam ${barcode_id}.maf | samtools view --threads ${task.cpus} -uT silva_SSU_tax.fasta | samtools sort -n --threads ${task.cpus} -O sam -o ${barcode_id}.sam
	"""
}
/*
process DAAConverter{
	label "big_cpus"
	label "big_mem"
	tag "$barcode_id"

	input:
	tuple val(selected_wf), val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.maf")

	output:
	tuple val(selected_wf), val(barcode_id), file("${barcode_id}.daa")

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

	input:
	tuple val(selected_wf), val(barcode_id), file("${barcode_id}.daa")
	path("SSURef_Nr99_tax_silva_to_NCBI_synonyms.map")

	output:
	tuple val(selected_wf), file("${barcode_id}.daa")

	when:
	params.stoptocheckparams == false

	shell:
	"""
	daa-meganizer -i ${barcode_id}.daa -p ${task.cpus} -s2t SSURef_Nr99_tax_silva_to_NCBI_synonyms.map --lcaAlgorithm ${params.megan_lcaAlgorithm} --lcaCoveragePercent ${params.megan_lcaCoveragePercent}
	"""
}
*/



process ComputeComparison{
	tag "$selected_wf"
	label "small_cpus"
	label "small_mem"
	
	publishDir "$params.outdir/Megan_Comparison", mode: "copy"

	input:
	tuple val(selected_wf), file("*")

	output:
	tuple val(selected_wf), file("${selected_wf}_comparison.megan")

	when:
	params.stoptocheckparams == false

	shell:
	"""
	xvfb-run --auto-servernum --server-num=1 /usr/local/bin/compute-comparison -i ./* -o ${selected_wf}_comparison.megan -n ${params.normalizeOtu}
	"""
}



process ExtractOtuTable {
	tag "$selected_wf"
	label "small_cpus"
	label "small_mem"

	publishDir "$params.outdir/Megan_Comparison", mode: "copy"

	input:
	tuple val(selected_wf), file("${selected_wf}_comparison.megan")

	output:
	path("${selected_wf}_OTU_Table.tsv")

	when:
	params.stoptocheckparams == false

	shell:
	"""
	#!/usr/bin/env Rscript
	rl <- readLines("${selected_wf}_comparison.megan")
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
	if (dim(mp)[1] == 1){
		mp <- t(mp)
	}
	colnames(mp) <- snam
	write.table(mp, file = "${selected_wf}_OTU_Table.tsv", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
	"""
}


process MakeMinimapDB {
	label "big_cpus"

	input:
	path("silva_SSU_tax.fasta")

	output:
	path("silva_k${params.minimap2_k}.mmi")

	shell:
	"""
	minimap2 -t ${task.cpus} -k ${params.minimap2_k} -d silva_k${params.minimap2_k}.mmi silva_SSU_tax.fasta
	"""
}

process Minimap2 {
	tag "$barcode_id"
	label "big_cpus"

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta")
	path("silva_k${params.minimap2_k}.mmi")

	output:
	tuple val(barcode_id), path("${barcode_id}.sam"), path("${barcode_id}.fasta")

	shell:
	"""
	minimap2 -K ${params.minimap2_KM}M -t ${task.cpus} -ax ${params.minimap2_x} silva_k${params.minimap2_k}.mmi ${barcode_id}.fasta > ${barcode_id}.sam
	"""
}

process Sam2Rma {
	tag "$barcode_id"
	label "big_cpus"

	publishDir "$params.outdir/Rma", mode: "copy"

	input:
	tuple val(barcode_id), path("${barcode_id}.sam"), path("${barcode_id}.fastq")
	path("SSURef_Nr99_tax_silva_to_NCBI_synonyms.map")
	val(selected_wf)

	output:
	tuple val(selected_wf), path("${selected_wf}_${barcode_id}.rma")

	shell:
	"""
	sam2rma -i ${barcode_id}.sam -r ${barcode_id}.fastq -o ${selected_wf}_${barcode_id}.rma -lg -alg ${params.megan_lcaAlgorithm} -lcp ${params.megan_lcaCoveragePercent} --topPercent ${params.megan_topPercent} -ram readCount -s2t SSURef_Nr99_tax_silva_to_NCBI_synonyms.map
	"""
}

process MakeBlastDB {
	label "big_cpus"

	input:
	path("silva_SSU_tax.fasta")

	output:
	file "silva_SSU_tax.*"

	shell:
	"""
	makeblastdb -in silva_SSU_tax.fasta -out silva_SSU_tax -parse_seqids -dbtype nucl
	"""
}

process MegaBlast {
	tag "$barcode_id"
	label "big_cpus"

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta")
	path("*")

	output:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.tab")


	shell:
	"""
	blastn -task "megablast" -num_threads ${task.cpus} -db silva_SSU_tax -query ${barcode_id}.fasta -out ${barcode_id}.tab -outfmt 6
	"""
}

process Blast2Rma {
	tag "$barcode_id"
	label "big_cpus"

	publishDir "$params.outdir/Rma", mode: "copy"

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.tab")
	path("SSURef_Nr99_tax_silva_to_NCBI_synonyms.map")
	val(selected_wf)

	output:
	path("${selected_wf}_${barcode_id}.rma")

	when:
	params.stoptocheckparams == false

	shell:
	"""
	blast2rma -i ${barcode_id}.tab -f BlastTab -bm BlastN -r ${barcode_id}.fasta -o ${selected_wf}_${barcode_id}.rma -lg -alg ${params.megan_lcaAlgorithm} -lcp ${params.megan_lcaCoveragePercent} --topPercent ${params.megan_topPercent} -ram readCount -s2t SSURef_Nr99_tax_silva_to_NCBI_synonyms.map
	"""
}

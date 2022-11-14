
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


process Porechop {
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

process NanoFilt {
	tag "$barcode_id"
	label "small_cpus"

	input:
	tuple val(barcode_id), file("${barcode_id}.fastq")

	output:
	tuple val(barcode_id), file("Filt_${barcode_id}.fastq"), optional: true

	shell:
	"""
	cat "${barcode_id}.fastq" | \
		NanoFilt \
		--quality ${params.nanofilt_quality} \
		--maxlength ${params.nanofilt_maxlength} \
		--length ${params.nanofilt_length} > \
		Filt_${barcode_id}.fastq

	[ -s Filt_${barcode_id}.fastq ] || rm Filt_${barcode_id}.fastq
	"""

}


process AutoMap {
	tag "$barcode_id"
	label "big_cpus"

	input:
	tuple val(barcode_id), file("Filt_${barcode_id}.fastq")

	output:
	tuple val(barcode_id), file("Filt_${barcode_id}.fastq"), file("overlap_${barcode_id}.paf"), optional: true

	shell:
	"""
	minimap2 \
		-x ava-ont \
		-t ${task.cpus} \
		-g 500 \
		-f${params.minimap2_f} \
		Filt_${barcode_id}.fastq Filt_${barcode_id}.fastq > overlap_${barcode_id}.paf
		[ -s overlap_${barcode_id}.paf ] || rm overlap_${barcode_id}.paf
	"""
}

process Yacrd {
	tag "$barcode_id"
	label "big_cpus"

	input:
	tuple val(barcode_id), file("Filt_${barcode_id}.fastq"), file("overlap_${barcode_id}.paf")

	output:
	tuple val(barcode_id), file("Filt_Scrubb_${barcode_id}.fastq")

	shell:
	"""
	yacrd \
		-i overlap_${barcode_id}.paf \
		-o report_${barcode_id}.yacrd \
		-c ${params.yacrd_c} \
		-n ${params.yacrd_n} \
		scrubb \
		-i Filt_${barcode_id}.fastq \
		-o Filt_Scrubb_${barcode_id}.fastq
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
	tuple val(barcode_id), file("Filt_Scrubb_${barcode_id}.fastq")

	output:
	tuple val(barcode_id), path("${barcode_id}_Filtered_Scrubbed_*"), emit: nanoplot, optional: true
	tuple val(barcode_id), path("count_Filt_Scrubb_${barcode_id}.txt"), emit: counts

	script:
	"""
	COUNT=\$(echo \$(cat Filt_Scrubb_${barcode_id}.fastq | wc -l)/4 | bc)
	TWO=2
	echo \$COUNT > count_Filt_Scrubb_${barcode_id}.txt
	if [ "\$COUNT" -gt "\$TWO" ]
	then
		NanoPlot -t ${task.cpus} --fastq Filt_Scrubb_${barcode_id}.fastq -p ${barcode_id}_Filtered_Scrubbed_
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
	filt_files <- list.files(pattern = "_Filt_Scrubb_")
	raw_files <- list.files(pattern = "_Raw_")
	filt_files <- setNames(filt_files, gsub("^count_Filt_Scrubb_|[.]txt\$", "", filt_files))
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

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta")
	path "*"

	output:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.tab")

	shell:
	"""
	lastal -f BlastTab -E ${params.last_E} -P${task.cpus} silva ${barcode_id}.fasta > ${barcode_id}.tab
	"""
}

process LastALPar {
	label "big_cpus"
	label "big_mem"
	tag "$barcode_id"

	input:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.par")
	path "*"

	output:
	tuple val(barcode_id), path("${barcode_id}.fasta"), path("${barcode_id}.tab")

	shell:
	"""
	lastal -f BlastTab -E ${params.last_E} -P${task.cpus} -p ${barcode_id}.par silva ${barcode_id}.fasta > ${barcode_id}.tab
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
	minimap2 \
		-t ${task.cpus} \
		-k ${params.minimap2_k} \
		-d silva_k${params.minimap2_k}.mmi \
		silva_SSU_tax.fasta
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
	minimap2 \
		-K ${params.minimap2_KM}M \
		-t ${task.cpus} \
		--secondary=no\
		-ax ${params.minimap2_x} \
		silva_k${params.minimap2_k}.mmi \
		${barcode_id}.fasta > ${barcode_id}.sam
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
	tuple val(selected_wf), val(barcode_id), path("${selected_wf}_${barcode_id}.rma")

	shell:
	"""
	sam2rma \
		-i ${barcode_id}.sam \
		-r ${barcode_id}.fastq \
		-o ${selected_wf}_${barcode_id}.rma \
		-lg \
		-alg ${params.megan_lcaAlgorithm} \
		-lcp ${params.megan_lcaCoveragePercent} \
		--topPercent ${params.megan_topPercent} \
		--minPercentReadCover ${params.megan_minPercentReadCover} \
		-ram readCount \
		-s2t SSURef_Nr99_tax_silva_to_NCBI_synonyms.map
	"""
}

process Rma2Info {
	tag "$barcode_id"
	label "small_cpus"

	input:
	tuple val(selected_wf), val(barcode_id), path("${selected_wf}_${barcode_id}.rma")

	output:
	tuple val(selected_wf), path("${selected_wf}_${barcode_id}.info")

	shell:
	"""
	rma2info \
		-i ${selected_wf}_${barcode_id}.rma \
		-c2c Taxonomy \
		-p -mro \
		-o ${selected_wf}_${barcode_id}.info
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
	blastn -task "megablast" -evalue ${params.megablast_evalue} -num_threads ${task.cpus} -db silva_SSU_tax -query ${barcode_id}.fasta -out ${barcode_id}.tab -outfmt 6
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
	tuple val(selected_wf), val(barcode_id), path("${selected_wf}_${barcode_id}.rma")

	shell:
	"""
	blast2rma \
		-i ${barcode_id}.tab \
		-f BlastTab \
		-bm BlastN \
		-r ${barcode_id}.fasta \
		-o ${selected_wf}_${barcode_id}.rma \
		-lg \
		-alg ${params.megan_lcaAlgorithm} \
		-lcp ${params.megan_lcaCoveragePercent} \
		--topPercent ${params.megan_topPercent} \
		--minPercentReadCover ${params.megan_minPercentReadCover} \
		-ram readCount \
		-s2t SSURef_Nr99_tax_silva_to_NCBI_synonyms.map
	"""
}

process MergeResults{
	tag "$selected_wf"
	label "small_cpus"
	label "small_mem"
	
	publishDir "$params.outdir/Merged_Results", mode: "copy"

	input:
	tuple val(selected_wf), file("*")

	output:
	tuple file("${selected_wf}_COUNTS.tsv"), file("${selected_wf}_TAXCLA.tsv")

	shell:
	"""
	#!/usr/bin/env Rscript
	
	# Read *.info files
	infos <- list.files(pattern = "[.]info\$")
	lp <- lapply(setNames(infos, infos), function(x) {
	a <- try(read.csv(x, sep = "\\t", header = FALSE))
	if (class(a) != "try-error"){
		setNames(a\$V2, a\$V1)
	}else{
		NULL
	}
	})

	# Get present taxa
	lvls <- unique(unlist(lapply(lp, names)))

	# Create new names (TAXA_*)
	lnlv <- length(lvls)
	wdth <- nchar(lnlv)
	bc <- paste0("TAXA_", formatC(seq_along(lvls),
								width = wdth,
								format = 'd',
								flag = '0'))

	# Parse tax paths for each TAXA
	rr <- data.frame(otu = bc, raw = lvls)
	spl <- strsplit(setNames(rr\$raw, rr\$otu), ";")
	rks <- c("[D]", "[P]", "[C]", "[O]", "[F]", "[G]", "[S]")
	tax <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
	rks <- setNames(rks, tax)
	tt <- lapply(spl, function(x) {
	unlist(lapply(rks, function(y){
		grep(y, x, fixed = TRUE, value = TRUE)
	}))
	})
	tt <- lapply(tt, function(y) gsub("[ ]*\\\\[\\\\w{1,2}\\\\][ ]", "", y))

	# Create TAX table
	taxt <- matrix(NA_character_, 
				nrow = length(tt), 
				ncol = length(tax), 
				dimnames = list(names(tt), tax))
	for (i in seq_along(tt)){
	taxt[names(tt)[i], names(tt[[i]])] <- tt[[i]]
	}


	re <- as.data.frame(cbind(rr, taxt))

	# Rename OTU counts
	lp <- lapply(lp, function(x){
		if (length(x)){
			setNames(x, re\$otu[match(names(x), re\$raw)])
		}
	})

	# Create "OTU" (TAXA) Table
	otut <- matrix(0L, nrow = length(re\$otu), ncol = length(lp),
				dimnames = list(re\$otu, names(lp)))
	for (i in seq_along(lp)){
	otut[names(lp[[i]]), names(lp)[i]] <- lp[[i]]
	}
	colnames(otut) <- sub("[.]info\$", "", colnames(otut))


	# Write Taxa Classification Table and Counts Table
	write.table(taxt, "${selected_wf}_TAXCLA.tsv", quote = FALSE, sep = "\\t")
	write.table(otut, "${selected_wf}_COUNTS.tsv", quote = FALSE, sep = "\\t")
	"""
}
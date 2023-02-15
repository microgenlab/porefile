
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


process reduceSilva {
	label 'small_cpus'
	label 'smal_mem'

	input:
	path("full_silva.fasta")

	output:
	path("reduced_silva.fasta")

	script:
	"""
	seqkit --threads 1 grep -r -i -n \
		-p Eukaryota \
		-p Phage \
		-v \
		full_silva.fasta > reduced_silva.fasta  
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

	publishDir "$params.outdir/", mode: "copy"

	input:
	file "tax_ncbi-species_ssu_ref_nr99_VERSION.txt"
	file "taxmap_slv_ssu_ref_nr_VERSION.txt"

	output:
	file "silva_to_NCBI_synonyms.map"

	script:
	"""
	generateSynonyms.R \
		--ncbisp tax_ncbi-species_ssu_ref_nr99_VERSION.txt \
		--taxmap taxmap_slv_ssu_ref_nr_VERSION.txt \
		--out silva_to_NCBI_synonyms.map
	"""
}

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
	porechop \
		-t ${task.cpus} \
		--extra_end_trim ${params.porechop_extra_end_trim} \
		-i ${fq} \
		-b porechop_results
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
		--headcrop ${params.nanofilt_headcrop} \
		--tailcrop ${params.nanofilt_tailcrop} \
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
	summarizeQC.R \
		--prefix_filt_scrubb count_Filt_Scrubb_ \
		--prefix_raw count_Raw_ \
		--out summary.tsv
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
	seqkit fq2fa --threads 1 ${barcode_id}.fastq -o ${barcode_id}.fasta
	"""
}

process MakeDB {
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
	tuple val(barcode_id), path("${barcode_id}.sam"), path("${barcode_id}.fasta"), optional: true

	shell:
	"""
	minimap2 \
		-K ${params.minimap2_KM}M \
		-t ${task.cpus} \
		--secondary=no\
		-ax ${params.minimap2_x} \
		silva_k${params.minimap2_k}.mmi \
		${barcode_id}.fasta > ${barcode_id}.sam

	if [ "\$(samtools view -c -F 0x04 ${barcode_id}.sam)" -eq "0" ]; then 
		rm ${barcode_id}.sam; 
	fi
	"""
}

process MeganLca {
	tag "$barcode_id"
	label "big_cpus"

	input:
	tuple val(barcode_id), path("${barcode_id}.sam"), path("${barcode_id}.fastq")
	path("SSURef_Nr99_tax_silva_to_NCBI_synonyms.map")

	output:
	tuple val(barcode_id), path("${barcode_id}.rma")

	shell:
	"""
	sam2rma \
		-i ${barcode_id}.sam \
		-r ${barcode_id}.fastq \
		-o ${barcode_id}.rma \
		-lg \
		-alg ${params.megan_lcaAlgorithm} \
		-lcp ${params.megan_lcaCoveragePercent} \
		--topPercent ${params.megan_topPercent} \
		--minPercentReadCover ${params.megan_minPercentReadCover} \
		-ram readCount \
		-s2t SSURef_Nr99_tax_silva_to_NCBI_synonyms.map
	"""
}


process GetReadInfo {
	tag "$barcode_id"
	label "small_cpus"

	input:
	tuple val(barcode_id), path("${barcode_id}.rma")

	output:
	tuple val(barcode_id), path("${barcode_id}.read_info")

	shell:
	"""
	rma2info \
		-i ${barcode_id}.rma \
		-r2c Taxonomy \
		-n False -r -mro \
		-o ${barcode_id}.ncbi

	rma2info \
		-i ${barcode_id}.rma \
		-r2c Taxonomy \
		-p -mro \
		-o ${barcode_id}.path

	join -t \$'\\t' \
		<(sort ${barcode_id}.ncbi) \
		<(sort ${barcode_id}.path) \
		> ${barcode_id}.read_info
	"""
}

process ComputeAbundances {
	label "small_cpus"
	label "small_mem"

	input:
	path("*")
	path("synonyms.txt")
	val(polish)

	output:
	path("COUNTS.tsv"), emit: counts
	path("TAXCLA.tsv"), emit: taxcla
	path("high_abundance_silva_ids.txt"), emit: silva_ids, optional: true
	path("*.read_txt"),  emit: read_ids, optional: true

	script:
	"""
	computeAbundances.R \
		--lowAbundanceThreshold ${params.lowAbundanceThreshold} \
		--in_suffix read_info \
		--synonyms synonyms.txt \
		--out_counts COUNTS.tsv \
		--out_taxcla TAXCLA.tsv \
		--polish ${polish} \
		--out_silva high_abundance_silva_ids.txt \
		--out_suffix read_txt
	"""
}

process SubsetSilva{
	label "small_cpus"
	label "small_mem"

	input:
	path("reduced_silva.fasta")
	path("high_abundance_silva_ids.txt")

	output:
	path("subset_silva.fasta")

	script:
	"""
	seqkit grep -f \
		high_abundance_silva_ids.txt reduced_silva.fasta > \
		subset_silva.fasta
	"""
}

process SubsetReads{
	tag "${barcode_id}"
	label "small_cpus"
	label "small_mem"

	input:
	tuple val(barcode_id), path("${barcode_id}.read_txt"), path("filetered_${barcode_id}.fasta")

	output:
	tuple val(barcode_id), path("lowAbundance_${barcode_id}.fasta")

	script:
	"""
	seqkit grep -f \
		${barcode_id}.read_txt filetered_${barcode_id}.fasta > \
		lowAbundance_${barcode_id}.fasta
	"""
}

process CorrectAssignment{
	tag "${barcode_id}"
	label "small_cpus"
	label "small_mem"

	input:
	tuple val(barcode_id), path("${barcode_id}_old.read_info"), path("${barcode_id}_new.read_info")

	output:
	tuple val(barcode_id), path("${barcode_id}_corrected.read_info")

	script:
	"""
	correctAssignments.R \
		--old ${barcode_id}_old.read_info \
		--new ${barcode_id}_new.read_info \
		--corrected ${barcode_id}_corrected.read_info
	"""
}
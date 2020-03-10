#!/usr/bin/env nextflow

params.fq = "$baseDir/data/*.fastq"
params.cpus = 4
params.stoptocheckparams = false
params.nanofilt_quality = 8
params.nanofilt_maxlength = 1500
params.silva_fasta = "$baseDir/SILVA138_SSU/SILVA_138_SSURef_NR99_tax_silva.fasta"
params.megandb_map = "$baseDir/MEGANDB/megan-map-Oct2019.db"
params.megandb_nucl = "$baseDir/MEGANDB/megan-nucl-Oct2019.db"
params.megan_lcaAlgorithm = "naive"
params.megan_lcaCoveragePercent = 100
params.help = false





def helpMessage() {
    log.info """
	--------------------------
	---> Long16S Pipeline <---
	--------------------------
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run iferres/long16S --fq './data/*.fastq' --silva_fasta './SILVA138_SSU/SILVA_138_SSURef_NR99_tax_silva.fasta' --megabdb_map './MEGANDB/megan-map-Oct2019.db' --megandb_nucl './MEGANDB/megan-nucl-Oct2019.db' --cpus 4 -profile local
    Mandatory arguments:
	--fq                          Path to input data (must be surrounded with quotes).
	--silva_fasta                 Path to silva fasta file.
	--megandb_map                 Path to megan-map-Oct2019.db file (date on filename may change).
	--megandb_nucl                Path to megan-nucl-Oct2019.db file (date on filename may change).
	-profile                      Configuration profile to use. Available: local, nagual.

    Other:
	--cpus                        The max number of cpus to use on each process (Default: 4).
	--stoptocheckparams           Whether to stop after Summary process to check parameters. Default: false.
                                      If true, then the pipeline stops to allow user to check parameters. If
                                      everything is ok, then this parameter should be set to false, and resumed
                                      by using the -resume flag. Previous steps will be cached. If some params 
                                      are modified, then those processes affected by them and their dependant
                                      processes will be re run.
	--nanofilt_quality            The '--quality' parameter of NanoFilt. Default: 8.
	--nanofilt_maxlength          The '--maxlength' parameter of NanoFilt. Default: 1500.
	--megan_lcaAlgorithm          The '--lcaAlgorithm' parameter of daa-meganizer (MEGAN). Default: naive.
	--megan_lcaCoveragePercent    The '--lcaCoveragePercent' parameter of daa-meganizer (MEGAN). Default: 100.

	Authors: Cecilia Salazar (csalazar@pasteur.edu.uy) & Ignacio Ferres (iferres@pasteur.edu.uy)
	Maintainer: Ignacio Ferres (iferres@pasteur.edu.uy)
	
	Microbial Genomics Laboratory 
	Institut Pasteur de Montevideo (Uruguay)

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


fqs = Channel.fromPath(params.fq)
//fqs = file(params.fq)

/*
process FAKE {
	cpus 1

	input:
	file ss from fqs

	output:
	file "all.fastq" into allfq

	"""
	cp ${ss} all.fastq
	"""

}
*/

process Concatenate { 
	cpus 1

	input:
	file "*" from fqs.collect()

	output:
	file "allfq.fastq" into allfq

	shell:
	"""	
	cat * > allfq.fastq
	"""
}


process Demultiplex {
	cpus params.cpus
	
	input:
	file fq from allfq

	output:
	file "porechop_results/*" into demultiplexed1, demultiplexed2

	shell:
	"""
	porechop -t ${params.cpus} -i ${fq} -b porechop_results
	rm -f porechop_results/none.fastq
	"""
}

process Filter {
	cpus 1

	input:
	file x from demultiplexed1.flatten()

	output:
	file "Filt_${x.fileName}" into nanofilted1, nanofilted2, nanofilted3
	// file "${x.baseName}" into nanofilted2

	shell:
	"""
	cat ${x} | NanoFilt --quality ${params.nanofilt_quality} --maxlength ${params.nanofilt_maxlength} > Filt_${x.fileName}
	"""

}


process NanoPlotNoFilt {
	cpus params.cpus

	publishDir "Nanoplots", mode: "copy"
	
	input:
	file y from demultiplexed2.flatten()

	output:
	file "Nanoplot.NoFilt.${y.baseName}" optional true into nofilt
	file "count_NoFilt_${y.baseName}.txt" into counts

	shell:
	"""
	COUNT=\$(echo \$(cat ${y.fileName} | wc -l)/4 | bc)
	TWO=2
	echo \$COUNT > count_NoFilt_${y.baseName}.txt
	if [ "\$COUNT" -gt "\$TWO" ]
	then
		NanoPlot -t ${params.cpus} --fastq ${y.fileName} -o Nanoplot.NoFilt.${y.baseName}
	fi
	"""
}


process NanoPlotFilt {
	cpus params.cpus

	publishDir "Nanoplots", mode: "copy"
	
	input:
	file y from nanofilted1.flatten()

	output:
	file "Nanoplot.${y.baseName}" optional true into filt
	file "count_${y.baseName}.txt" into countsfilts

	shell:
	"""
	COUNT=\$(echo \$(cat ${y.fileName} | wc -l)/4 | bc)
	TWO=2
	echo \$COUNT > count_${y.baseName}.txt
	if [ "\$COUNT" -gt "\$TWO" ]
	then
		NanoPlot -t ${params.cpus} --fastq ${y.fileName} -o Nanoplot.${y.baseName}
	fi
	"""
}



process SummaryTable{
	cpus 1

	publishDir "Nanoplots", mode: "copy"

	input:
	file y from counts.collect()
	file x from countsfilts.collect()

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
	cpus 1

	input:
	file ff from nanofilted2

	output:
	file "${ff.baseName}.fasta" into fastas1
	
	when:
	params.stoptocheckparams == false

	shell:
	"""
	paste - - - - < ${ff} | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > ${ff.baseName}.fasta
	"""
}

silva = Channel.fromPath(params.silva_fasta)

process MakeLastDB {
	cpus params.cpus

	input:
	file fasta from silva

	output:
	file "*" into silvadb1, silvadb2

	shell:
	"""
	lastdb -cR01 -P${params.cpus} silva ${fasta}
	"""
}

/*
process LastTrain {
	cpus params.cpus

	input:
	file fa from fastas1
	file idx from silvadb1.collect()

	output:
	file "${fa}.par" into lastpars
	file "${fa.fileName}" into fastas2
		
	when:
	params.stoptocheckparams == false

	shell:
	"""
	last-train -P${params.cpus} -Q0 silva ${fa.fileName} > ${fa.baseName}.par 
	"""
}
*/

process LastAL {
	cpus params.cpus
//	publishDir "LastAlignment", mode = "copy"
	publishDir "Last_Alignment"

	input:
	//file par from lastpars
	file fa from fastas1
	file idx from silvadb2.collect()

	output:
	file "${fa.baseName}.maf" into alignment
	file "${fa.fileName}" into fastas3
		
	when:
	params.stoptocheckparams == false

	shell:
	"""
	lastal -P${params.cpus} silva ${fa.fileName} > ${fa.baseName}.maf
	
	"""
}
// lastal -P${params.cpus} -p ${par} silva ${fa.fileName} > ${fa.baseName}.maf
/*
process LastAl_DAAConverter {

	cpus params.cpus

	input:
	file fa from fastas1
	file idx from silvadb.collect()

	output:
	file "${faa.baseName}.daa" into daaconv
	file "${fa.fileName}" into fastas3

	"""
	lastal -P${params.cpus} silva ${fa.fileName} | java -jar /opt/DAA_Converter_v0.9.0.jar -top 20 -p ${params.cpus} -r ${fa.fileName} -o ${fa.baseName}.daa
	"""

}
*/

process DAAConverter{
	cpus params.cpus

	input:
	file fa from fastas3
	file maf from alignment

	output:
	file "${fa.baseName}.daa" into daaconv
	
	when:
	params.stoptocheckparams == false

	shell:
	"""
	java -jar /opt/DAA_Converter_v0.9.0.jar -top 20 -p ${params.cpus} -i ${maf.fileName} -r ${fa.fileName} -o ${fa.baseName}.daa
	"""
}


meganMap = Channel.fromPath(params.megandb_map)
meganNucl = Channel.fromPath(params.megandb_nucl)

process DAAMeganizer{
	cpus params.cpus

	input:
	file daa from daaconv.flatten()
	file mmap from meganMap.collect()
	file nucl from meganNucl.collect()

	output:
	file "${daa.fileName}" into meganized
	
	when:
	params.stoptocheckparams == false

	shell:
	"""
	daa-meganizer -i ${daa.fileName} -p ${params.cpus} -a2t ${nucl} -s2t ${mmap} --lcaAlgorithm ${params.megan_lcaAlgorithm} --lcaCoveragePercent ${params.megan_lcaCoveragePercent}
	"""
}

process ComputeComparison{
	cpus 1

	publishDir "Megan_Comparison", mode: "copy"

	input:
	file "*" from meganized.collect()

	output:
	file "comparison.megan"
	
	when:
	params.stoptocheckparams == false
	
	shell:
	"""
	xvfb-run --auto-servernum --server-num=1 /usr/local/bin/compute-comparison -i ./* -o comparison.megan
	"""
}

/*

process ExtractBiome {}

*/

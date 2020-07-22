#!/usr/bin/env nextflow

params.fq = "$baseDir/data/*.fastq"
params.outdir = "results"
params.cpus = 4
params.keepmaf = false
params.stoptocheckparams = false
params.nanofilt_quality = 8
params.nanofilt_maxlength = 1500
params.megan_lcaAlgorithm = "naive"
params.megan_lcaCoveragePercent = 100
params.help = false

nextflow.preview.dsl = 2



def helpMessage() {
    log.info """
    --------------------------
    ---> Long16S Pipeline <---
    --------------------------
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run iferres/long16S --fq 'data/*.fastq' --cpus 4 -profile local

    Note:
    SILVA and MEGAN databases are must be provided. Provide those parameters between quotes.

    Mandatory arguments:
        --fq                          Path to input data (must be surrounded with quotes).
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

include Concatenate from './modules/processes'
include Demultiplex from './modules/processes'
include Filter from './modules/processes'
include NanoPlotNoFilt from './modules/processes'
include NanoPlotFilt from './modules/processes'
include {SummaryTable} from './modules/processes'
include {Fastq2Fasta} from './modules/processes'
include {LastAL} from './modules/processes'
include {DAAConverter} from './modules/processes'
include {DAAMeganizer} from './modules/processes'
include {ComputeComparison} from './modules/processes'
include {ExtractOtuTable} from './modules/processes'

workflow {
  Channel.fromPath(params.fq)
    .set{ fqs_ch }
  Concatenate( fqs_ch.collect() )
  Demultiplex( Concatenate.out )
  Demultiplex.out
    .flatten()
    .map { file -> tuple(file.baseName, file) }
    .set{ barcode_ch }
  Filter( barcode_ch )
  NanoPlotNoFilt( barcode_ch )
  NanoPlotFilt( Filter.out )
  NanoPlotNoFilt.out.counts
    .mix( NanoPlotFilt.out.counts )
    .set{ counts_ch }
  SummaryTable( counts_ch.collect() )
  Fastq2Fasta( Filter.out )
  LastAL( Fastq2Fasta.out )
  DAAConverter( LastAL.out )
  DAAMeganizer( DAAConverter.out )
  ComputeComparison( DAAConverter.out.collect() )
  ExtractOtuTable( ComputeComparison.out )
  /*
  */
}

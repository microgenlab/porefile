#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fq = "$baseDir/data/*.fastq"
params.silvaDir = "$baseDir/silva"
params.downloadSilvaFiles = false
params.outdir = "results"
//params.cpus = 4
params.minimap2 = false
params.last = false
params.keepmaf = false
params.stoptocheckparams = false
params.nanofilt_quality = 8
params.nanofilt_maxlength = 1500
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
    nextflow run microgenlab/long16S --fq 'data/*.fastq' --minimap2 --downloadSilvaFiles

    Mandatory arguments:
        --fq                          Path to input data (must be surrounded with quotes).
        --minimap2 or --last          One, or both flags, to select workflow.

    Other:
        -profile                      Configuration profile to use. Available: standard (default), nagual, gcp.
        --outdir                      Name of the results directory. Default: "results".
      
    Boolean control:    
        --downloadSilvaFiles          Whether to download SILVA files. Requires internet connection.
        --stoptocheckparams           Whether to stop after Summary process to check parameters. Default: false.
                                      If true, then the pipeline stops to allow user to check parameters. If
                                      everything is ok, then this parameter should be set to false, and resumed
                                      by using the -resume flag. Previous steps will be cached. If some params
                                      are modified, then those processes affected by them and their dependant
                                      processes will be re run.
        --keepmaf                     
        

    Process specific parameters:
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

if ( ! (params.minimap2 || params.last) ){
  println("You must specify either --minimap2 or --last, or both.")
  System.exit(1)
}

// include modules
include {Concatenate} from './modules/processes'
include {Demultiplex} from './modules/processes'
include {Filter} from './modules/processes'
include {NanoPlotRaw} from './modules/processes'
include {NanoPlotFilt} from './modules/processes'
include {SummaryTable} from './modules/processes'

// include sub-workflows
include {DownloadSilva} from './workflows/Download'
include {LastWorkflow} from './workflows/Last'
include {Minimap2Workflow} from './workflows/Minimap2'

workflow {
  if ( params.downloadSilvaFiles ){
    DownloadSilva()
    DownloadSilva.out.fasta
      .set{ silva_fasta_ch }
    DownloadSilva.out.acctax
      .set{ silva_acctax_ch }
  } /*else {
    Channel.fromPath( "${params.silvaDir}" )
      .set{ raw_silva_ch }
  }*/
  Channel.fromPath(params.fq)
    .set{ fqs_ch }
  Concatenate( fqs_ch.collect() )
  Demultiplex( Concatenate.out )
  Demultiplex.out
    .flatten()
    .map { file -> tuple(file.baseName, file) }
    .set{ barcode_ch }
  Filter( barcode_ch )
  Filter.out
    .set{ filtered_ch }
  NanoPlotRaw( barcode_ch )
  NanoPlotFilt( Filter.out )
  NanoPlotRaw.out.counts
    .mix( NanoPlotFilt.out.counts )
    .set{ counts_ch }
  SummaryTable( counts_ch.collect() )
  if ( params.minimap2 ) {
    Minimap2Workflow( filtered_ch, silva_fasta_ch, silva_acctax_ch )
  }
  if ( params.last ) {
    LastWorkflow( filtered_ch, silva_fasta_ch, silva_acctax_ch )
  }
}

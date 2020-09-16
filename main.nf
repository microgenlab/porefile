#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fq = "$baseDir/data/*.fastq"
params.outdir = "results"
params.minimap2 = false
params.last = false
params.lasttrain = false
params.isDemultiplexed = false
params.keepmaf = false
params.stoptocheckparams = false
params.nanofilt_quality = 8
params.nanofilt_length = 1000
params.nanofilt_maxlength = 1700
params.megan_lcaAlgorithm = "naive"
params.megan_lcaCoveragePercent = 100
params.minimap2_k = 15
params.minimap2_x = "map-ont"
params.normalizeOtu = false
params.help = false



def helpMessage() {
    log.info """
    --------------------------------------------------------
    ---> porefile: a full-length 16S profiling Pipeline <---
    --------------------------------------------------------
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run microgenlab/long16S --fq 'data/*.fastq' --minimap2

    Mandatory arguments:
        --fq                          Path to input data (must be surrounded with quotes).
        --minimap2 or --last          One or both flags to select workflow.

    Other:
        --silvaFasta                  Path to SILVA_*_SSURef_NR99_tax_silva.fasta.gz file. You can provide it 
                                      either compressed (.gz) or not. If not provided, the workflow automatically
                                      adds a download step (you must have internet connection).
        --silvaFastaURL               URL to SILVA_*_SSURef_NR99_tax_silva.fasta.gz file. It will be used if you
                                      don't provide the --silvaFasta parameter (above). Default is:
                                      'https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz',
                                      change it if you want to download another SILVA version.
        --silvaAccTaxID               Path to tax_slv_ssu_*.acc_taxid.gz file. You can provide it either
                                      compressed (.gz) or not. If not provided, the workflow automatically adds
                                      a download step (you must have internet connection).
        --silvaAccTaxIDURL            URL to tax_slv_ssu_*.acc_taxid.gz file. It will be used if you don't 
                                      provide the --silvaAccTaxID parameter (above). Default is:
                                      'https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_138.acc_taxid.gz',
                                      change it if you want to download another SILVA version.
        --outdir                      Name of the results directory. Default: "results".
      
    Boolean control:    
        --stoptocheckparams           Set this option to stop the pipeline after Summary process to check parameters.
                                      If set, then the pipeline stops to allow user to check parameters. If
                                      everything is ok, then this parameter should be removed, and resumed
                                      by using the -resume flag (only one '-'!). Previous steps will be cached. 
                                      If some params are modified, then those processes affected by them and their
                                      dependant processes will be re run. Specially useful for checking NanoFilt
                                      parameters. Don't forget to allways use the -resume option!
        --keepmaf                     Use this flag if you want to copy .maf files (last workflow) to the results
                                      directory. Since generating the .maf files is time consuming and resource 
                                      intensive, and that people generally discard the 'work' directory, we decide
                                      to give this option to the users. Take into account that maf files are usually
                                      HUGE, and the file copy operation can take some time. 
        

    Process specific parameters:
        --nanofilt_quality            The '--quality' parameter of NanoFilt. Default: 8.
        --nanofilt_maxlength          The '--maxlength' parameter of NanoFilt. Default: 1500.
        --megan_lcaAlgorithm          The '--lcaAlgorithm' parameter of daa-meganizer (MEGAN). Default: naive. Possible
                                      values are: 'naive', 'weighed', or 'longReads'.
        --megan_lcaCoveragePercent    The '--lcaCoveragePercent' parameter of daa-meganizer (MEGAN). Default: 100.
        --minimap2_k                  The '-k' parameter of minimap2. Default: 15.
        --minimap2_x                  The '-x' parameter of minimap2. Default: 'map-ont'. Possible values: 'map-ont', 
                                      'asm5', 'asm10', 'asm20'.

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

if (! params.fq ){
  println("You must provide at least 1 fastq file using --fq flag.")
  System.exit(1)
}

Channel
  .fromPath(params.fq, checkIfExists: true)
  .ifEmpty { exit 1, "Non fastq files found: ${params.fq}" }
  .set{ fqs_ch }

// include modules
include {Concatenate} from './modules/processes'
include {Demultiplex} from './modules/processes'
include {Filter} from './modules/processes'
include {NanoPlotRaw} from './modules/processes'
include {NanoPlotFilt} from './modules/processes'
include {SummaryTable} from './modules/processes'
include {ExtractOtuTable} from './modules/processes'
include {ComputeComparison} from './modules/processes'

// include sub-workflows
include {SetSilva} from './workflows/Silva'
include {LastWorkflow} from './workflows/LastWorkflow'
include {Minimap2Workflow} from './workflows/Minimap2'

workflow {
  SetSilva()
    SetSilva.out.fasta
      .set{ silva_fasta_ch }
    SetSilva.out.synonyms
      .set{ silva_synonyms_ch }
  if (! params.isDemultiplexed ){
    Concatenate( fqs_ch.collect() )
    Demultiplex( Concatenate.out )
    Demultiplex.out
      .flatten()
      .map { file -> tuple(file.baseName, file) }
      .set{ barcode_ch }
  } else {
    fqs_ch
      .flatten()
      .map { file -> tuple(file.baseName, file) }
      .set{ barcode_ch }
  }
  Filter( barcode_ch )
  Filter.out
    .set{ filtered_ch }
  NanoPlotRaw( barcode_ch )
  NanoPlotFilt( Filter.out )
  NanoPlotRaw.out.counts
    .mix( NanoPlotFilt.out.counts )
    .set{ counts_ch }
  SummaryTable( counts_ch.collect() )
  Channel.empty()
    .set{ stage_to_comprare_ch }
  if ( params.minimap2 ) {
    Minimap2Workflow( filtered_ch, silva_fasta_ch, silva_synonyms_ch )
      stage_to_comprare_ch.mix( Minimap2Workflow.out )
        .set{ stage_to_comprare_ch }
      
  }
  if ( params.last || params.lasttrain  ) {
    LastWorkflow( filtered_ch, silva_fasta_ch, silva_synonyms_ch )
      stage_to_comprare_ch.mix( LastWorkflow.out )
        .set{ stage_to_comprare_ch }
  }
  ComputeComparison( stage_to_comprare_ch )
  ExtractOtuTable( ComputeComparison.out )
}

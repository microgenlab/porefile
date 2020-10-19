#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fq = "$baseDir/data/*.fastq"
params.outdir = "results"
params.minimap2 = false
params.last = false
params.lasttrain = false
params.megablast = false
params.isDemultiplexed = false
params.noNanoplot = false
params.keepmaf = false
params.stoptocheckparams = false
params.nanofilt_quality = 8
params.nanofilt_length = 1000
params.nanofilt_maxlength = 1700
params.megan_lcaAlgorithm = "naive"
params.megan_lcaCoveragePercent = 100
params.megan_topPercent = 10
params.megan_minPercentReadCover = 70
params.megan_minPercentReferenceCover = 70
params.minimap2_k = 15
params.minimap2_x = "map-ont"
params.minimap2_KM = 200
params.normalizeOtu = false
params.help = false

def sayHi(){
  log.info """
 .____   __  ____  ____  ____  __  __    ____ 
(  _ \\ /  \\(  _ \\(  __)(  __)(  )(  )  (  __)
 ) __/(  O ))   / ) _)  ) _)  )( / (_/\\ ) _) 
(__)   \\__/(__\\_)(____)(__)  (__)\\____/(____)
---------------------------------------------
... A full-length 16S profiling Pipeline ....
---------------------------------------------
  """
}

sayHi()

def helpMessage() {
    log.info """
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

if ( ! (params.minimap2 || params.last || params.lasttrain || params.megablast) ){
  println("You must specify one or more workflows to run. Implemented: --minimap2, --last, --lasttrain, and --megablast")
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
include {Fastq2Fasta} from './modules/processes'
include {MergeResults} from './modules/processes'

// include sub-workflows
include {SetSilva} from './workflows/Silva'
include {Demultiplex} from './workflows/Demultiplex'
include {QFilt} from './workflows/QFiltWorkflow'
include {QCheck} from './workflows/QCheckWorkflow'
include {LastWorkflow} from './workflows/LastWorkflow'
include {Minimap2Workflow} from './workflows/Minimap2'
include {MegaBlastWorkflow} from './workflows/MegaBlastWorkflow'

workflow {
  SetSilva()
    SetSilva.out.fasta
      .set{ silva_fasta_ch }
    SetSilva.out.acctax
      .set{ silva_acctax_ch }
  if (! params.isDemultiplexed ){
    Demultiplex( fqs_ch )
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
  QFilt( barcode_ch )
  QFilt.out.set{ filtered_scrubbed_ch }
  if (! params.noNanoplot ) {
    QCheck( barcode_ch, filtered_scrubbed_ch )
  }
  Fastq2Fasta( filtered_scrubbed_ch )
  Fastq2Fasta.out.set{ fasta_ch }
  Channel.empty()
    .set{ stage_to_comprare_ch }
  if ( params.minimap2 ) {
    Minimap2Workflow( fasta_ch, silva_fasta_ch, silva_acctax_ch )
      stage_to_comprare_ch.mix( Minimap2Workflow.out )
        .set{ stage_to_comprare_ch }
      
  }
  if ( params.last || params.lasttrain  ) {
    LastWorkflow( fasta_ch, silva_fasta_ch, silva_acctax_ch )
      stage_to_comprare_ch.mix( LastWorkflow.out )
        .set{ stage_to_comprare_ch }
  }
  if (params.megablast ){
    MegaBlastWorkflow( fasta_ch, silva_fasta_ch, silva_acctax_ch )
    stage_to_comprare_ch.mix( MegaBlastWorkflow.out )
        .set{ stage_to_comprare_ch }
  }
  MergeResults( stage_to_comprare_ch )
}

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
params.nanofilt_quality = 8
params.nanofilt_length = 1000
params.nanofilt_maxlength = 1700
params.yacrd_c = 4
params.yacrd_n = 0.4
params.megan_lcaAlgorithm = "naive"
params.megan_lcaCoveragePercent = 100
params.megan_topPercent = 10
params.megan_minPercentReadCover = 70
params.minimap2_k = 15
params.minimap2_x = "map-ont"
params.minimap2_KM = 200
params.megablast_evalue = 1e-50
params.last_E = 1e-50
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
    nextflow run microgenlab/porefile --fq 'data/*.fastq' --minimap2

    Mandatory arguments:
        --fq                          Path to input data (must be surrounded with quotes).

        One or more than one of the following available workflows:
        --minimap2                    Run Minimap2Workflow (Fast and accurate enough).
        --last                        Run LastWorkflow:Last (Not so fast but accurate).
        --lasttrain                   Run LastWorkflow:Train (Slow but more accurate).
        --megablast                   Run MegablastWorkflow (Very slow and not so accurate).

    Other:
        --silvaFasta                  Path to SILVA_132_SSURef_NR99_tax_silva.fasta.gz file. You can provide it 
                                      either compressed (.gz) or not. If not provided, the workflow automatically
                                      adds a download step (you must have internet connection).
        --silvaFastaURL               URL to SILVA_*_SSURef_NR99_tax_silva.fasta.gz file. It will be used if you
                                      don't provide the --silvaFasta parameter (above). Default is:
                                      'https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz'.
        --meganSynMap                 Path to MEGAN's SSURef_Nr99_132_tax_silva_to_NCBI_synonyms.map.gz file. You 
                                      can provide it either compressed (.gz) or not. If not provided, the workflow 
                                      automatically adds a download step (you must have internet connection).
        --meganSynMapURL              URL to SSURef_Nr99_132_tax_silva_to_NCBI_synonyms.map.gz file. It will be 
                                      used if you don't provide the --meganSynMap parameter (above). Default is:
                                      'https://software-ab.informatik.uni-tuebingen.de/download/megan6/SSURef_Nr99_132_tax_silva_to_NCBI_synonyms.map.gz'.
        --outdir                      Name of the results directory. Default: "results".
        

    Process specific parameters:
        NanoFilt parameters:
        --nanofilt_quality            The '--quality' parameter of NanoFilt. Default: 8.
        --nanofilt_length             The '--length' parameter of NanoFilt (minimum length). Default: 1000.
        --nanofilt_maxlength          The '--maxlength' parameter of NanoFilt. Default: 1500.
        Yacrd parameters:
        --yacrd_c                     The '-c' parameter of Yacrd (minimum coverage). Default: 4 .
        --yacrd_n                     The '-n' parameter of Yacrd (minimum coverage of read). Default: 0.4 .
        Minimap2 parameters:
        --minimap2_k                  The '-k' parameter of minimap2. Default: 15.
        --minimap2_x                  The '-x' parameter of minimap2. Default: 'map-ont'. Possible values: 'map-ont', 
                                      'asm5', 'asm10', 'asm20'.
        Last parameters:
        --last_E                      The '-E' parameter of lastal (e-value). Default: 1e-50.
        Megablast parameters:
        --megablast_evalue            The '-evalue' parameter of megablast. Default: 1e-50.
        Megan6 parameters:
        --megan_lcaAlgorithm          The '--lcaAlgorithm' parameter of sam2rma and blast2rma tools (Megan6). Default: 
                                      naive. Possible values are: 'naive', 'weighed', or 'longReads'.
        --megan_lcaTopPercent         The '--topPercent' parameter of sam2rma and blast2rma tools (Megan6). Default: 10.
        --megan_minPercentReadCover   The '--minPercentReadCover' parameter of sam2rma and blast2rma tools (Megan6).
                                      Default: 70.
        --megan_lcaCoveragePercent    The '--lcaCoveragePercent' parameter of sam2rma and blast2rma tools (Megan6). 
                                      Default: 100.


    Other control options:
        --isDemultiplexed             Set this flag to avoid Demultiplex sub-workflow. If set, each fastq file is 
                                      processed as a different barcode.
        --noNanoplot                  Set this flag to avoid QCheck sub-workflow. 

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

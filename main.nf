#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fq = "*.fastq"
params.outdir = "results"
params.minimap2 = false
params.last = false
params.lasttrain = false
params.megablast = false
params.isDemultiplexed = false
params.porechop_extra_end_trim = 0
params.noNanoplot = false
params.nanofilt_quality = 8
params.nanofilt_length = 1000
params.nanofilt_maxlength = 1700
params.nanofilt_headcrop = 0
params.nanofilt_tailcrop = 0
params.yacrd_c = 4
params.yacrd_n = 0.4
params.megan_lcaAlgorithm = "naive"
params.megan_lcaCoveragePercent = 100
params.megan_topPercent = 10
params.megan_minPercentReadCover = 70
params.minimap2_k = 15
params.minimap2_f = 1000
params.minimap2_x = "map-ont"
params.minimap2_KM = 200
params.megablast_evalue = 1e-50
params.last_E = 1e-50
params.help = false

params.silvaFasta = "./silvadb/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
params.silvaTaxNcbiSp = "./silvadb/Exports/taxonomy/ncbi/tax_ncbi-species_ssu_ref_nr99_138.1.txt.gz"
params.silvaTaxmap = "./silvadb/Exports/taxonomy/ncbi/taxmap_slv_ssu_ref_nr_138.1.txt.gz"

params.silvaFastaURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
params.silvaTaxNcbiSpURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/tax_ncbi-species_ssu_ref_nr99_138.1.txt.gz"
params.silvaTaxmapURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz"

params.fullSilva = false

def sayHi(){
  log.info """
 ____   __  ____  ____  ____  __  __    ____ 
(  _ \\ /  \\(  _ \\(  __)(  __)(  )(  )  (  __)
 ) __/(  O ))   / ) _)  ) _)  )( / (_/\\ ) _) 
(__)   \\__/(__\\_)(____)(__)  (__)\\____/(____)
---------------------------------------------
... A full-length 16S profiling Pipeline ....
---------------------------------------------"""
}

sayHi()

def helpMessage() {
    log.info """
    Usage:
    A typical command for running the pipeline would be as follows:

        nextflow run microgenlab/porefile --fq 'data/*.fastq' --minimap2

    Mandatory arguments:
        --fq                          Path to input data (must be surrounded with quotes).

        One or more than one of the following available workflows:
        --minimap2                    Run Minimap2Workflow (Fast and accurate enough).
        --last                        Run LastWorkflow:Last (Not so fast but accurate).
        --lasttrain                   Run LastWorkflow:Train (Slow but more accurate).
        --megablast                   Run MegablastWorkflow (Very slow and not so accurate).

    Other:
        --silvaFasta                  Path to SILVA_*_SSURef_NR99_tax_silva.fasta.gz file. You can provide it 
                                      either compressed (.gz) or not. If not provided, the workflow automatically
                                      adds a download step (you must have internet connection).
        --silvaFastaURL               URL to SILVA_*_SSURef_NR99_tax_silva.fasta.gz file. It will be used if you
                                      don't provide the --silvaFasta parameter (above). Default is:
                                      'https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz'.

        --silvaTaxNcbiSp              Path to tax_ncbi-species_ssu_ref_nr99_*.txt.gz file. You can provide it
                                      either compressed (.gz) or not. If not provided, the workflow automatically
                                      adds a download step.
        ---silvaTaxNcbiSpURL          URL to tax_ncbi-species_ssu_ref_nr99_*.txt.gz file. It will be used if you
                                      don't provide the --silvaFasta parameter (above). Default is:
                                      'https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/tax_ncbi-species_ssu_ref_nr99_138.1.txt.gz'.

        --silvaTaxmap                 Path to taxmap_slv_ssu_ref_nr_*.txt.gz file. You can provide it
                                      either compressed (.gz) or not. If not provided, the workflow automatically
                                      adds a download step.
        --silvaTaxmapURL              URL to taxmap_slv_ssu_ref_nr_*.txt.gz file. It will be used if you
                                      don't provide the --silvaFasta parameter (above). Default is:
                                      'https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz'.

        --fullSilva                   By default, porefile reduces SILVA to prokatyote SSU (16S). Use this flag
                                      to deactivate the reducing step and use the full SILVA database.

        --outdir                      Name of the results directory. Default: "results".
        

    Process specific parameters:
        Porechop parameters:
        --porechop_extra_end_trim      The '--extra_end_trim' parameter of Porechop. Default: 0.

        NanoFilt parameters:
        --nanofilt_quality            The '--quality' parameter of NanoFilt. Default: 8.
        --nanofilt_length             The '--length' parameter of NanoFilt (minimum length). Default: 1000.
        --nanofilt_maxlength          The '--maxlength' parameter of NanoFilt. Default: 1700.
        --nanofilt_headcrop           The '--headcrop' parameter of NanoFilt. Default: 0.
        --nanofilt_tailcrop           The '--tailcrop' parameter of NanoFilt. Default: 0.

        Yacrd parameters:
        --yacrd_c                     The '-c' parameter of Yacrd (minimum coverage). Default: 4 .
        --yacrd_n                     The '-n' parameter of Yacrd (minimum coverage of read). Default: 0.4 .

        Minimap2 parameters:
        --minimap2_k                  The '-k' parameter of minimap2. Default: 15.
        --minimap2_x                  The '-x' parameter of minimap2. Default: 'map-ont'. Possible values: 'map-ont', 
                                      'asm5', 'asm10', 'asm20'.
        --minimap2_f                  The '-f' parameter of minimap2. Default: 1000. Only applied in the Automap module.
        --minimap2_KM                 The '-K' parameter of minimap2, in Megabases. Default: 200.

        Last parameters:
        --last_E                      The '-E' parameter of lastal (e-value). Default: 1e-50.
        
        Megablast parameters:
        --megablast_evalue            The '-evalue' parameter of megablast. Default: 1e-50.
        
        Megan6 parameters:
        --megan_lcaAlgorithm          The '--lcaAlgorithm' parameter of sam2rma and blast2rma tools (Megan6). Default: 
                                      naive. Possible values are: 'naive', 'weighted', or 'longReads'.
        --megan_topPercent            The '--topPercent' parameter of sam2rma and blast2rma tools (Megan6). Default: 10.
        --megan_minPercentReadCover   The '--minPercentReadCover' parameter of sam2rma and blast2rma tools (Megan6).
                                      Default: 70.
        --megan_lcaCoveragePercent    The '--lcaCoveragePercent' parameter of sam2rma and blast2rma tools (Megan6). 
                                      Default: 100.


    Other control options:
        --isDemultiplexed             Set this flag to avoid Demultiplex sub-workflow. If set, each fastq file is 
                                      processed as a different barcode.
        --noNanoplot                  Set this flag to avoid QCheck sub-workflow. 

    Container options (note single dash usage!):
        -profile docker               Use docker as container engine (default).
        -profile singularity          Use singularity as container engine.
        -profile podman               Use podman as container engine.

    Help:
        --help                        Print this help and exit.

    Authors: Cecilia Salazar (csalazar@pasteur.edu.uy) & Ignacio Ferrés (iferres@pasteur.edu.uy)
    Maintainer: Ignacio Ferrés (iferres@pasteur.edu.uy)

    Microbial Genomics Laboratory
    Institut Pasteur de Montevideo (Uruguay)

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validation of parameters
def parameters_expected = [
  'help',
  'fq', 
  'outdir',
  'minimap2',
  'last',
  'lasttrain',
  'megablast',
  'isDemultiplexed', 'is-demultiplexed', // This is because  https://github.com/nextflow-io/nextflow/issues/2061
  'porechop_extra_end_trim',
  'noNanoplot', 'no-nanoplot',
  'nanofilt_quality',
  'nanofilt_length',
  'nanofilt_maxlength',
  'nanofilt_headcrop',
  'nanofilt_tailcrop',
  'yacrd_c',
  'yacrd_n',
  'megan_lcaAlgorithm', 'megan_lca-algorithm',
  'megan_lcaCoveragePercent', 'megan_lca-coverage-percent',
  'megan_topPercent', 'megan_top-percent',
  'megan_minPercentReadCover', 'megan_min-percent-read-cover',
  'minimap2_k',
  'minimap2_f',
  'minimap2_x',
  'minimap2_KM',
  'megablast_evalue',
  'last_E',
  'silvaFasta', 'silva-fasta',
  'silvaTaxNcbiSp', 'silva-tax-ncbi-sp',
  'silvaTaxmap', 'silva-taxmap',
  'silvaFastaURL', 'silva-fasta-URL',
  'silvaTaxNcbiSpURL', 'silva-tax-ncbi-sp-URL',
  'silvaTaxmapURL', 'silva-taxmap-URL',
  'fullSilva', 'full-silva'
  ] as Set

def parameter_diff = params.keySet() - parameters_expected
if (parameter_diff.size() != 0){
   exit 1, "[Pipeline error] Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
}

helloParameters()

if ( ! (params.minimap2 || params.last || params.lasttrain || params.megablast) ){
  println("You must specify one or more workflows to run. Implemented: --minimap2, --last, --lasttrain, and --megablast")
  System.exit(1)
}

if (! params.fq ){
  println("You must provide at least 1 fastq file using --fq flag.")
  System.exit(1)
}

if (!["naive", "weighted", "longReads"].contains(params.megan_lcaAlgorithm)){
  println("lcaAlgorithm not valid, must be one of 'naive', 'weighted', or 'longReads'")
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
    SetSilva.out.synonyms
      .set{ silva_synonyms_ch }
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
    Minimap2Workflow( fasta_ch, silva_fasta_ch, silva_synonyms_ch )
      stage_to_comprare_ch.mix( Minimap2Workflow.out )
        .set{ stage_to_comprare_ch }
      
  }
  if ( params.last || params.lasttrain  ) {
    LastWorkflow( fasta_ch, silva_fasta_ch, silva_synonyms_ch )
      stage_to_comprare_ch.mix( LastWorkflow.out )
        .set{ stage_to_comprare_ch }
  }
  if (params.megablast ){
    MegaBlastWorkflow( fasta_ch, silva_fasta_ch, silva_synonyms_ch )
    stage_to_comprare_ch.mix( MegaBlastWorkflow.out )
        .set{ stage_to_comprare_ch }
  }
  MergeResults( stage_to_comprare_ch )
}



def helloParameters(){

  log.info """  Nextflow-version:             $nextflow.version
  Porefile-version:             $workflow.revision
  Profile:                      $workflow.profile
  Work directory:               $workflow.workDir
  Container:                    $workflow.container
  Input directory:              $params.fq
  Output directory:             $params.outdir
  _____________________________________________""".stripIndent()
  if ( params.minimap2 ) {
    log.info """  Minimap2 aligner selected. Related parameters:
  --minimap2_k:                 $params.minimap2_k
  --minimap2_f:                 $params.minimap2_f
  --minimap2_x:                 $params.minimap2_x
  --minimap2_KM:                $params.minimap2_KM
  _____________________________________________""".stripIndent()
  }
  if ( params.last || params.lasttrain ) {
    log.info """  LAST aligner selected. Related parameters:
  --last:                       $params.last
  --lasttrain:                  $params.lasttrain
  --last_E:                     $params.last_E
  _____________________________________________""".stripIndent()
  }
  if (params.megablast){
    log.info """ Megablast aligner selected. Related parameters:
  --megablast_evalue:          $params.megablast_evalue
  _____________________________________________""".stripIndent()
  }
  log.info """ SILVAdb related parameters: """
  if (file(params.silvaFasta).exists()){
    log.info """ --silvaFasta:                 $params.silvaFasta""".stripIndent()
  }else{
    log.info """SILVAdb fasta file not provided. Download URL:
    --silvaFastaURL           $params.silvaFastaURL""".stripIndent()
  }

  if (file(params.silvaTaxNcbiSp).exists()){
    log.info """ --silvaTaxNcbiSp:             $params.silvaTaxNcbiSp""".stripIndent()
  }else{
    log.info """SILVAdb tax_ncbi-species file not provided. Download URL:
    --silvaTaxNcbiSpURL       $params.silvaTaxNcbiSpURL""".stripIndent()
  }

  if (file(params.silvaTaxmap).exists()){
    log.info """ --silvaTaxmap:                $params.silvaTaxmap""".stripIndent()
  }else{
    log.info """SILVAdb taxmap file not provided. Download URL:
    --silvaTaxmapURL          $params.silvaTaxmapURL""".stripIndent()
  }

  if (params.fullSilva){
    log.info """ Full SILVAdb selected:""".stripIndent()
  }else{
    log.info """ Reduce SILVAdb selected:""".stripIndent()
  }
  log.info """--fullSilva:                  $params.fullSilva
_____________________________________________""".stripIndent()

  log.info """Other process parameters:
Data is already demultiplexed?
--isDemultiplexed:            $params.isDemultiplexed
Porechop
--porechop_extra_end_trim:    $params.porechop_extra_end_trim
NanoFilt
--nanofilt_quality:           $params.nanofilt_quality
--nanofilt_length:            $params.nanofilt_length
--nanofilt_maxlength:         $params.nanofilt_maxlength
--nanofilt_headcrop:          $params.nanofilt_headcrop
--nanofilt_tailcrop:          $params.nanofilt_tailcrop
Yacrd
--yacrd_c:                    $params.yacrd_c
--yacrd_n:                    $params.yacrd_n
NanoPlot
--noNanoplot:                 $params.noNanoplot
_____________________________________________""".stripIndent()
}
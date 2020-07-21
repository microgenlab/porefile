#!/usr/bin/env nextflow

params.fq = "$baseDir/data/*.fastq"
params.cpus = 4
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

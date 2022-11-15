[![Nextflow](https://img.shields.io/badge/Nextflow-20.01.0-brightgreen)](https://www.nextflow.io/)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/iferres/porefile)](https://hub.docker.com/repository/docker/iferres/porefile/general)

# Porefile: a Nextflow full-length 16S profiling pipeline for ONT reads
`Porefile` is a Nextflow pipeline that wraps a bunch of third-party software to process and classify full length 16S (SSU) long reads generated using Oxford Nanopore sequencing. Users can select among a set of sub-workflows implemented (see below), or all at once.

Each sub-workflow uses different software to align ONT amplicon reads against the [SILVA](https://www.arb-silva.de/) SSU NR99 database, which is downloaded on the fly if not provided by the user.

Reads are then classified by [MEGAN6 CE](https://software-ab.informatik.uni-tuebingen.de/download/megan6/welcome.html) tools, and using a SILVA-to-NCBI accession mapping file generated at runtime if not provided by the user. 

Porefile uses SILVA SSU NR99 version 138.1 by default, which is the latest available up to this date (Nov 2022). If a new version were released, users can manually provide the new links to tell `Porefile` to download it.

![Porefile Scheme](./docs/images/scheme.png)

## Running Porefile
A typical command for running the pipeline would be as follows (don't run):
```sh
nextflow run microgenlab/porefile --fq 'path/to/*.fastq' --minimap2
```
The above command would run the `Minimap2Workflow` sub-workflow (recomended). Other available sub-workflows can be selected using these flags: `--last`, `--lasttrain`, and/or `--megablast`.

## Help
Run the following for more details about parameter tuning:
```
nextflow run microgenlab/porefile --help
```

## Dependencies
Install [Nextflow](https://www.nextflow.io/) and at least one of the following container engines: Docker, Singularity, Podman.

All workflow dependencies have been packaged into a [docker container](https://hub.docker.com/repository/docker/iferres/porefile), which is automatically downloaded when the pipeline is executed. That's it, you don't need to install any other software on your own.

Porefile has been tested with each three mencioned container technologies.

#### Dependencies included in the container

Dependencies used by the pipeline and included in the container are:
 * [Porechop](https://github.com/rrwick/Porechop) (Demultiplex)
 * [NanoFilt](https://github.com/wdecoster/nanofilt/) (Quality filtering)
 * [Yacrd](https://github.com/natir/yacrd) (Chimera removal)
 * [NanoPlot](https://github.com/wdecoster/NanoPlot) (Quality check)
 * [seqtk](https://github.com/lh3/seqtk) (fastq/fasta manipulation)
 * [last](http://last.cbrc.jp/doc/last.html) (Alignment)
 * [minimap2](https://github.com/lh3/minimap2) (Alignment)
 * [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (Alignment)
 * [r-base](https://www.r-project.org/) (Processing)
 * [MEGAN6](https://software-ab.informatik.uni-tuebingen.de/download/megan6/welcome.html) (Taxonomy assignment)

If you use Porefile, please cite them accordingly.

## Profiles
Porefile comes with a minimal set of configuration profiles. Please, refer to [Nextflow](https://www.nextflow.io/) documentation to create a configuration file for your HPC infrastructure.

#### Container engines
 * Use `-profile docker` to run the pipeline using [Docker](https://www.docker.com/). 
 * Use `-profile singularity` to run the pipeline using [Singularity](https://sylabs.io/). 
 * Use `-profile podman` to run the pipeline using [Podman](https://podman.io/). 

 #### Other configuration (for dev mostly)
  * `-profile test`: Tests the pipeline on a local machine with low resources using a toy dataset (5K ONT reads) included in the repo. Mostly used to develop on my desktop machine. Assigns at most 16Gb of RAM and 4 cpus per process. To run the test using (say) Singularity as container engine (takes about ~5min on a Intel Core i7-4790, 32Gb RAM):
  `nextflow run microgenlab/porefile --minimap2 -profile test,singularity`
  * `-profile nagual`: Configuration to use at IPMont servers.

## Usage

```
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

        --outdir                      Name of the results directory. Default: "results".
        

    Process specific parameters:
        Porechop parameters:
        --porechop_extra_end_trim      The '--extra_end_trim' parameter of Porechop. Default: 0.

        NanoFilt parameters:
        --nanofilt_quality            The '--quality' parameter of NanoFilt. Default: 8.
        --nanofilt_length             The '--length' parameter of NanoFilt (minimum length). Default: 1000.
        --nanofilt_maxlength          The '--maxlength' parameter of NanoFilt. Default: 1700.

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
        --megan_lcaTopPercent         The '--topPercent' parameter of sam2rma and blast2rma tools (Megan6). Default: 10.
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
```

## Citation
A manuscript is under preparation.


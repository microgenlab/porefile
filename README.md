[![Nextflow](https://img.shields.io/badge/Nextflow-20.01.0-brightgreen)](https://www.nextflow.io/)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/iferres/porefile)](https://hub.docker.com/repository/docker/iferres/porefile/general)

# Porefile: a Nextflow full-length 16S profiling pipeline
![Porefile Scheme](./docs/images/scheme.png)

The typical command for running the pipeline is as follows (Not run):
```
nextflow run microgenlab/porefile --fq 'path/to/*.fastq' --minimap2
```
The above command will run the Minimap2Workflow sub-workflow. Other available sub-workflows are: `--last`, `--lasttrain`, and `--megablast`.

## Dependencies
Install [Nextflow](https://www.nextflow.io/) and at least one of the following container engines: Docker, Singularity, Podman, Shifter.

All workflow dependencies have been packaged into a docker container, which is automatically downloaded and transformed into a Singularity Image File by default when the pipeline is executed. Please, refer to nextflow documentation for pipeline configuration and usage of other container engine alternatives.

## Help
Run the following for more details:
```
nextflow run microgenlab/porefile --help
```


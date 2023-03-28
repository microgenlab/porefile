# porefile 1.0.0.9000 `(xx/03/2023)`
* Refactoring. Rscripts are now in bin/ folder.
* Added seqkit to environment.yml
* Added reduceSilva process to remove eukaryota and phage sequences. Deactivate it with `--fullSilva` flag. This reduction step improves 16S classification and possibly reduces the running time a bit.
* Removed seqtk in favour of seqkit for fasta/fastq manipulation. 
* Exposed the `--headcrop` and `--tailcrop` parameters of NanoFilt.
* Synonyms mapping file now gets published in results folder. 
* Parameters are displayed at the begining of the run.
* Nanofilt quality default is now 10 (was 8 before).
* Added the Polishing workflow. Related: added megan_topPercentPolish params.
* Now read assignments get published.
* Removed Last, LastTrain, and Megablast workflows since added complexity to the code and any of them wasn't good enough as Minimap2 workflow.
* New `--help`.
* New `README.md`.
* New container, now based on micromamba instead of conda.
* New process `GetVersions` to print a file with the versions of the dependencies used (installed in the container).
* Added `--removeChimeras` to activate chimera-removing step with (Yacrd).
* Added PacBio management if user provides --minimap2_x map-pb or --minimap2_x map-hifi

# porefile 0.2.2 `14/12/2022`
* Minor fix in "nagual" config, which was using a "dev" tagged docker container. Now it uses the "latest".

# porefile 0.2.1 `24/11/2022`
* Fixed wrong parameter documentation. "--megan_topPercent" was documented as "--megan_lcaTopPercent", so this parameter was completely ignored and the default was always used instead. 
* Added parameters sanity checking. This is to avoid ignoring misspelled parameters.

# porefile 0.2.0 `15/11/2022`
* Merge s138: Now Porefiles only uses SILVA files, not needed MEGAN mapping files any more. Uses latest SILVA version (138.1).
* Fix #12: If AutoMap returns empty, this barcode is discarded.
* Added "--porechop_extra_end_trim" parameter to allow trimming n bases from the end of the amplicon in case adapters were added. 
* Minor improvements in error handling.
* Updated documentation and scheme figure.
* Updated --help

# porefile 0.1.1 `29/03/2022`
* General changes
* Updated environment.yml
* Updated Dockerfile
* Updated container to verion 0.1.1

# porefile 0.1.0 `29/10/2020`
* Initial release

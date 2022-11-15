# porefile 0.2.0
* Merge s138: Now Porefiles only uses SILVA files, not needed MEGAN mapping files any more. Uses latest SILVA version (138.1).
* Fix #12: If AutoMap returns empty, this barcode is discarded.
* Added "--porechop_extra_end_trim" parameter to allow trimming n bases from the end of the amplicon in case adapters were added. 
* Minor improvements in error handling.
* Updated documentation and scheme figure.
* Updated --help

# porefile 0.1.1
* General changes
* Updated environment.yml
* Updated Dockerfile
* Updated container to verion 0.1.1

# porefile 0.1.0
* Initial release

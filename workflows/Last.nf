nextflow.preview.dsl = 2

include MakeLastDB from '../modules/processes'
include Fastq2Fasta from '../modules/processes'
include LastAL from '../modules/processes'
include DAAConverter from '../modules/processes'
include DAAMeganizer from '../modules/processes'

workflow LastWorkflow {
    take:
        filtered_ch
        silva
    main:
        MakeLastDB( silva )
        Fastq2Fasta( filtered_ch )
        LastAL( Fastq2Fasta.out )
        DAAConverter( LastAL.out )
        DAAMeganizer( DAAConverter.out )
    emit:
        DAAMeganizer.out
}
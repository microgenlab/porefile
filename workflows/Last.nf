nextflow.preview.dsl = 2

include MakeLastDB from '../modules/processes'
include Fastq2Fasta from '../modules/processes'
include LastAL from '../modules/processes'
include DAAConverter from '../modules/processes'
include DAAMeganizer from '../modules/processes'

workflow LastWorkflow {
    take:
        filtered_ch
        silva_fasta_ch
        acctax
    main:
        MakeLastDB( silva_fasta_ch )
        Fastq2Fasta( filtered_ch )
        LastAL( Fastq2Fasta.out , MakeLastDB.out)
        DAAConverter( LastAL.out )
        DAAMeganizer( DAAConverter.out , acctax)
    emit:
        DAAMeganizer.out
}
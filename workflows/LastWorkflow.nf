nextflow.enable.dsl = 2

// include modules
include {MakeLastDB} from '../modules/processes'
include {Fastq2Fasta} from '../modules/processes'

// include sub-sub-workflows
include {Last} from './Last'
include {Train} from './Train'

workflow LastWorkflow {
    take:
        filtered_ch
        silva_fasta_ch
        acctax
    main:
        MakeLastDB( silva_fasta_ch )
        MakeLastDB.out.collect()
            .set{ lastdb_ch }
        Fastq2Fasta( filtered_ch )
        if ( params.last ){
            Last( Fastq2Fasta.out, lastdb_ch )
        }
        if ( params.lasttrain ){
            Train( Fastq2Fasta.out, lastdb_ch )
        }
        
}
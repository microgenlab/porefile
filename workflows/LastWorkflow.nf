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
        /*Fastq2Fasta( filtered_ch )
        Fastq2Fasta.out 
            .set{ fasta_ch }*/
        Channel.empty()
            .set{ to_compare_ch }
        if ( params.last ){
            Last( filtered_ch, lastdb_ch, acctax, silva_fasta_ch )
            to_compare_ch.mix( Last.out )
                .set{ to_compare_ch }
        }
        if ( params.lasttrain ){
            Train( filtered_ch, lastdb_ch, acctax, silva_fasta_ch )
            to_compare_ch.mix( Train.out )
                .set{ to_compare_ch }
        }

    emit:
        to_compare_ch
        
}
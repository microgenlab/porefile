nextflow.enable.dsl = 2

include {MakeLastDB} from '../modules/processes'
include {Fastq2Fasta} from '../modules/processes'
include {LastTrain} from '../modules/processes'
include {LastAL} from '../modules/processes'
include {DAAConverter} from '../modules/processes'
include {DAAMeganizer} from '../modules/processes'
include {ComputeComparison} from '../modules/processes'
include {ExtractOtuTable} from '../modules/processes'

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
        LastTrain( Fastq2Fasta.out, lastdb_ch )
        Fastq2Fasta.out.join( LastTrain.out )
            .set{ to_align }
        LastAL( to_align , lastdb_ch )
        DAAConverter( LastAL.out )
        DAAMeganizer( DAAConverter.out , acctax )
        Channel.value( "last" )
            .set{ workflow_ch }
        ComputeComparison( DAAMeganizer.out.collect(), workflow_ch )
        ExtractOtuTable( ComputeComparison.out, workflow_ch)
        
}
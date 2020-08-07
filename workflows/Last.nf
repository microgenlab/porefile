nextflow.enable.dsl = 2

// include modules
include {LastTrain} from '../modules/processes'
include {LastAL} from '../modules/processes'
include {DAAConverter} from '../modules/processes'
include {DAAMeganizer} from '../modules/processes'
include {ComputeComparison} from '../modules/processes'
include {ExtractOtuTable} from '../modules/processes'

workflow Last {
    take:

    main:
        Channel.value( "last" )
            .set{ workflow_ch }
        LastTrain( Fastq2Fasta.out, lastdb_ch )
        Fastq2Fasta.out.join( LastTrain.out )
            .set{ to_align }
        LastAL( to_align , lastdb_ch )
        DAAConverter( LastAL.out )
        DAAMeganizer( DAAConverter.out , acctax )
        ComputeComparison( DAAMeganizer.out.collect(), workflow_ch )
        ExtractOtuTable( ComputeComparison.out, workflow_ch)
}
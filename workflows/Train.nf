nextflow.enable.dsl = 2

// include modules
include {LastTrain} from '../modules/processes'
include {LastALPar} from '../modules/processes'
include {Blast2Rma} from '../modules/processes'
include {Rma2Info} from '../modules/processes'

workflow Train {
    take:
        filtered_ch
        lastdb_ch
        acctax
        silva_fasta_ch

    main:
        selected_wf =  "lasttrain" 
        LastTrain( filtered_ch, lastdb_ch )
        filtered_ch.join( LastTrain.out )
            .set{ to_align }
        LastALPar( to_align , lastdb_ch )
        Blast2Rma( LastALPar.out, acctax, selected_wf )
        Rma2Info( Blast2Rma.out )
        Rma2Info.out
            .groupTuple()
            .set{ train_out_ch }
    
    emit:
        train_out_ch
}

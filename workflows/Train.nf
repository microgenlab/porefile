nextflow.enable.dsl = 2

// include modules
include {LastTrain} from '../modules/processes'
include {LastALPar} from '../modules/processes'
include {DAAConverter} from '../modules/processes'
include {DAAMeganizer} from '../modules/processes'

workflow Train {
    take:
        fasta_ch
        lastdb_ch
        acctax

    main:
        selected_wf =  "lasttrain" 
        LastTrain( fasta_ch, lastdb_ch )
        fasta_ch.join( LastTrain.out )
            .set{ to_align }
        LastALPar( to_align , lastdb_ch )
        ali_ch = Channel.from(selected_wf)
            .combine( LastALPar.out )
        DAAConverter( ali_ch )
        DAAMeganizer( DAAConverter.out , acctax )
        DAAMeganizer.out
            .groupTuple()
            .set{ train_out_ch }
    
    emit:
        train_out_ch
}

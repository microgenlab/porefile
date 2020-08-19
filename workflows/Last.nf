nextflow.enable.dsl = 2

// include modules
include {LastAL} from '../modules/processes'
include {DAAConverter} from '../modules/processes'
include {DAAMeganizer} from '../modules/processes'

workflow Last {
    take:
        fasta_ch
        lastdb_ch
        acctax

    main:
        selected_wf = "last"
        LastAL( fasta_ch , lastdb_ch )
        ali_ch = Channel.from(selected_wf)
            .combine( LastAL.out )
        DAAConverter( ali_ch )
        DAAMeganizer( DAAConverter.out , acctax )
        DAAMeganizer.out
            .groupTuple()
            .set{ last_out_ch }
    
    emit:
        last_out_ch
}
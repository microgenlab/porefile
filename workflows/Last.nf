nextflow.enable.dsl = 2

// include modules
include {LastAL} from '../modules/processes'
include {Blast2Rma} from '../modules/processes'

workflow Last {
    take:
        filtered_ch
        lastdb_ch
        acctax
        silva_fasta_ch

    main:
        selected_wf = "last"
        LastAL( filtered_ch , lastdb_ch )
        Blast2Rma( LastAL.out, acctax, selected_wf )
        Blast2Rma.out
            .groupTuple()
            .set{ train_out_ch }
    
    emit:
        train_out_ch
}
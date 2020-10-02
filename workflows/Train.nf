nextflow.enable.dsl = 2

// include modules
include {LastTrain} from '../modules/processes'
include {LastALPar} from '../modules/processes'
include {Maf2Sam} from '../modules/processes'
include {Sam2Rma} from '../modules/processes'

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
        Maf2Sam( LastALPar.out, silva_fasta_ch )
        Sam2Rma( Maf2Sam.out, acctax, selected_wf )
        Sam2Rma.out
            .groupTuple()
            .set{ train_out_ch }
    
    emit:
        train_out_ch
}

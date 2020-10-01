nextflow.enable.dsl = 2

include {MakeBlastDB} from '../modules/processes'
include {MegaBlast} from '../modules/processes'
include {Blast2Rma} from '../modules/processes'

workflow MegaBlastWorkflow {
    take:
        filtered_ch
        silva_fasta_ch
        acctax

    main:
        selected_wf = "megablast"
        MakeBlastDB( silva_fasta_ch )
        MegaBlast( filtered_ch, MakeBlastDB.out )
        Blast2Rma( MegaBlast.out, acctax, selected_wf )
        Blast2Rma.out
            .groupTuple()
            .set{ to_compare_ch }

    emit:
        to_compare_ch

}

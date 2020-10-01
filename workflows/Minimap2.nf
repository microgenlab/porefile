nextflow.enable.dsl = 2

include {MakeMinimapDB} from '../modules/processes'
include {Minimap2} from '../modules/processes'
include {Sam2Rma} from '../modules/processes'
include {ComputeComparison} from '../modules/processes'

workflow Minimap2Workflow {
    take:
        filtered_ch
        silva_fasta_ch
        silva_acctax_ch

    main:
        selected_wf = "minimap2"
        MakeMinimapDB( silva_fasta_ch )
        Minimap2( filtered_ch, MakeMinimapDB.out )
        Sam2Rma( Minimap2.out, silva_acctax_ch, selected_wf )
        Sam2Rma.out
            .groupTuple()
            .set{ to_compare_ch }

    emit:
        to_compare_ch

}

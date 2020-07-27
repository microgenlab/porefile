nextflow.preview.dsl = 2

include {MakeMinimapDB} from '../modules/processes'
include {Minimap2} from '../modules/processes'
include {Sam2Rma} from '../modules/processes'

workflow Minimap2Workflow {
    take:
        filtered_ch
        silva_fasta_ch
        acctax

    main:
        MakeMinimapDB( silva_fasta_ch )
        Minimap2( filtered_ch, MakeMinimapDB.out )
        Sam2Rma( Minimap2.out, acctax )
    emit:
        Sam2Rma.out
}

nextflow.enable.dsl = 2

include {MakeMinimapDB} from '../modules/processes'
include {Minimap2} from '../modules/processes'
include {Sam2Rma} from '../modules/processes'
include {ComputeComparison} from '../modules/processes'
include {ExtractOtuTable} from '../modules/processes'

workflow Minimap2Workflow {
    take:
        filtered_ch
        silva_fasta_ch
        acctax

    main:
        MakeMinimapDB( silva_fasta_ch )
        Minimap2( filtered_ch, MakeMinimapDB.out )
        Sam2Rma( Minimap2.out, acctax )
        Channel.value( "minimap2" )
            .set{ workflow_ch }
        ComputeComparison( Sam2Rma.out.collect(), workflow_ch )
        ExtractOtuTable( ComputeComparison.out, workflow_ch)
}

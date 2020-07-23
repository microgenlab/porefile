nextflow.preview.dsl = 2

include MakeMinimapDB from '../modules/processes'
include Minimap2 from '../modules/processes'
include Sam2Rma from '../modules/processes'

workflow Minimap2Workflow {
    take:
        filtered_ch
        silva
    main:
        MakeMinimapDB( silva )
        Minimap2( filtered_ch )
        Sam2Rma( Minimap2.out )
    emit:
        Sam2Rma.out
}

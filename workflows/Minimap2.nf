nextflow.enable.dsl = 2

include {MakeMinimapDB} from '../modules/processes'
include {Minimap2} from '../modules/processes'
include {Sam2Rma} from '../modules/processes'
include {Rma2InfoR2C} from '../modules/processes'
include {MergeResults} from '../modules/processes'
include {PolishMinimap} from './PolishMinimap'

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
        Rma2InfoR2C( Sam2Rma.out )
        Rma2InfoR2C.out
            .set{ base_read_assingments_ch }

        // Publish read taxonomy assignments
        base_read_assingments_ch
            .collectFile(storeDir: "$params.outdir/Read_Assignments") {
                val, file -> 
                [ "${val}.read_info" , file ]
            }
        
        if (!params.noSpeciesPolishing){
            PolishMinimap( base_read_assingments_ch, filtered_ch, silva_fasta_ch, silva_acctax_ch )
            PolishMinimap.out
                .set{ result }
        } else {
            base_read_assingments_ch.set{ result }
        }

    emit:
        result
}

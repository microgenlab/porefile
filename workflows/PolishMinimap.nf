nextflow.enable.dsl = 2

include {Rma2InfoR2C; Rma2InfoR2C as GetReadInfo} from '../modules/processes'
include {ComputeAbundances} from '../modules/processes'
include {SubsetSilva} from '../modules/processes'
include {SubsetReads} from '../modules/processes'
include {MakeMinimapDB as MakeDB} from '../modules/processes'
include {Minimap2 as Search} from '../modules/processes'
include {Sam2Rma as Lca} from '../modules/processes' addParams(megan_topPercent: params.megan_topPercentPolish)
include {CorrectAssignment} from '../modules/processes'/*
include {MergeResults} from '../modules/processes'*/

workflow PolishMinimap {
    take:
        filtered_ch
        silva_fasta_ch
        silva_acctax_ch
        base_read_assingments_ch
        silva_ids
        read_ids

    main:
        SubsetSilva(silva_fasta_ch, silva_ids)
        MakeDB ( SubsetSilva.out )
        read_ids                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
            .flatten()
            .map {file -> tuple(file.baseName, file)}
            .set { low_ch }
        low_ch
            .join ( filtered_ch )
            .set{ low_abundance_ch }
        SubsetReads( low_abundance_ch )
        Search(SubsetReads.out, MakeDB.out)
        Lca( Search.out , silva_acctax_ch )
        GetReadInfo( Lca.out )
        base_read_assingments_ch
            .join( GetReadInfo.out )
            .set{ to_correct_ch }
        CorrectAssignment( to_correct_ch )
        CorrectAssignment.out
            .mix( base_read_assingments_ch )
            .groupTuple()
            // If there weren't corrections, return original assignment, 
            // otherwise return the correction:
            .map{ val, list -> 
                list.size() == 1 ? 
                    tuple(val, list[0]) : 
                    list[0] =~ /corrected/ ?
                        tuple(val, list[0]) :
                        tuple(val, list[1])
            }
            .set{ final_assignments_ch }


        // Publish polished read taxonomy assignments
        final_assignments_ch
            .collectFile(storeDir: "$params.outdir/Read_Assignments_Polished") {
                val, file -> 
                [ "${val}_polished.read_info" , file ]
            }
            .collect()
            .set{ all_polished_assignments }

        // Compute polished taxa counts and classification
        ComputeAbundances( all_polished_assignments, silva_acctax_ch, false ) 


    emit:
        counts = ComputeAbundances.out.counts
        taxcla = ComputeAbundances.out.taxcla

}


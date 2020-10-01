nextflow.enable.dsl = 2

// include modules
include {LastAL} from '../modules/processes'
include {Maf2Sam} from '../modules/processes'
include {Sam2Rma} from '../modules/processes'

workflow Last {
    take:
        fasta_ch
        lastdb_ch
        acctax
        silva_fasta_ch

    main:
        selected_wf = "last"
        LastAL( fasta_ch , lastdb_ch )
        Maf2Sam( LastAL.out, silva_fasta_ch )
        Sam2Rma( Maf2Sam.out, acctax )
        Channel.from(selected_wf)
            .combine( Sam2Rma.out )
            .groupTuple()
            .set{ last_out_ch }
    
    emit:
        last_out_ch
}
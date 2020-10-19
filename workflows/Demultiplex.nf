nextflow.enable.dsl = 2

include {Concatenate} from '../modules/processes'
include {Porechop} from '../modules/processes'

workflow Demultiplex {
    take:
        fqs

    main:
        Concatenate( fqs.collect() )
        Porechop( Concatenate.out )

    emit:
        Porechop.out
}

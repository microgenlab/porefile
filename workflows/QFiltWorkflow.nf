nextflow.enable.dsl = 2

include {NanoFilt} from '../modules/processes'
include {AutoMap} from '../modules/processes'
include {Yacrd} from '../modules/processes' 

workflow QFilt {
    take:
        barcode_ch

    main:
        NanoFilt( barcode_ch )
        AutoMap( NanoFilt.out )
        Yacrd( AutoMap.out )

    emit:
        Yacrd.out
}

nextflow.enable.dsl = 2

include {NanoFilt} from '../modules/processes'
include {AutoMap} from '../modules/processes'
include {Yacrd} from '../modules/processes' 

workflow QFilt {
    take:
        barcode_ch

    main:
        NanoFilt( barcode_ch )
        if (params.removeChimeras){
            AutoMap( NanoFilt.out )
            Yacrd( AutoMap.out )
            Yacrd.out
                .set{ qfilt_out }
        } else{
            NanoFilt.out
                .set{ qfilt_out }
        }

    emit:
        qfilt_out
}

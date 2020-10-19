nextflow.enable.dsl = 2

include {NanoPlotRaw} from '../modules/processes'
include {NanoPlotFilt} from '../modules/processes'
include {SummaryTable} from '../modules/processes'

workflow QCheck {
    take:
        barcode_ch
        filtered_ch

    main:
        NanoPlotRaw( barcode_ch )
        NanoPlotFilt( filtered_ch )
        NanoPlotRaw.out.counts
           .mix( NanoPlotFilt.out.counts )
            .set{ counts_ch }
        SummaryTable( counts_ch.collect() )
}

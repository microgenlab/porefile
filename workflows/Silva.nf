nextflow.enable.dsl = 2

params.silvaFasta = "./silvadb/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
params.silvaTaxNcbiSp = "./silvadb/Exports/taxonomy/ncbi/tax_ncbi-species_ssu_ref_nr99_138.1.txt.gz"
params.silvaTaxmap = "./silvadb/Exports/taxonomy/ncbi/taxmap_slv_ssu_ref_nr_138.1.txt.gz"

params.silvaFastaURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
params.silvaTaxNcbiSpURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/tax_ncbi-species_ssu_ref_nr99_138.1.txt.gz"
params.silvaTaxmapURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz"

include {downloadSilvaFasta} from '../modules/processes'
include {gunzip as gunzipFasta} from '../modules/processes'
include {downloadSilvaTaxNcbiSp} from '../modules/processes'
include {gunzip as gunzipTaxNcbiSp} from '../modules/processes'
include {downloadSilvaTaxmap} from '../modules/processes'
include {gunzip as gunzipTaxmap} from '../modules/processes'
include {generateSynonyms} from '../modules/processes'

workflow SetSilva {
    main:

    parfasta = file( params.silvaFasta )
    partaxncbisp = file( params.silvaTaxNcbiSp )
    partaxmap = file( params.silvaTaxmap )
    //paracctax = file( params.meganSynMap )

    if ( ! parfasta.exists() ){

        downloadSilvaFasta()
        downloadSilvaFasta.out
            .set{ silva_fasta_ch }

    } else{
        
        if ( parfasta.getExtension() == "gz" ){

            gunzipFasta( parfasta )
            gunzipFasta.out
                .set{ silva_fasta_ch }

        } else if ( parfasta.getExtension() == "fasta" ){

            Channel.value( parfasta ) // use value to allow recycle
                .set{ silva_fasta_ch } 

        }else{
            println("Unrecognized silvaFasta extension (not gz nor fasta).")
            System.exit(1)
        }

    }

    if ( ! partaxncbisp.exists() ){

        downloadSilvaTaxNcbiSp()
        downloadSilvaTaxNcbiSp.out
            .set{ tax_ncbi_ch }

    } else {

        if ( partaxncbisp.getExtension() == "gz" ){

            gunzipTaxNcbiSp( partaxncbisp )
            gunzipTaxNcbiSp.out
                .set{ tax_ncbi_ch }

        } else if ( partaxncbisp.getExtension() == "txt" ) {

            Channel.value( partaxncbisp )
                .set{ tax_ncbi_ch }

        } else{
            println("Unrecognized silvaTaxNcbiSp extension (not gz nor txt).")
            System.exit(1)
        }

    }

    if ( ! partaxmap.exists() ){

        downloadSilvaTaxmap()
        downloadSilvaTaxmap.out
            .set{ tax_map_ch }
        
    } else {

        if ( partaxmap.getExtension() == "gz" ){

            gunzipTaxmap( partaxmap )
            gunzipTaxmap.out
                .set{ tax_map_ch }

        } else if ( partaxncbisp.getExtension() == "txt" ) {

            Channel.value( partaxncbisp )
                .set{ tax_map_ch }

        } else{
            println("Unrecognized silvaTaxmap extension (not gz nor txt).")
            System.exit(1)
        }
    }

    generateSynonyms( tax_ncbi_ch, tax_map_ch)
    generateSynonyms.out
        .set{ silva_synonyms_ch }

/*
    if ( ! paracctax.exists() ){

        downloadMeganSynMap()
        downloadMeganSynMap.out
            .set{ to_trim }

    } else{

        if ( paracctax.getExtension() == "gz" ){
            
            gunzipMeganSynMap( paracctax )
            gunzipMeganSynMap.out
                .set{ to_trim }

        } else if ( paracctax.getExtension() == "map" ){
            
            Channel.value( paracctax ) // use value to allow recycle
                .set{ to_trim }

        }else{
            println("Unrecognized silva extension (not gz nor map).")
            System.exit(1)
        }

    }

*/
    //trimAccTaxID( to_trim )
    //trimAccTaxID.out
    //        .set{ silva_acctax_ch }
    
    emit: 
    fasta = silva_fasta_ch
    //acctax = silva_acctax_ch
    synonyms = silva_synonyms_ch
}
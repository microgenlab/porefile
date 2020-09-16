nextflow.preview.dsl = 2

params.silvaFasta = "$baseDir/silvadb/Exports/SILVA_132_SSURef_NR99_tax_silva.fasta.gz"
params.meganSynMap = "$baseDir/megan6/SSURef_Nr99_132_tax_silva_to_NCBI_synonyms.map.gz"
params.

params.silvaFastaURL = "https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz"
params.meganSynMapURL = "https://software-ab.informatik.uni-tuebingen.de/download/megan6/SSURef_Nr99_132_tax_silva_to_NCBI_synonyms.map.gz"
params.silvaTaxNcbiSp = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/tax_ncbi-species_ssu_ref_nr99_138.1.txt.gz"
params.silvaTaxmap = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/taxmap_ncbi_ssu_ref_nr99_138.1.txt.gz"


include {gunzip as gunzipFasta} from '../modules/processes'
include {gunzip as gunzipMeganSynMap} from '../modules/processes'
include {downloadFasta} from '../modules/processes'
include {downloadMeganSynMap} from '../modules/processes' 
//include {trimAccTaxID} from '../modules/processes'


workflow SetSilva {
    main:

    parfasta = file( params.silvaFasta )
    paracctax = file( params.meganSynMap )

    if ( ! parfasta.exists() ){

        downloadFasta()
        downloadFasta.out
            .set{ silva_fasta_ch }

    } else{
        
        if ( parfasta.getExtension() == "gz" ){

            gunzipFasta( parfasta )
            gunzipFasta.out
                .set{ silva_fasta_ch }

        } else if ( parfasta.getExtension() == "fasta" ){

            Channel.from( parfasta )
                .set{ silva_fasta_ch }

        }else{
            println("Unrecognized silva extension (not gz nor fasta).")
            System.exit(1)
        }

    }

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
            
            Channel.from( paracctax )
                .set{ to_trim }

        }else{
            println("Unrecognized silva extension (not gz nor map).")
            System.exit(1)
        }

    }

    //trimAccTaxID( to_trim )
    //trimAccTaxID.out
    //        .set{ silva_acctax_ch }
    
    emit: 
    fasta = silva_fasta_ch
    //acctax = silva_acctax_ch
    acctax = to_trim
}
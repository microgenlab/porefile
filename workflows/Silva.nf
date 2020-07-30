nextflow.preview.dsl = 2

params.silva-ssuref-nr99-fasta = "$baseDir/silvadb/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz"
params.silva-tax-slv-ssu-acc-taxid = "$baseDir/silvadb/Exports/taxonomy/tax_slv_ssu_138.acc_taxid.gz"

params.silvaFastaURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz"
params.silvaAccTaxIDURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_138.acc_taxid.gz"


include {downloadFasta} from '../modules/processes'
include {downloadAccTaxID} from '../modules/processes' 
include {trimAccTaxID} from '../modules/processes'


workflow SetSilva {
    main:

    parfasta = file( params.silva-ssuref-nr99-fasta )
    paracctax = file( params.silva-tax-slv-ssu-acc-taxid )
    Channel.empty().set{ stage_to_gunzip }

    if ( ! parfasta.exists() ){

        downloadFasta()
        downloadFasta.out
            .set{ silva_fasta_ch }

    } else{
        
        if ( parfasta.getExtension() == "gz" ){
             stage_to_gunzip.mix( Channel.from(parfasta) )
                .set{ stage_to_gunzip }
        } else if ( parfasta.getExtension() == "fasta" ){

        }else{
            println("Unrecognized silva extension (not gz nor fasta).")
            System.exit(1)
        }

    }

    if ( ! paracctax.exists() ){

        downloadAccTaxID()
        downloadAccTaxID.out
        trimAccTaxID( downloadAccTaxID.out )
            .set{ silva_acctax_ch }

    } else{

        if ( paracctax.getExtension() == "gz" ){
            stage_to_gunzip.mix( Channel.from(paracctax) )
                .set{ stage_to_gunzip }
        } else if ( paracctax.getExtension() == "acc_taxid" ){

        }else{
            println("Unrecognized silva extension (not gz nor acc_taxid).")
            System.exit(1)
        }

    }
    
    
    
    
    emit: 
    fasta = silva_fasta_ch
    acctax = silva_acctax_ch
}
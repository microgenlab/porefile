nextflow.preview.dsl = 2

params.silvaFastaURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz"
//params.silvaTaxMapURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.txt.gz"
params.silvaAccTaxIDURL = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_138.acc_taxid.gz"


include downloadFasta from '../modules/processes'
//include downloadTaxmap from '../modules/processes'
include downloadAccTaxID from '../modules/processes' 
include trimAccTaxID from '../modules/processes'


workflow DownloadSilva {
    main:
        downloadFasta()
        //downloadTaxmap()
        downloadAccTaxID()
        trimAccTaxID( downloadAccTaxID.out )
    
    emit: 
    fasta = downloadFasta.out
    //taxmap = downloadTaxmap.out
    acctax = trimAccTaxID.out
}
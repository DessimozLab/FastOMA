
// NXF_WRAPPER_STAGE_FILE_THRESHOLD='50000'

params.input_folder = "./in_folder/"
params.output_folder = "./out_folder/"
params.proteome_folder = params.input_folder + "/proteome"
params.proteomes = params.proteome_folder + "/*"
params.hogmap_in = params.input_folder + "/hogmap_in"

params.hogmap_folder = params.output_folder + "/hogmap"
params.splice_folder = params.output_folder + "/splice"
params.species_tree = params.input_folder + "/species_tree.nwk"
params.species_tree_checked = params.output_folder + "/species_tree_checked.nwk"

params.pickles_temp = params.output_folder + "/pickles_temp"
params.genetrees_folder = params.output_folder + "/genetrees"

params.temp_omamer_rhogs = params.output_folder +"/temp_omamer_rhogs"


process check_input{
    publishDir params.output_folder, mode: 'copy'
    input:
        path proteome_folder
        path hogmap_folder
        path species_tree
        path omamerdb
        path splice_folder
    output:
        path "species_tree_checked.nwk"
        val true
    script:
        """
        check-fastoma-input
        """
}


process omamer_run{
    publishDir params.hogmap_folder, mode: 'copy'
    input:
        path proteomes_omamerdb_inputhog
        val ready_input_check_c
    output:
        path "*.hogmap"
        val true
    script: //todo this if condition can be done as part of nextflow, so it won't submit job for cp
        """
            if [ -f ${proteomes_omamerdb_inputhog[2]}/${proteomes_omamerdb_inputhog[0]}.hogmap ]
        then
            cp ${proteomes_omamerdb_inputhog[2]}/${proteomes_omamerdb_inputhog[0]}.hogmap  ${proteomes_omamerdb_inputhog[0]}.hogmap
        else
            omamer search -n 10 --db ${proteomes_omamerdb_inputhog[1]} --query ${proteomes_omamerdb_inputhog[0]} --out ${proteomes_omamerdb_inputhog[0]}.hogmap
        fi
        """  // --nthreads 10
}



process infer_roothogs{
    publishDir = [
        path: params.temp_omamer_rhogs,
        mode:  'copy', // pattern: "temp_omamer_rhogs", saveAs: { filename -> filename.equals('temp_omamer_rhogs') ? null : filename }
        ]
        //publishDir params.temp_omamer_rhogs , mode: 'copy'
    input:
        val ready_omamer_run
        path hogmap_folder
        path proteome_folder
        path splice_folder
    output:
        path "temp_omamer_rhogs"
        path "gene_id_dic_xml.pickle"
        val true         // nextflow-io.github.io/patterns/state-dependency/
    script:
        """
        infer-roothogs  --logger-level DEBUG
        """
}


process batch_roothogs{
    input:
        val ready_infer_roothogs
        path "temp_omamer_rhogs"
    output:
        path "rhogs_rest/*", optional: true
        path "rhogs_big/*" , optional: true
        val true
    script:
        """
        batch-roothogs
        """
}

process hog_big{
    cpus 2
    publishDir params.pickles_temp
    input:
        val rhogsbig_tree_ready
    output:
        path "*.pickle"
        path "*.fa", optional: true   // msa         if write True
        path "*.nwk", optional: true  // gene trees  if write True
        val true
    script:
        """
        infer-subhogs  --input-rhog-folder ${rhogsbig_tree_ready[0]} --species-tree ${rhogsbig_tree_ready[1]} --parallel
        """
}

process hog_rest{
    publishDir params.pickles_temp
    input:
        val rhogsrest_tree_ready
    output:
        path "*.pickle"
        path "*.fa" , optional: true   // msa         if write True
        path "*.nwk" , optional: true  // gene trees  if write True
        val true
    script:
        """
        infer-subhogs  --input-rhog-folder ${rhogsrest_tree_ready[0]}  --species-tree ${rhogsrest_tree_ready[1]}
        """  // --parrallel False
}


process collect_subhogs{
    publishDir params.output_folder, mode: 'copy'
    input:
        val ready_hog_rest
        val ready_hog_big
        path "pickles_temp"   // this is the folder includes pickles_rhogs
        path "gene_id_dic_xml.pickle"
        path "temp_omamer_rhogs"
    output:
        path "output_hog.orthoxml"
        path "OrthologousGroupsFasta"
        path "OrthologousGroups.tsv"
        path "rootHOGs.tsv"
    script:
        """
        collect-subhogs
        """
}




workflow {
    proteomes = Channel.fromPath(params.proteomes,  type:'any' ,checkIfExists:true)
    proteome_folder = Channel.fromPath(params.proteome_folder)
    hogmap_folder = Channel.fromPath(params.hogmap_folder)
    splice_folder = Channel.fromPath(params.splice_folder)
    temp_omamer_rhogs = Channel.fromPath(params.temp_omamer_rhogs)
    genetrees_folder = Channel.fromPath(params.genetrees_folder)
    hogmap_in = Channel.fromPath(params.hogmap_in)

    pickles_temp =  Channel.fromPath(params.pickles_temp)
    omamerdb = Channel.fromPath(params.input_folder+"/omamerdb.h5")
    species_tree = Channel.fromPath(params.species_tree)
    species_tree_checked = Channel.fromPath(params.species_tree_checked)

    (species_tree_checked_, ready_input_check) = check_input(proteome_folder,hogmap_in,species_tree,omamerdb,splice_folder)
    ready_input_check_c = ready_input_check.collect()
    //species_tree_checked.view{"species_tree_checked ${it}"}

    proteomes_omamerdb = proteomes.combine(omamerdb)
    proteomes_omamerdb_inputhog = proteomes_omamerdb.combine(hogmap_in) // proteomes_omamerdb_inputhog.view{" rhogsbig ${it}"}
    //proteomes_omamerdb_inputhog_inputcheck =  proteomes_omamerdb_inputhog.combine(ready_input_check_c)
    (hogmap, ready_omamer_run)= omamer_run(proteomes_omamerdb_inputhog,ready_input_check_c)
    ready_omamer_run_c = ready_omamer_run.collect()

    (temp_omamer_rhogs, gene_id_dic_xml, ready_infer_roothogs) = infer_roothogs(ready_omamer_run_c, hogmap_folder, proteome_folder, splice_folder)
    ready_infer_roothogs_c = ready_infer_roothogs.collect()

    (rhogs_rest_list, rhogs_big_list, ready_batch_roothogs) = batch_roothogs(ready_infer_roothogs_c, temp_omamer_rhogs)
    ready_batch_roothogs_c = ready_batch_roothogs.collect()


    rhogsbig = rhogs_big_list.flatten()
    rhogsbig_tree =  rhogsbig.combine(species_tree_checked)
    rhogsbig_tree_ready = rhogsbig_tree.combine(ready_batch_roothogs)   //     rhogsbig_tree_ready.view{"rhogsbig_tree_ready ${it}"}
    (pickle_big_rhog, msas_out, genetrees_out, ready_hog_big) = hog_big(rhogsbig_tree_ready)

    rhogsrest = rhogs_rest_list.flatten()
    rhogsrest_tree =  rhogsrest.combine(species_tree_checked)
    rhogsrest_tree_ready = rhogsrest_tree.combine(ready_batch_roothogs_c)
    (pickle_rest_rhog,  msas_out_rest, genetrees_out_test, ready_hog_rest) = hog_rest(rhogsrest_tree_ready)

    (orthoxml_file, OrthologousGroupsFasta, OrthologousGroups_tsv, rootHOGs_tsv)  = collect_subhogs(ready_hog_rest.collect(), ready_hog_big.collect(), pickles_temp, gene_id_dic_xml, temp_omamer_rhogs)
    temp_omamer_rhogs.view{" output omamer_rhogs ${it}"}
    orthoxml_file.view{" output orthoxml file ${it}"}

}



// todo: check input files very beginning (before omamer starts) e.g all species are in the species tree. No species chars in fasta record.

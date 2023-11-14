
// NXF_WRAPPER_STAGE_FILE_THRESHOLD='50000'

params.input_folder = "./in_folder/"
params.output_folder = "./out_folder/"
params.proteome_folder = params.input_folder + "/proteome"
params.proteomes = params.proteome_folder + "/*"
params.hogmap_in = params.input_folder + "/hogmap_in"

params.hogmap_folder = params.output_folder + "/hogmap"
params.splice_folder = params.output_folder + "/splice"
params.species_tree = params.input_folder + "/species_tree.nwk"
params.pickles_temp = params.output_folder + "/pickles_temp"
params.genetrees_folder = params.output_folder + "/genetrees"


process omamer_run{
  time {4.h}
  publishDir params.hogmap_folder

  input:
  tuple path(proteome), path(omamer_db), path(precomputed_hogmap_folder)

  output:
  path "*.hogmap"
  val true

  script:
  """
    if [ -f ${precomputed_hogmap_folder}/${proteome}.hogmap ] ; then
        cp ${precomputed_hogmap_folder}/${proteome}.hogmap  ${proteome}.hogmap
    else
        omamer search --db ${omamer_db} --query ${proteome} --out ${proteome}.hogmap
    fi
  """
}


process infer_roothogs{
  input:
    val ready_omamer_run
    path hogmap_folder
    path proteome_folder
    path splice_folder
  output:
    path "omamer_rhogs/*"
    path "gene_id_dic_xml.pickle"
  script:
    """
       infer-roothogs  --logger-level DEBUG
    """
}


process batch_roothogs{
  input:
    path rhogs, stageAs: "omamer_rhogs/*"
  output:
    path "rhogs_rest/*", optional: true
    path "rhogs_big/*" , optional: true
  script:
    """
        batch-roothogs --input-roothogs omamer_rhogs/ --out-big rhogs_big --out-rest rhogs_rest -vvv
    """
}

process hog_big{
  cpus  2
  time {20.h}     // for very big rhog it might need more, or you could re-run and add `-resume`
  input:
    val rhogsbig_tree_ready
  output:
    path "*.pickle"
    path "*.fa", optional: true   // msa         if write True
    path "*.nwk", optional: true  // gene trees  if write True
    val true
  script:
    """
        infer-subhogs  --input-rhog-folder ${rhogsbig_tree_ready[0]} --species-tree ${rhogsbig_tree_ready[1]} --parallel --fragment-detection --low-so-detection
    """
}

process hog_rest{
  input:
    each rhogsrest
    path species_tree
  output:
    path "pickle_hogs"
    path "msa/*.fa" , optional: true   // msa         if write True
    path "gene_trees/*.nwk" , optional: true  // gene trees  if write True
    val true
  script:
    """
        infer-subhogs  --input-rhog-folder ${rhogsrest}  \
                       --species-tree ${species_tree} \
                       --fragment-detection \
                       --low-so-detection
                       #--out pickle_hogs
    """
}


process collect_subhogs{
  publishDir params.output_folder, mode: 'copy'
  input:
    path pickles, stageAs: "pickle_folders/?"
    path "gene_id_dic_xml.pickle"
    path rhogs, stageAs: "omamer_rhogs/*"
  output:
    path "output_hog.orthoxml"
    path "OrthologousGroupsFasta"
    path "OrthologousGroups.tsv"
    path "rootHOGs.tsv"
  script:
    """
        collect-subhogs --pickle-folder pickle_folders/ \
                        --roothogs-folder omamer_rhogs/ \
                        --gene-id-pickle-file gene_id_dic_xml.pickle \
                        --out output_hog.orthoxml \
                        -vv
    """
}

workflow {
    proteomes = Channel.fromPath(params.proteomes, type:'any', checkIfExists:true)
    proteome_folder = Channel.fromPath(params.proteome_folder)
    hogmap_folder = Channel.fromPath(params.hogmap_folder)
    splice_folder = Channel.fromPath(params.splice_folder)
    species_tree = Channel.fromPath(params.species_tree)

    genetrees_folder = Channel.fromPath(params.genetrees_folder)
    hogmap_in = Channel.fromPath(params.hogmap_in, type:'dir')

    pickles_temp =  Channel.fromPath(params.pickles_temp)
    omamerdb = Channel.fromPath(params.omamer_db)
    proteomes_omamerdb = proteomes.combine(omamerdb)
    proteomes_omamerdb_inputhog = proteomes_omamerdb.combine(hogmap_in)

    (hogmap, ready_omamer_run) = omamer_run(proteomes_omamerdb_inputhog)
    ready_omamer_run_c = ready_omamer_run.collect()

    (omamer_rhogs, gene_id_dic_xml, ready_infer_roothogs) = infer_roothogs(ready_omamer_run_c, hogmap_folder, proteome_folder, splice_folder)

    (rhogs_rest_batches, rhogs_big_list) = batch_roothogs(omamer_rhogs)
    rhogs_rest_batches.flatten().view{ "batch $it" }


    //rhogsbig = rhogs_big_list.flatten()
    //rhogsbig_tree =  rhogsbig.combine(species_tree)
    //rhogsbig_tree_ready = rhogsbig_tree.combine(ready_batch_roothogs)   //     rhogsbig_tree_ready.view{"rhogsbig_tree_ready ${it}"}
    //(pickle_big_rhog, msas_out, genetrees_out, ready_hog_big) = hog_big(rhogsbig_tree)

    //rhogsrest_tree =  rhogs_rest_list.combine(species_tree)
    //rhogsrest_tree_ready = rhogsrest_tree.combine(ready_batch_roothogs_c)
    (pickle_rest_rhog,  msas_out_rest, genetrees_out_test, ready_hog_rest) = hog_rest(rhogs_rest_batches.flatten(), species_tree)
    pickle_rest_rhog.view()

    (orthoxml_file, OrthologousGroupsFasta, OrthologousGroups_tsv, rootHOGs_tsv)  = collect_subhogs(pickle_rest_rhog.collect(), gene_id_dic_xml, omamer_rhogs)
    orthoxml_file.view{" output orthoxml file ${it}"}

}

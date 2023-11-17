
// NXF_WRAPPER_STAGE_FILE_THRESHOLD='50000'

params.input_folder = "./in_folder/"
params.output_folder = "./out_folder/"
params.proteome_folder = params.input_folder + "/proteome"
params.hogmap_in = params.input_folder + "/hogmap_in"
params.splice_folder = params.input_folder + "/splice"
params.species_tree = params.input_folder + "/species_tree.nwk"




// output subfolder definition
params.genetrees_folder = params.output_folder + "/genetrees"
params.hogmap_folder = params.output_folder + "/hogmap"

params.temp_output = params.output_folder +"/temp_output" //"/temp_omamer_rhogs"



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
    script:
        """
        check-fastoma-input --proteomes ${proteome_folder} \
                            --species-tree ${species_tree} \
                            --out-tree species_tree_checked.nwk \
                            --splice ${splice_folder} \
                            --hogmap ${hogmap_folder} \
                            --omamer_db ${omamerdb} \
                            -vv
        """
}


process omamer_run{
  time {4.h}
  publishDir params.hogmap_folder, mode: 'copy'

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
  publishDir = [
    path: params.temp_output,
    mode: 'copy', // pattern: "temp_output", saveAs: { filename -> filename.equals('temp_omamer_rhogs') ? null : filename }
    ]

  input:
    val ready_omamer_run
    path hogmap_folder
    path proteome_folder
    path splice_folder
  output:
    path "omamer_rhogs/*"
    path "gene_id_dic_xml.pickle"
    path "selected_isoforms" , optional: true
  script:
    """
       infer-roothogs  --proteomes ${proteome_folder} \
                       --hogmap ${hogmap_folder} \
                       --splice ${splice_folder} \
                       --out-rhog-folder "omamer_rhogs/" \
                       -vv
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
    each rhogsbig
    path species_tree
  output:
    path "pickle_hogs"
    path "msa/*.fa" , optional: true          // msa         if write True
    path "gene_trees/*.nwk" , optional: true  // gene trees  if write True
  script:
    """
        infer-subhogs  --input-rhog-folder ${rhogsbig}  \
                       --species-tree ${species_tree} \
                       --fragment-detection \
                       --low-so-detection \
                       --parallel
    """
}

process hog_rest{
  input:
    each rhogsrest
    path species_tree
  output:
    path "pickle_hogs"
    path "msa/*.fa" , optional: true          // msa         if write True
    path "gene_trees/*.nwk" , optional: true  // gene trees  if write True
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
                        --marker-groups-fasta OrthologousGroups.tsv \
                        --roothog-tsv rootHOGs.tsv \
                        -vv
    """
}

workflow {
    proteome_folder = Channel.fromPath(params.proteome_folder, type: "dir", checkIfExists:true).first()
    proteomes = Channel.fromPath(params.proteome_folder + "/*", type:'any', checkIfExists:true)
    species_tree = Channel.fromPath(params.species_tree, type: "file", checkIfExists:true).first()

    hogmap_folder = Channel.fromPath(params.hogmap_folder, type: "dir")
    splice_folder = Channel.fromPath(params.splice_folder, type: "dir")

    genetrees_folder = Channel.fromPath(params.genetrees_folder, type: 'dir')
    hogmap_in = Channel.fromPath(params.hogmap_in, type:'dir')

    omamerdb = Channel.fromPath(params.omamer_db)
    proteomes_omamerdb = proteomes.combine(omamerdb)
    proteomes_omamerdb_inputhog = proteomes_omamerdb.combine(hogmap_in)

    (species_tree_checked_, ready_input_check) = check_input(proteome_folder, hogmap_in, species_tree, omamerdb, splice_folder)
    (hogmap, ready_omamer_run) = omamer_run(proteomes_omamerdb_inputhog)
    ready_omamer_run_c = ready_omamer_run.collect()

    (omamer_rhogs, gene_id_dic_xml, ready_infer_roothogs) = infer_roothogs(ready_omamer_run_c, hogmap_folder, proteome_folder, splice_folder)

    (rhogs_rest_batches, rhogs_big_batches) = batch_roothogs(omamer_rhogs)

    (pickle_big_rhog, msa_out_big, genetrees_out_rest) = hog_big(rhogs_big_batches.flatten(), species_tree)
    (pickle_rest_rhog,  msas_out_rest, genetrees_out_test) = hog_rest(rhogs_rest_batches.flatten(), species_tree)
    channel.empty().concat(pickle_big_rhog, pickle_rest_rhog).set{ all_rhog_pickle }

    (orthoxml_file, OrthologousGroupsFasta, OrthologousGroups_tsv, rootHOGs_tsv)  = collect_subhogs(all_rhog_pickle.collect(), gene_id_dic_xml, omamer_rhogs)
    orthoxml_file.view{" output orthoxml file ${it}"}

}

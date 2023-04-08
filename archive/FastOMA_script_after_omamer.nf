

params.input_folder = "./in_folder/"
params.output_folder = "./out_folder/"
params.proteome_folder = params.input_folder + "/proteome"
params.proteomes = params.proteome_folder + "/*"

//params.hogmap_folder = params.output_folder + "/hogmap"
params.hogmap_folder = "./hogmap"
params.rhogs_folder = params.output_folder + "/rhogs_all"
params.species_tree = params.input_folder + "/species_tree.nwk"
params.pickles_rhogs_folder = params.output_folder + "/pickle_rhogs"


process omamer_run{
  time {4.h}
  memory {50.GB}
  cpus  10
  publishDir params.hogmap_folder
  input:
  path proteomes_omamerdb
  output:
  path "*.hogmap"
  val true      // ready_omamer_run
  script:
  //   omamer search --db ${proteomes_omamerdb[1]} --query ${proteomes_omamerdb[0]} --nthreads 1  --out ${proteomes_omamerdb[0]}.hogmap
  // cp /work/FAC/FBM/DBC/cdessim2/default/smajidi1/qfo_hogmap/${proteomes_omamerdb[0]}.hogmap .
  """
  omamer search --db ${proteomes_omamerdb[1]} --query ${proteomes_omamerdb[0]} --nthreads 10  --out ${proteomes_omamerdb[0]}.hogmap
  """

}


process infer_roothogs{
  publishDir  params.rhogs_folder // "${params.output_folder}/rhogs_all"
  input:
  val ready_omamer_run
  path hogmap_folder
  path proteome_folder
  output:
  path "*.fa"
  path "gene_id_dic_xml.pickle"
  val true     // ready_infer_roothogs   nextflow-io.github.io/patterns/state-dependency/
  script:
  """
   infer-roothogs  --logger-level DEBUG
  """
}

process batch_roothogs{
  publishDir params.output_folder
  input:
  val ready_infer_roothogs
  path rhogs_folder //"${params.output_folder}/rhogs_all"
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
  cpus  8
  time {10.h}    // for very big rhog it might need more, or you could re-run and add `-resume`
  memory {80.GB}
  publishDir params.pickles_rhogs_folder
  input:
  // val ready_batch_roothogs
  // path rhogsbig_tree // = rhogsbig.combine(species_tree)
  // rhogs_big_i  //"$rhogs_big/*.fa"
  // path "species_tree.nwk"
  val rhogsbig_tree_ready
  output:
  path "*.pickle"
  val true
  // path "pi_big_subhog/*"
  // pi_big rhogs_big
  // params.species_tree

  script:
  """
  infer-subhogs  --input-rhog-folder ${rhogsbig_tree_ready[0]} --parrallel True  --species-tree ${rhogsbig_tree_ready[1]}
  """
}


process hog_rest{
  publishDir params.pickles_rhogs_folder

  input:
  // val ready_batch_roothogs
  //path rhogsrest_tree // = rhogsrest.combine(species_tree)
  val rhogsrest_tree_ready

  output:
  path "*.pickle"

  val true
  script:
  """
  infer-subhogs  --input-rhog-folder ${rhogsrest_tree_ready[0]} --parrallel False --species-tree ${rhogsrest_tree_ready[1]}
  """
}

process collect_subhogs{

  memory {50.GB}
  publishDir params.output_folder, mode: 'copy'
  input:
  val ready_hog_rest
  val ready_hog_big
  // path pickle_rhogs   // this is for depenedcy
  path "pickle_rhogs" // this is the folder includes pickles_rhogs
  path "gene_id_dic_xml.pickle"

  output:
  path "output_hog_.orthoxml"

  script:
  """
   collect-subhogs
  """
}



workflow {
    proteomes = Channel.fromPath(params.proteomes,  type:'any' ,checkIfExists:true)
    proteome_folder = Channel.fromPath(params.proteome_folder)
    hogmap_folder = Channel.fromPath(params.hogmap_folder)
    // params.hogmap_folder = params.output_folder + "/hogmap"


    rhogs_folder = Channel.fromPath(params.rhogs_folder)
    pickles_rhogs_folder =  Channel.fromPath(params.pickles_rhogs_folder)
    // omamerdb = Channel.fromPath(params.input_folder+"/omamerdb.h5")
    // proteomes.view{"prot ${it}"}
    // proteomes_omamerdb = proteomes.combine(omamerdb)
    // proteomes_omamerdb.view{"proteomes_omamerdb ${it}"}
    // (hogmap, ready_omamer_run)= omamer_run(proteomes_omamerdb)
    // ready_omamer_run_c = ready_omamer_run.collect()
    ready_omamer_run_c = true
    // hogmaps.view{"hogmap ${it}"}

    // proteome_folder.view{"proteome_folder ${it} "}
    // (rhogs, gene_id_dic_xml) = infer_roothogs(hogmaps, hogmap_folder, proteome_folder)
    (rhogs, gene_id_dic_xml, ready_infer_roothogs) = infer_roothogs(ready_omamer_run_c, hogmap_folder, proteome_folder)
    // rhogs.view{"rhogs ${it}"}
    // rhogs_folder.view{"rhogs_folder xx ${it}"}

    ready_infer_roothogs_c = ready_infer_roothogs.collect()
    (rhogs_rest_list, rhogs_big_list, ready_batch_roothogs) = batch_roothogs(ready_infer_roothogs_c, rhogs_folder)
    ready_batch_roothogs_c = ready_batch_roothogs.collect()

    ready_batch_roothogs_c.view{" ready_batch_roothogs_c 44  ${it}"}

    species_tree = Channel.fromPath(params.species_tree)
    rhogsbig = rhogs_big_list.flatten()
    // rhogsbig.view{" rhogsbig ${it}"}
    rhogsbig_tree =  rhogsbig.combine(species_tree)
    rhogsbig_tree_ready = rhogsbig_tree.combine(ready_batch_roothogs)
    rhogsbig_tree_ready.view{"rhogsbig_tree_ready ${it}"}
    (pickle_big_rhog, ready_hog_big) = hog_big(rhogsbig_tree_ready)

    rhogsrest = rhogs_rest_list.flatten()
//     rhogsrest.view{" rhogs rest ${it}"}
    rhogsrest_tree =  rhogsrest.combine(species_tree)




    rhogsrest_tree_ready = rhogsrest_tree.combine(ready_batch_roothogs_c)
//     rhogsrest_tree_ready.view{"rhogsrest_tree_ready ${it}"}

    (pickle_rest_rhog, ready_hog_rest) = hog_rest(rhogsrest_tree_ready)

//     pickle_rest_rhog.flatten().view{" pickle_rest_rhog rest ${it}"}
//     pickle_big_rhog.flatten().view{" pickle_big_rhog rest ${it}"}
    prb = pickle_big_rhog.collect()
    prr = pickle_rest_rhog.collect()
    all_pickles = prb.mix(prr)
//     gene_id_dic_xml = Channel.fromPath("gene_id_dic_xml.pickle")
    pickle_rhogs_folder = Channel.fromPath(params.output_folder+"/pickle_rhogs")
//     orthoxml_file = collect_subhogs(all_pickles.collect(), pickle_rhogs_folder, gene_id_dic_xml)

    orthoxml_file = collect_subhogs(ready_hog_rest.collect(), ready_hog_big.collect(), pickles_rhogs_folder, gene_id_dic_xml)
    orthoxml_file.view{" output orthoxml file ${it}"}



}

// memory {12.GB * (2*task.attempt - 1)}
//    time {24.hour}
//    errorStrategy {
//      task.exitStatus in [1,99,143,137,104,134,139,145,140] ? ‘retry’ : ‘terminate’
//    }
//    maxRetries 4

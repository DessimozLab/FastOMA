
// nextflow /work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/gethog3/gethog3_script.nf  --input_folder test2/   --output_folder t5/  -resume

params.input_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/test/testdata/working_folder/"
params.output_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/test/testdata/"
params.proteome_folder = params.input_folder + "/proteome"
params.proteomes = params.proteome_folder + "/*"

params.hogmap_folder = params.output_folder + "/hogmap"
params.rhogs_folder = params.output_folder + "/rhogs_all"

// params.rhogs_big_folder = params.input_folder + "rhogs_big"

process omamer_run{

  memory {50.GB}

  publishDir params.hogmap_folder
  input:
  path proteomes_omamerdb
  output:
  path "*.hogmap"
  script:
  """
  omamer search --db ${proteomes_omamerdb[1]} --query ${proteomes_omamerdb[0]} --nthreads 1  --out ${proteomes_omamerdb[0]}.hogmap
  """
}


process infer_roothogs{

  publishDir  params.rhogs_folder // "${params.output_folder}/rhogs_all"
  input:
  path hogmaps
  path hogmap_folder
  path proteome_folder
  output:
  path "*.fa"
  path "gene_id_dic_xml.pickle"

  script:
  """
   infer-roothogs  --logger-level DEBUG
  """
}

process batch_roothogs{
  publishDir params.output_folder
  input:
  path rhogs
  path rhogs_folder //"${params.output_folder}/rhogs_all"

  output:
  path "rhogs_rest/*", optional: true
  path "rhogs_big/*" , optional: true
  script:
  """
   batch-roothogs
  """
}

process hog_big{

  cpus  8
  time {10.h}    // for very big rhog it might need more, or you could re-run and add `-resume`
  memory {50.GB}

  publishDir params.output_folder+"/pickle_rhogs"

  input:
  path rhogsbig_tree // = rhogsbig.combine(species_tree)
  // rhogs_big_i  //"$rhogs_big/*.fa"
  // path "species_tree.nwk"
  output:
  path "*.pickle"
  // path "pi_big_subhog/*"
  // pi_big rhogs_big
  script:
  """
  infer-subhogs  --input-rhog-folder ${rhogsbig_tree[0]} --parrallel True
  """
}


process hog_rest{
  publishDir params.output_folder+"/pickle_rhogs"

  input:
  path rhogsrest_tree // = rhogsrest.combine(species_tree)

  output:
  path "*.pickle"
  script:
  """
  infer-subhogs  --input-rhog-folder ${rhogsrest_tree[0]} --parrallel False
  """
}




process collect_subhogs{
  publishDir params.output_folder, mode: 'copy'
  input:
  path pickle_rhogs   // this is for depenedcy
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
    rhogs_folder = Channel.fromPath(params.rhogs_folder)

    // hogmap_folder.view{"hogmap_folder ${it}"}  // using provided argument

    omamerdb = Channel.fromPath(params.input_folder+"/omamerdb.h5")
    // proteomes.view{"prot ${it}"}
    proteomes_omamerdb = proteomes.combine(omamerdb)
    // proteomes_omamerdb.view{"proteomes_omamerdb ${it}"}

    hogmap = omamer_run(proteomes_omamerdb)
    hogmaps = hogmap.collect()
//     hogmaps.view{"hogmap ${it}"}

    // proteome_folder.view{"proteome_folder ${it} "}
    (rhogs, gene_id_dic_xml) = infer_roothogs(hogmaps, hogmap_folder, proteome_folder)
    rhogs.view{"rhogs ${it}"}
    // rhogs_folder.view{"rhogs_folder xx ${it}"}

    (rhogs_rest_list, rhogs_big_list) = batch_roothogs(rhogs, rhogs_folder)
    rhogs_rest_list.view{"rhogs_rest_list ${it}"}

    rhogsrest=rhogs_rest_list.flatten()
    rhogsrest.view{" rhogs rest ${it}"}

    rhogsbig = rhogs_big_list.flatten()
    rhogsbig.view{" rhogs big ${it}"}

    species_tree = Channel.fromPath(params.input_folder + "/species_tree.nwk")
    rhogsbig_tree =  rhogsbig.combine(species_tree)
    rhogsbig_tree.view{"rhogsbig_tree ${it}"}

    rhogsrest_tree =  rhogsrest.combine(species_tree)
    rhogsrest_tree.view{"rhogsrest_tree ${it}"}

    pickle_big_rhog = hog_big(rhogsbig_tree)
    pickle_rest_rhog = hog_rest(rhogsrest_tree)

    pickle_rest_rhog.flatten().view{" pickle_rest_rhog rest ${it}"}
    pickle_big_rhog.flatten().view{" pickle_big_rhog rest ${it}"}

    prb = pickle_big_rhog.collect()
    prr = pickle_rest_rhog.collect()
    all_pickles = prb.mix(prr)

    // gene_id_dic_xml = Channel.fromPath("gene_id_dic_xml.pickle")

    pickle_rhogs_folder = Channel.fromPath(params.output_folder+"/pickle_rhogs")
    orthoxml_file = collect_subhogs(all_pickles.collect(), pickle_rhogs_folder, gene_id_dic_xml)
    orthoxml_file.view{" output orthoxml file ${it}"}

}



includeConfig "nextflow_slurm.config"

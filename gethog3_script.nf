
params.input_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/test/testdata/working_folder/"
params.proteome_folder = params.input_folder + "proteome"
params.hogmap_folder = params.input_folder + "hogmap"
params.rhogs_folder = params.input_folder + "rhogs_all"
params.proteomes = params.proteome_folder+"/*"
// params.rhogs_big_folder = params.input_folder + "rhogs_big"

process omamer_run{
  publishDir params.hogmap_folder
  input:
  path proteomes_omamerdb
  output:
  path "*.hogmap"
  script:
  """
  omamer search --db ${proteomes_omamerdb[1]} --query ${proteomes_omamerdb[0]} --nthreads 2  --out ${proteomes_omamerdb[0]}.hogmap
  """
}


process infer_roothogs{
  publishDir "rhogs_all"
  input:
  path hogmaps
  path hogmap_folder
  path proteome_folder
  output:
  path "*.fa"
  script:
  """
   infer-roothogs  --logger-level DEBUG
  """
}

process batch_roothogs{

  input:
  path rhogs
  path "rhogs_all"

  output:
  path "rhogs_rest/*", optional: true
  path "rhogs_big/*" , optional: true
  script:
  """
   batch-roothogs
  """
}

process hog_big{
  publishDir "pickle_rhogs"
  input:
  path rhogs_big_i  //"$rhogs_big/*.fa"

  output:
  path "*.pickle"
  // path "pi_big_subhog/*"
  // pi_big rhogs_big
  script:
  """
  infer-subhogs  --input-rhog-folder $rhogs_big_i --parrallel False
  """
}


workflow {

    proteomes = Channel.fromPath(params.proteomes,  type:'any' ,checkIfExists:true)
    proteome_folder = Channel.fromPath(params.proteome_folder)
    hogmap_folder = Channel.fromPath(params.hogmap_folder)
    rhogs_folder = Channel.fromPath(params.rhogs_folder)

    omamerdb = Channel.fromPath("omamerdb.h5")
    proteomes.view{"prot ${it}"}
    proteomes_omamerdb = proteomes.combine(omamerdb)
    proteomes_omamerdb.view{"proteomes_omamerdb ${it}"}

    hogmap = omamer_run(proteomes_omamerdb)
    hogmaps = hogmap.collect()
    hogmaps.view{"hogmap ${it}"}

    proteome_folder.view{"proteome_folder ${it} "}
    rhogs = infer_roothogs(hogmaps, hogmap_folder, proteome_folder)
    rhogs.view{"rhogs ${it}"}

    (rhogs_rest_list, rhogs_big_list) = batch_roothogs(rhogs, rhogs_folder)
    // rhogs_rest_list.view{"rhogs_rest_list ${it}"}

    rhogs_rest=rhogs_rest_list.flatten()
    rhogs_rest.view{" rhogs rest ${it}"}

    rhogs_big=rhogs_big_list.flatten()
    rhogs_big.view{" rhogs big ${it}"}

    hog_big(rhogs_big)


}



//
//

// process hog_rest{
//   publishDir "${params.outputdir}/pickle_rhogs/"
//   input:
//   path rhogs_rest_i   //"$rhogs_big/*.fa"
//   val gethog3
//   output:
//   path "*.pickle"
//
//   script:
//   """
//    python ${gethog3}//infer_folder.py  $rhogs_rest_i False pi_rest rhogs_rest
//   """
// }
// process collect_orthoxml{
//   publishDir "${params.outputdir}"
//   input:
//   path pickle_rhogs
//   // path gene_id_dic_xml
//   val gethog3
//   output:
//   path "output_hog_.orthoxml"
//
//   script:
//   """
//    python ${gethog3}/collect_orthoxml.py
//   """
// }
//
//     hogmap = omamer_run(proteomes, omamer_db, num_threads_omamer,omamer,outputdir)
//     rhogs = inferrhog(hogmap.collect(), gethog3)
//     rhogs.flatten().view{"rhogs ${it}"}

//     pickle_rest_rhog = hog_rest(rhogs_rest, gethog3)
//     pickle_rest_rhog.flatten().view{" pickle_rest_rhog rest ${it}"}
//
//     pickle_big_rhog = hog_big(rhogs_big, gethog3)
//     pickle_big_rhog.flatten().view{" pickle_big_rhog rest ${it}"}
//
//
//     prb = pickle_big_rhog.collect()
//     prr = pickle_rest_rhog.collect()
//     all_pickles = prb.mix(prr)
//     ortho = collect_orthoxml(all_pickles.collect(), gethog3)
//     ortho.view{" output orthoxml file ${it}"}



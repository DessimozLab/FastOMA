params.proteomes = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/proteome/*"
params.omamer_db= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/archive/Primates.h5" // LUCA.h5"
params.species_tree= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/tree.nwk"
params.num_threads_omamer= 2
params.gethog3= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/gethog3/gethog3/"
params.omamer = "omamer" // if installed, otherwise address to the executable
params.outputdir = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/"


process omamer_run{
  publishDir "${params.outputdir}/hogmap/"
  input:
  path proteomes
  val omamer_db
  val num_threads_omamer
  val omamer
  val outputdir
  output:
  path "*.hogmap"
  script:
  """
   ${omamer} search --db $omamer_db  --query $proteomes --nthreads $num_threads_omamer  --out ${proteomes}.hogmap
  """
}

process inferrhog{
  publishDir "${params.outputdir}/rhogs_all/"
  input:
  path hogmap
  // path proteomes
  val gethog3
  output:
  path "*.fa"
  // path "gene_id_dic_xml.pickle"
  script:
  """
   python ${gethog3}/infer_rhog.py ./
  """
}

process rhog_distributor{
  publishDir "${params.outputdir}/"

  input:
  path rhogs
  val gethog3
  output:
  path "rhogs_rest/*", optional: true
  path "rhogs_big/*" , optional: true
  script:
  """
   python ${gethog3}/rhog_distributor.py
  """
}


process hog_big{
  publishDir "${params.outputdir}/pickle_rhogs/"
  input:
  path rhogs_big_i //"$rhogs_big/*.fa"
  val gethog3
  output:
  path "*.pickle"
  //path "pi_big_subhog/*"

  script:
  """
   python ${gethog3}//infer_folder.py  $rhogs_big_i False pi_big rhogs_big
  """
}


process hog_rest{
  publishDir "${params.outputdir}/pickle_rhogs/"
  input:
  path rhogs_rest_i //"$rhogs_big/*.fa"
  val gethog3
  output:
  path "*.pickle"

  script:
  """
   python ${gethog3}//infer_folder.py  $rhogs_rest_i False pi_rest rhogs_rest
  """
}


process collect_orthoxml{
  publishDir "${params.outputdir}"
  input:
  path pickle_rhogs
  // path gene_id_dic_xml
  val gethog3
  output:
  path "hog__3.orthoxml"

  script:
  """
   python ${gethog3}/collect_orthoxml.py
  """
}




workflow {

    proteomes = Channel.fromPath(params.proteomes,  type:'any' ,checkIfExists:true)
    omamer_db = Channel.value(params.omamer_db)
    num_threads_omamer = Channel.value(params.num_threads_omamer)
    species_tree = Channel.value(params.species_tree)
    gethog3 = Channel.value(params.gethog3)
    omamer = Channel.value(params.omamer)
    outputdir = Channel.value(params.outputdir)

    hogmap = omamer_run(proteomes, omamer_db, num_threads_omamer,omamer,outputdir)
    // hogmap = Channel.fromPath("./hogmap/*", type: 'any')
    // hogmap.view{ "${it}"}

    rhogs = inferrhog(hogmap.collect(), gethog3)
    // rhogs.flatten().view{" rhogs ${it}"}

    //rhogs_check = Channel.fromPath("./rhogs_all/*", type: 'any')
    // rhogs_check.view{"rhogs ${it}"}

    (rhogs_rest, rhogs_big) = rhog_distributor(rhogs, gethog3)

    rhogs_rest.flatten().view{" rhogs rest ${it}"}
    rhogs_big.flatten().view{" rhogs big ${it}"}


    //rhogs_batch = rhog_distributor(rhogs, gethog3)
    //rhogs_batch.view{"rhogs_batch ${it}"}


    //rhogs_rest_collect = rhogs_rest_.collect()
   // if ( rhogs_rest_collect.ifEmpty()){
    //rhogs_rest = Channel.fromPath("./rhogs_rest/*", type: 'any')
    //rhogs_rest.view{"rhogsrest ${it}"}
    pickle_rest_rhog = hog_rest(rhogs_rest, gethog3)
    pickle_rest_rhog.flatten().view{" pickle_rest_rhog rest ${it}"}

    pickle_big_rhog = hog_big(rhogs_big, gethog3)
    pickle_big_rhog.flatten().view{" pickle_big_rhog rest ${it}"}
    // pickle_rest_rhog.view{"pi_rest_rhog ${it}"}


    prb = pickle_big_rhog.collect()
    prr = pickle_rest_rhog.collect()
    all_pickles = prb.mix(prr)
    ortho = collect_orthoxml(all_pickles.collect(), gethog3)

}


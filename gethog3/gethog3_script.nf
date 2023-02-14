params.proteomes = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/proteome/*"
params.omamer_db= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/archive/Primates.h5" // LUCA.h5"
params.species_tree= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/tree.nwk"
params.num_threads_omamer= 2

process omamer{
  publishDir "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/hogmap/"   // , mode: 'copy'
  input:
  path proteomes
  val omamer_db
  val num_threads_omamer
  output:
  path "*.hogmap"
  script:
  """
   omamer search --db $omamer_db  --query $proteomes --nthreads $num_threads_omamer  --out ${proteomes}.hogmap
  """
}

process inferrhog{
  publishDir "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/rhogs_all/"
  input:
  path hogmap
  output:
  path "*.fa"
  script:
  """
   python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/gethog3/gethog3/infer_rhog.py
  """
}

process rhog_distributor{
  publishDir "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/"

  input:
  path rhogs
  output:
  path "rhogs_rest/*"
  path "rhogs_big/*" //"${rhogs_big}/*.fa"
  script:
  """
   python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/gethog3/gethog3/rhog_distributor.py
  """
}


process hog_big{
  publishDir "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/pi_big_rhog/"
  input:
  path rhogs_big_i //"$rhogs_big/*.fa"
  output:
  path "*.pickle"
  //path "pi_big_subhog/*"

  script:
  """
   python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/gethog3/gethog3/infer_folder.py  $rhogs_big_i False pi_big rhogs_big
  """
}


process hog_rest{
  publishDir "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/pi_rest_rhog/"
  input:
  path rhogs_rest_i //"$rhogs_big/*.fa"
  output:
  path "*.pickle"

  script:
  """
   python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/gethog3/gethog3/infer_folder.py  $rhogs_rest_i False pi_rest rhogs_rest
  """
}




workflow {
  proteomes = Channel.fromPath(params.proteomes,  type:'any' ,checkIfExists:true)
  omamer_db = Channel.value(params.omamer_db)
  num_threads_omamer = Channel.value(params.num_threads_omamer)
  hogmap = omamer(proteomes, omamer_db, num_threads_omamer)
  // hogmap = Channel.fromPath("./hogmap/*", type: 'any')
  // hogmap.view{ "${it}"}

  rhogs = inferrhog(hogmap.collect())
  //rhogs_check = Channel.fromPath("./rhogs_all/*", type: 'any')
  // rhogs_check.view{"rhogs ${it}"}

  println "rhog_distributor started"
  (rhogs_rest_, rhogs_big_) = rhog_distributor(rhogs)
  rhogs_rest = Channel.fromPath("./rhogs_rest/*", type: 'any')
  rhogs_rest.view{"rhogsrest ${it}"}
  pi_rest_rhog = hog_rest(rhogs_rest)
  pi_rest_rhog.view{"pi_rest_rhog ${it}"}


  rhogs_big = Channel.fromPath("./rhogs_big/*", type: 'any')
  rhogs_big.view{"rhogsbig ${it}"}
  pi_big_rhog = hog_big(rhogs_big)
  pi_big_rhog.view{"pi_big_rhog ${it}"}

// collector

}


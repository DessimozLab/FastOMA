
params.proteomes = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/proteome/*"
params.omamer_db= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/archive/Primates.h5" // LUCA.h5"
params.species_tree= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/tree.nwk"
params.num_threads_omamer= 2

process omamer{
  publishDir "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/hogmap/"  // , mode: 'copy'

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
  //output:
  // path "*.fa"

  script:
  """
   python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/gethog3/gethog3/_infer_rhog.py
  """
}


workflow {
  proteomes = Channel.fromPath(params.proteomes,  type:'any' ,checkIfExists:true)
  omamer_db = Channel.value(params.omamer_db)
  num_threads_omamer = Channel.value(params.num_threads_omamer)
  hogmap = omamer(proteomes, omamer_db, num_threads_omamer)
  hogmap.view{ "${it}"}
  // hogmap = Channel.fromPath("/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/hogmap/*")

  inferrhog(hogmap.collect())
  //rhogs = inferrhog(hogmap)
  //rhogs.view{ "${it}"}
}




params.proteomes = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/proteome/*"
params.omamer_db= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/archive/Primates.h5"
params.species_tree= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/tree.nwk"

process omamer{
  publishDir "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/hogmap/"

  input:
  path proteomes
  path omamer_db

  output:
  path "*.hogmap"

  script:
  // wc -l $proteomes > ${proteomes}.txt
  """
   omamer search --db $omamer_db  --query $proteomes --nthreads 4 --out ${proteomes}.hogmap
  """
}


workflow {
  proteomes = Channel.fromPath(params.proteomes,  type:'any' ,checkIfExists:true)
  hogmap = omamer(proteomes)
  hogmap.view{ "${it}"}
}



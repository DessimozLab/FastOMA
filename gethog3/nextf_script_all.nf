


params.proteomes = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/proteome/*"
// params.hogmap = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/hogmap/"

process omamer{
  publishDir "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/hogmap/"

  input:
  path proteomes

  output:
  path "*.txt"

  script:
  """
   wc -l $proteomes > ${proteomes}.txt
  """
}


workflow {
  proteomes = Channel.fromPath(params.proteomes,  type:'any' ,checkIfExists:true)
  hogmap = omamer(proteomes)
  hogmap.view{ "${it}"}
}



params.inputs = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf/rhogs/*"


process gHOG {

  cpus 1      // 8
  time {2.h}  // 8.h
  memory {20.GB}


  debug true
  input:
  path rhog_fa

  script:
  """
  python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/gethog3_nf/gethog3_30nov/infer_folder.py  $rhog_fa  True
  """
}

workflow {
  rhog_files = Channel.fromPath(params.inputs,  type:'any' ,checkIfExists:true)
  gHOG(rhog_files)
}
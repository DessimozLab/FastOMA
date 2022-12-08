
// improve : add rhog inference, rhog grouping
// check folder exitst
// resubmit some jobs

params.inputs_rest = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf/rhogs_rest/*"
params.inputs_big = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf/rhogs_big/*"

process start {

  script:
  """
  mkdir /work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf/pick_big_subhog/
  mkdir /work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf/pick_big_rhog/
  mkdir /work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf/pick_rest_subhog/
  mkdir /work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf/pick_rest_rhog/
  """

}

process qfhogbg {
  cpus 7
  time {8.h}
  memory {30.GB}

  debug true
  input:
  path rhog_fa

  script:
  """
  python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/gethog3/infer_folder.py  $rhog_fa  True pick_big rhogs_big
  """
}

process qfhogrs {
  cpus 1      // 8
  time {2.h}  // 8.h
  memory {20.GB}

  debug true
  input:
  path rhog_fa

  script:
  """
  python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/gethog3/infer_folder.py  $rhog_fa  False pick_rest rhogs_rest
  """
}

workflow {
 // start()
  rhog_files_rest = Channel.fromPath(params.inputs_rest,  type:'any' ,checkIfExists:true)
  qfhogrs(rhog_files_rest)
  rhog_files_big = Channel.fromPath(params.inputs_big,  type:'any' ,checkIfExists:true)
  qfhogbg(rhog_files_big)
}

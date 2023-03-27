

// improve : add rhog inference, rhog grouping
// check folder exitst
// resubmit some jobs

params.inputs_rest = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota//working_nf/rhogs_rest/*"
params.inputs_big = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota//working_nf/rhogs_big/*"

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
  memory {40.GB}

  debug true
  input:
  path rhog_fa

  script:
  """
  python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/gethog3/infer_folder.py  $rhog_fa  True pi_big rhogs_big
  """
}

process qfhogrs {
  cpus 1      // 8
  time {1.h}  // 8.h
  memory {20.GB}

  debug true
  input:
  path rhog_fa
  output:
  path pi_rhogs_rest/*.pickl



  script:
  """
  python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/gethog3/infer_folder.py  $rhog_fa  False pi_rest rhogs_rest
  """
}

process collect_ortho{
  input:

  script:
  """
  cp /work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/working_nf/pi_big_rhog/* /work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/working_nf/pi_rest_rhog/
  python /work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/gethog3/collect_orthoxml.py
  """
}

workflow {
 // start()
 // omamer
 // main.py rhog
 input rhogs_all

  rhog_files_rest = Channel.fromPath(params.inputs_rest,  type:'any' ,checkIfExists:true)
  qfhogrs(rhog_files_rest)
  rhog_files_big = Channel.fromPath(params.inputs_big,  type:'any' ,checkIfExists:true)
  qfhogbg(rhog_files_big)
  collect_ortho()

}



input :   fasta files
omaer
ouptu




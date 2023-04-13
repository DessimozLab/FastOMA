

params.input_folder = "./in_folder/"
params.output_folder = "./out_folder/"
params.proteome_folder = params.input_folder + "/proteome"
params.proteomes = params.proteome_folder + "/*"
params.hogmap_input_folder = params.input_folder + "/hogmap_input_folder"


params.hogmap_folder = params.output_folder + "/hogmap"
params.rhogs_folder = params.output_folder + "/rhogs_all"
params.species_tree = params.input_folder + "/species_tree.nwk"
params.pickles_rhogs_folder = params.output_folder + "/pickle_rhogs"
params.genetrees_folder = params.output_folder + "/genetrees"


process omamer_run{
  time {4.h}
  memory {50.GB}
  cpus  10
  publishDir params.hogmap_folder
  input:
  path proteomes_omamerdb_inputhog
  output:
  path "*.hogmap"
  val true      // ready_omamer_run
  script:
  //   omamer search --db ${proteomes_omamerdb[1]} --query ${proteomes_omamerdb[0]} --nthreads 1  --out ${proteomes_omamerdb[0]}.hogmap
  // cp /work/FAC/FBM/DBC/cdessim2/default/smajidi1/qfo_hogmap/${proteomes_omamerdb[0]}.hogmap .
  //
  """
    if [ -f ${proteomes_omamerdb_inputhog[2]}/${proteomes_omamerdb_inputhog[0]}.hogmap ]
    then
        cp ${proteomes_omamerdb_inputhog[2]}/${proteomes_omamerdb_inputhog[0]}.hogmap  ${proteomes_omamerdb_inputhog[0]}.hogmap
    else
        omamer search --db ${proteomes_omamerdb_inputhog[1]} --query ${proteomes_omamerdb_inputhog[0]} --nthreads 10  --out ${proteomes_omamerdb_inputhog[0]}.hogmap
    fi
  """

}


process infer_roothogs{
  publishDir  params.rhogs_folder // "${params.output_folder}/rhogs_all"
  input:
  val ready_omamer_run
  path hogmap_folder
  path proteome_folder
  output:
  path "*.fa"
  path "gene_id_dic_xml.pickle"
  val true     // ready_infer_roothogs   nextflow-io.github.io/patterns/state-dependency/
  script:
  """
   infer-roothogs  --logger-level DEBUG
  """
}


workflow {
    proteomes = Channel.fromPath(params.proteomes,  type:'any' ,checkIfExists:true)
    proteome_folder = Channel.fromPath(params.proteome_folder)
    hogmap_folder = Channel.fromPath(params.hogmap_folder)
    rhogs_folder = Channel.fromPath(params.rhogs_folder)

    genetrees_folder = Channel.fromPath(params.genetrees_folder)
    hogmap_input_folder = Channel.fromPath(params.hogmap_input_folder)

    pickles_rhogs_folder =  Channel.fromPath(params.pickles_rhogs_folder)
    omamerdb = Channel.fromPath(params.input_folder+"/omamerdb.h5")
    // proteomes.view{"prot ${it}"}

    proteomes_omamerdb = proteomes.combine(omamerdb)
    proteomes_omamerdb_inputhog = proteomes_omamerdb.combine(hogmap_input_folder)
    proteomes_omamerdb_inputhog.view{" rhogsbig ${it}"}

    (hogmap, ready_omamer_run)= omamer_run(proteomes_omamerdb_inputhog)
    ready_omamer_run_c = ready_omamer_run.collect()

    (rhogs, gene_id_dic_xml, ready_infer_roothogs) = infer_roothogs(ready_omamer_run_c, hogmap_folder, proteome_folder)


 }
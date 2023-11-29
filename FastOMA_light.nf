
// NXF_WRAPPER_STAGE_FILE_THRESHOLD='50000'

params.input_folder = "testdata/in_folder"
params.output_folder = "out_folder/"
params.proteome_folder = params.input_folder + "/proteome"
params.hogmap_in = params.input_folder + "/hogmap_in"
params.splice_folder = params.input_folder + "/splice"
params.species_tree = params.input_folder + "/species_tree.nwk"



// output subfolder definition
params.genetrees_folder = params.output_folder + "/genetrees"
params.hogmap_folder = params.output_folder + "/hogmap"

params.temp_output = params.output_folder +"/temp_output" //"/temp_omamer_rhogs"



if (params.help) {
    log.info """
    ===========================================
      FastOMA -- PIPELINE
    ===========================================
    Usage:
    Run the pipeline with default parameters:
    nexflow run FastOMA.nf

    Run with user parameters:
    nextflow run FastOMA.nf --input_folder {input.dir}  --output_folder {results.dir}

    Mandatory arguments:
        --input_folder          Input data folder. Defaults to ${params.input_folder}. This folder
                                must contain the proteomes (in a subfolder named 'proteome') and
                                a species tree file. Optionally the folder might contain
                                 - a sub-folder 'splice' containing splicing variant mappings
                                 - a sub-folder 'hogmap_in' containing precomputed OMAmer
                                   placement results for all proteomes

                                All sub-folders and sub-files can also be placed in orther
                                locations if you provide alternative values for them (see below on
                                optional arguments section).

        --output_folder         Path where all the output should be stored. Defaults to
                                ${params.output_folder}


    Profile selection:
        -profile                FastOMA can be run using several execution profiles. The default
                                set of available profiles is
                                 - docker       Run pipeline using docker containers. Docker needs
                                                to be installed on your system. Containers will be
                                                fetched automatically from dockerhub. See also
                                                additional options '--container_version' and
                                                '--container_name'.

                                 - singlularity Run pipeline using singularity. Singularity needs
                                                to be installed on your system. On HPC clusters,
                                                it often needs to be loaded as a seperate module.
                                                Containers will be fetched automatically from
                                                dockerhub. See also additional options
                                                '--container_version' and '--container_name'.

                                 - conda        Run pipeline in a conda environment. Conda needs
                                                to be installed on your system. The environment
                                                will be created automatically.

                                 - standard     Run pipeline on your local system. Mainly intended
                                                for development purpose. All dependencies must be
                                                installed in the calling environment.

                                 - slurm_singularity
                                                Run pipeline using SLURM job scheduler and
                                                singularity containers. This profile can also be a
                                                template for other HPC clusters that use different
                                                schedulers.

                                 - slurm_conda  Run pipeline using SLURM job scheduler and conda
                                                environment.

                                Profiles are defined in nextflow.config and can be extended or
                                adjusted according to your needs.


    Additional options:
        --proteome_folder       Overwrite location of proteomes (default ${params.proteome_folder})
        --species_tree          Overwrite location of species tree file (newick format).
                                Defaults to ${params.species_tree}
        --splice_folder         Overwrite location of splice file folder. The splice files must be
                                named <proteome_file>.splice.
                                Defaults to ${params.splice_folder}
        --omamer_db             Path or URL to download the OMAmer database from.
                                Defaults to ${params.omamer_db}
        --hogmap_in             Optional path where precomputed omamer mapping files are located.
                                Defaults to ${params.hogmap_in}

    Flags:
        --help                  Display this message
        --debug_enabled         Store addtional information that might be helpful to debug in case
                                of a problem with FastOMA.

    """.stripIndent()

    exit 1
}


log.info """
===========================================
  FastOMA -- PIPELINE
===========================================

 Project : ${workflow.projectDir}
 Git info: ${workflow.repository} - ${workflow.revision} [${workflow.commitId}]
 Cmd line: ${workflow.commandLine}
 Manifest's pipeline version: ${workflow.manifest.version}

Parameters:
   input_folder              ${params.input_folder}
   proteome folder           ${params.proteome}
   species_tree              ${params.species_tree}
   splice_folder             ${params.splice_folder}        (optional)
   omamer_db                 ${params.omamer_db}
   hogmap_in                 ${params.hogmap_in}            (optional)
   
   debug_enabled             ${params.debug_enabled}
""".stripIndent()


process check_input{
    publishDir params.output_folder, mode: 'copy'
    input:
        path proteome_folder
        path hogmap_folder
        path species_tree
        path omamerdb
        path splice_folder
    output:
        path "species_tree_checked.nwk"
        val "check_completed"
    script:
        """
        fastoma-check-input --proteomes ${proteome_folder} \
                            --species-tree ${species_tree} \
                            --out-tree species_tree_checked.nwk \
                            --splice ${splice_folder} \
                            --hogmap ${hogmap_folder} \
                            --omamer_db ${omamerdb} \
                            -vv
        """
}


process omamer_run{
  time {4.h}
  publishDir params.hogmap_folder, mode: 'copy'

  input:
  tuple path(proteome), path(omamer_db), path(precomputed_hogmap_folder), val(ready)

  output:
  path "*.hogmap"

  script:
  """
    if [ -f ${precomputed_hogmap_folder}/${proteome}.hogmap ] ; then
        cp ${precomputed_hogmap_folder}/${proteome}.hogmap  ${proteome}.hogmap
    else
        omamer search --db ${omamer_db} --query ${proteome} --out ${proteome}.hogmap
    fi
  """
}


process infer_roothogs{
  publishDir = [
    path: params.temp_output,
    enabled: params.debug_enabled,
    //mode: 'copy', saveAs: { filename -> filename.equals('temp_omamer_rhogs') ? null : filename }
    ]

  input:
    path hogmaps, stageAs: "hogmaps/*"
    path proteome_folder
    path splice_folder
  output:
    path "omamer_rhogs/*"
    path "gene_id_dic_xml.pickle"
    path "selected_isoforms" , optional: true
  script:
    """
       fastoma-infer-roothogs  --proteomes ${proteome_folder} \
                               --hogmap hogmaps \
                               --splice ${splice_folder} \
                               --out-rhog-folder "omamer_rhogs/" \
                               -vv
    """
}


process batch_roothogs{
  input:
    path rhogs, stageAs: "omamer_rhogs/*"
  output:
    path "rhogs_rest/*", optional: true
    path "rhogs_big/*" , optional: true
  script:
    """
        fastoma-batch-roothogs --input-roothogs omamer_rhogs/ \
                               --out-big rhogs_big \
                               --out-rest rhogs_rest \
                               -vv
    """
}

process hog_big{
  cpus  2
  time {20.h}     // for very big rhog it might need more, or you could re-run and add `-resume`
  input:
    each rhogsbig
    path species_tree
  output:
    path "pickle_hogs"
    path "msa/*.fa" , optional: true          // msa         if write True
    path "gene_trees/*.nwk" , optional: true  // gene trees  if write True
  script:
    """
        fastoma-infer-subhogs  --input-rhog-folder ${rhogsbig}  \
                               --species-tree ${species_tree} \
                               --fragment-detection \
                               --low-so-detection \
                               --parallel
    """
}

process hog_rest{
  input:
    each rhogsrest
    path species_tree
  output:
    path "pickle_hogs"
    path "msa/*.fa" , optional: true          // msa         if write True
    path "gene_trees/*.nwk" , optional: true  // gene trees  if write True
  script:
    """
        fastoma-infer-subhogs --input-rhog-folder ${rhogsrest}  \
                              --species-tree ${species_tree} \
                              --fragment-detection \
                              --low-so-detection
                              #--out pickle_hogs
    """
}


process collect_subhogs{
  publishDir params.output_folder, mode: 'copy'
  input:
    path pickles, stageAs: "pickle_folders/?"
    path "gene_id_dic_xml.pickle"
    path rhogs, stageAs: "omamer_rhogs/*"
  output:
    path "output_hog.orthoxml"
    path "OrthologousGroupsFasta"
    path "OrthologousGroups.tsv"
    path "rootHOGs.tsv"
  script:
    """
        fastoma-collect-subhogs --pickle-folder pickle_folders/ \
                                --roothogs-folder omamer_rhogs/ \
                                --gene-id-pickle-file gene_id_dic_xml.pickle \
                                --out output_hog.orthoxml \
                                --marker-groups-fasta OrthologousGroups.tsv \
                                --roothog-tsv rootHOGs.tsv \
                                -vv
    """
}

workflow {
    proteome_folder = Channel.fromPath(params.proteome_folder, type: "dir", checkIfExists:true).first()
    proteomes = Channel.fromPath(params.proteome_folder + "/*", type:'any', checkIfExists:true)
    species_tree = Channel.fromPath(params.species_tree, type: "file", checkIfExists:true).first()

    splice_folder = Channel.fromPath(params.splice_folder, type: "dir")
    genetrees_folder = Channel.fromPath(params.genetrees_folder, type: 'dir')
    hogmap_in = Channel.fromPath(params.hogmap_in, type:'dir')

    omamerdb = Channel.fromPath(params.omamer_db)
    (species_tree_checked_, ready_input_check) = check_input(proteome_folder, hogmap_in, species_tree, omamerdb, splice_folder)
    omamer_input_channel = proteomes.combine(omamerdb).combine(hogmap_in).combine(ready_input_check)
    hogmap = omamer_run(omamer_input_channel)

    (omamer_rhogs, gene_id_dic_xml, ready_infer_roothogs) = infer_roothogs(hogmap.collect(), proteome_folder, splice_folder)

    (rhogs_rest_batches, rhogs_big_batches) = batch_roothogs(omamer_rhogs)

    (pickle_big_rhog, msa_out_big, genetrees_out_rest) = hog_big(rhogs_big_batches.flatten(), species_tree)
    (pickle_rest_rhog,  msas_out_rest, genetrees_out_test) = hog_rest(rhogs_rest_batches.flatten(), species_tree)
    channel.empty().concat(pickle_big_rhog, pickle_rest_rhog).set{ all_rhog_pickle }

    (orthoxml_file, OrthologousGroupsFasta, OrthologousGroups_tsv, rootHOGs_tsv)  = collect_subhogs(all_rhog_pickle.collect(), gene_id_dic_xml, omamer_rhogs)

}

//workflow.onComplete {
//    println "Completed at    : $workflow.complete"
//    println "Duration        : $workflow.duration"
//    println "Output in       : $params.output_folder"
//    println "Nextflow report : $params.statsdir"
//	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
//}


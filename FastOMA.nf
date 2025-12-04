// Import nf-schema functions
include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'


//Set dynamic defaults for input/output paths before validation
params.input_folder    = params.input_folder ?: "${projectDir}/testdata/in_folder"
params.proteome_folder = params.proteome_folder ?: "${params.input_folder}/proteome"
params.hogmap_in       = params.hogmap_in ?: "${params.input_folder}/hogmap_in"
params.splice_folder   = params.splice_folder ?: "${params.input_folder}/splice"
params.species_tree    = params.species_tree ?: "${params.input_folder}/species_tree.nwk"

/*

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
        --fasta_header_id_transformer
                                choice of transformers of input proteome fasta header
                                to reported IDs in output files (e.g. orthoxml files)
                                Defaults to '${params.fasta_header_id_transformer}', and can be set to
                                  - noop         : no transformation (input header == output header)
                                  - UniProt      : extract accession from uniprot header
                                                   e.g. '>sp|P68250|1433B_BOVIN' --> 'P68250'

    Algorithmic parameters:
        --nr_repr_per_hog       The maximum number of representatives per subhog to keep during the
                                inference. Higher values lead to slighlty higher runtime.
                                Default to ${params.nr_repr_per_hog}.
        --filter_method         The applied filtering method on the MSAs before tree building.
                                must be one of "col-row-threshold", "col-elbow-row-threshold", "trimal".
                                Defaults to ${params.filter_method}.
        --min_sequence_length   Minimum length of a sequence to be considered for orthology
                                inference. Too short sequences tend to be problematic.
                                Defaults to ${params.min_sequence_length}.


    Flags:
        --help                  Display this message
        --debug_enabled         Store addtional information that might be helpful to debug in case
                                of a problem with FastOMA.
        --write_msas            MSAs used during inference of subhogs will be stored at
                                every taxonomic level.
        --write_genetrees       Inferred gene trees will be stored at every taxonomic level.
        --force_pairwise_ortholog_generation
                                Force producing the pairwise orthologs.tsv.gz file even if the
                                dataset contains many proteomes. By default, FastOMA produces the
                                pairwise ortholog file only if there are at most 25 proteomes in
                                the dataset.
        --report                Produce nextflow report and timeline and store in in
                                $params.statdir

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
   proteome folder           ${params.proteome_folder}
   species_tree              ${params.species_tree}
   splice_folder             ${params.splice_folder}
   omamer_db                 ${params.omamer_db}
   hogmap_in                 ${params.hogmap_in}
   fasta_header_id_transformer    ${params.fasta_header_id_transformer}

   filter_method             ${params.filter_method}
   filter_gap_ratio_row      ${params.filter_gap_ratio_row}
   filter_gap_ratio_col      ${params.filter_gap_ratio_col}
   nr_repr_per_hog           ${params.nr_repr_per_hog}
   min_sequence_length       ${params.min_sequence_length}

   debug_enabled             ${params.debug_enabled}
   report                    ${params.report}
   force_pairwise_ortholog_generation    ${params.force_pairwise_ortholog_generation}
""".stripIndent()
*/

process fetchTestData {
    // Cache in a dedicated cache directory, not input_folder
    storeDir "${params.test_data_cache ?: "${launchDir}/.test-datasets"}"
    tag "fetch data from ${url}"

    input:
    val url

    output:
    path "${dataset_name}", emit: testDataDir

    script:
    dataset_name=url.tokenize('/').last().replaceAll(/\.(tar\.gz|tgz|zip)\?.*$/, '')
    """
    python3 - <<'EOF'
import os
import requests
import tarfile

url = "${url}"
outfile = "dataset.tar.gz"
outdir = "${dataset_name}"

os.makedirs(outdir, exist_ok=True)

# Download dataset
print("Downloading dataset from:", url)
with requests.get(url, stream=True) as r:
    r.raise_for_status()
    with open(outfile, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)

# Extract tar.gz
print("Extracting dataset to:", outdir)
with tarfile.open(outfile, "r:gz") as tar:
    tar.extractall(outdir)

print("Download and extraction completed.")
EOF
    """
}


process check_input{
  label "process_single"
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
  label "process_single"
  tag "$proteome"

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
        omamer search -n 10 --db ${omamer_db} --query ${proteome} --out ${proteome}.hogmap
    fi
  """
}


process infer_roothogs{
  label "process_medium"
  
  publishDir params.temp_output, enabled: params.debug_enabled

  input:
    path hogmaps, stageAs: "hogmaps/*"
    path proteome_folder
    path splice_folder

  output:
    path "omamer_rhogs"
    path "gene_id_dic_xml.pickle"
    path "selected_isoforms" , optional: true

  script:
    """
       fastoma-infer-roothogs  --proteomes ${proteome_folder} \
                               --hogmap hogmaps \
                               --splice ${splice_folder} \
                               --out-rhog-folder "omamer_rhogs" \
                               --min-sequence-length ${params.min_sequence_length} \
                               -vv
    """
}


process batch_roothogs{
  label "process_single"

  input:
    path rhogs
  
  output:
    path "rhogs_rest/*", optional: true, emit: rhogs_rest_batches
    path "rhogs_big/*" , optional: true, emit: rhogs_big_batches
  
  script:
    """
        fastoma-batch-roothogs --input-roothogs omamer_rhogs/ \
                               --out-big rhogs_big \
                               --out-rest rhogs_rest \
                               -vv
    """
}


process hog_big{
  cpus { 4 }
  memory { 
    def max_filesize = Utils.getMaxFileSize(rhogsbig)
    def mem_base = Utils.mem_cat( max_filesize, nr_species as int )
    return [6.GB, mem_base].max() * params.memory_multiplier * task.attempt
  }
  time {
    def max_filesize = Utils.getMaxFileSize(rhogsbig)
    def time_base = Utils.time_cat(max_filesize, nr_species as int)
    return [2.h, time_base].max() * params.time_multiplier * task.attempt 
  }

  publishDir path: params.temp_output, enabled: params.debug_enabled, pattern: "pickle_hogs"
  publishDir path: params.msa_folder, enabled: params.write_msas, pattern: "*fa"
  publishDir path: params.genetrees_folder, enabled: params.write_genetrees, pattern: "*nwk"
  publishDir path: params.genetrees_folder, enabled: params.write_genetrees, pattern: "*tsv"
  publishDir path: params.genetrees_folder, enabled: params.write_genetrees, pattern: "*tsv.gz"

  input:
    each rhogsbig
    path species_tree
    val nr_species

  output:
    path "pickle_hogs"
    path "*.fa" , optional: true          // msa         if write True
    path "*.nwk" , optional: true  // gene trees  if write True
    path "*.tsv", optional: true
    path "*.tsv.gz", optional: true

  script:
    """
        fastoma-infer-subhogs  --input-rhog-folder ${rhogsbig}  \
                               --species-tree ${species_tree} \
                               --output-pickles pickle_hogs \
                               --parallel  \
                               -vv \
                               --msa-filter-method ${params.filter_method} \
                               --gap-ratio-row ${params.filter_gap_ratio_row} \
                               --gap-ratio-col ${params.filter_gap_ratio_col} \
                               --number-of-samples-per-hog ${params.nr_repr_per_hog} \
                               ${ params.write_msas ? "--msa-write" : ""} \
                               ${ params.write_genetrees ? "--gene-trees-write" : ""}
    """
}

process hog_rest{
  label "process_single"

  publishDir path: params.temp_output, enabled: params.debug_enabled, pattern: "pickle_hogs"
  publishDir path: params.msa_folder, enabled: params.write_msas, pattern: "*fa"
  publishDir path: params.genetrees_folder, enabled: params.write_genetrees, pattern: "*nwk"
  publishDir path: params.genetrees_folder, enabled: params.write_genetrees, pattern: "*tsv"
  publishDir path: params.genetrees_folder, enabled: params.write_genetrees, pattern: "*tsv.gz"

  input:
    each rhogsrest
    path species_tree
  output:
    path "pickle_hogs"
    path "*.fa" , optional: true   // msa         if write True
    path "*.nwk" , optional: true  // gene trees  if write True
    path "*.tsv", optional: true
    path "*.tsv.gz", optional: true
  script:
    """
        fastoma-infer-subhogs --input-rhog-folder ${rhogsrest}  \
                              --species-tree ${species_tree} \
                              --output-pickles pickle_hogs \
                              -vv \
                              --msa-filter-method ${params.filter_method} \
                              --gap-ratio-row ${params.filter_gap_ratio_row} \
                              --gap-ratio-col ${params.filter_gap_ratio_col} \
                              --number-of-samples-per-hog ${params.nr_repr_per_hog} \
                              ${ params.write_msas ? "--msa-write" : ""} \
                              ${ params.write_genetrees ? "--gene-trees-write" : ""}
    """
}


process collect_subhogs{
  label "process_high"
  
  publishDir params.output_folder, mode: 'copy'

  input:
    path pickles, stageAs: "pickle_folders/?"
    path "gene_id_dic_xml.pickle"
    path rhogs
    path species_tree
    val  id_transform

  output:
    path "FastOMA_HOGs.orthoxml"
    path "OrthologousGroupsFasta"
    path "OrthologousGroups.tsv"
    path "RootHOGs.tsv"
    path "RootHOGsFasta"

  script:
    """
        fastoma-collect-subhogs --pickle-folder pickle_folders/ \
                                --roothogs-folder omamer_rhogs/ \
                                --gene-id-pickle-file gene_id_dic_xml.pickle \
                                --out FastOMA_HOGs.orthoxml \
                                --marker-groups-fasta OrthologousGroups.tsv \
                                --roothog-tsv RootHOGs.tsv \
                                --species-tree ${species_tree} \
                                --id-transform $id_transform \
                                -vv
    """
}

process extract_pairwise_ortholog_relations {
  label "process_medium"
  publishDir params.output_folder, mode: 'copy'
  input:
    path orthoxml
    val nr_species
  output:
    path "orthologs.tsv.gz"
  script:
    """
        fastoma-helper -vv pw-rel --orthoxml $orthoxml \
                                  --out orthologs.tsv.gz \
                                  --type ortholog

    """
}


process fastoma_report {
  label "process_medium"

  publishDir params.output_folder, mode: 'copy'

  input:
    path notebook
    path orthoxml
    path proteome_folder
    path species_tree_checked

  output:
    path "report.ipynb"
    path "report.html"
    path "*.html"

  script:
    """
    if ! papermill --version ; then 
        >&2 echo "papermill dependency not found!"
        >&2 echo "Ensure you have installed fastoma with the 'report' feature\n (e.g. pip install fastoma[report])"
        exit 1
    fi
    papermill $notebook \
              report.ipynb \
              -p output_folder "./" \
              -p min_sequence_length ${params.min_sequence_length} \
              -p proteome_folder "$proteome_folder"

    jupyter nbconvert --to html report.ipynb
    """
}

workflow {
    // Print help message if requested
    if (params.help) {
        log.info paramsHelp("nextflow run FastOMA.nf")
        exit 0
    }

    // Validate input parameters
    validateParameters()

    // Print parameter summary
    log.info paramsSummaryLog(workflow)

    // Handle data source
    if (params.test_data_url) {
        // Fetch test dataset from remote URL
        log.info "Fetching test dataset from: ${params.test_data_url}"
        input_folder_path = fetchTestData(params.test_data_url)
    } else {
        input_folder_path = Channel.value(params.input_folder)
    }

    // Set up all channels based on the single input folder
    proteome_folder = input_folder_path.map { "${it}/proteome" }
    proteomes = input_folder_path.flatMap{ dir ->
        file("${dir}/proteome").listFiles().findAll {
            it.name.endsWith('.fa') || it.name.endsWith('.fasta')
        }
    }
    species_tree = input_folder_path.map { "${it}/species_tree.nwk" }
    splice_folder = input_folder_path.map { "${it}/splice" }
    hogmap_in = input_folder_path.map { "${it}/hogmap_in" }

    // Static channels
    omamerdb = Channel.fromPath(params.omamer_db)
    notebook = Channel.fromPath("$workflow.projectDir/FastOMA/fastoma_notebook_stat.ipynb", type: "file", checkIfExists: true).first()

    // Run the pipeline
    (species_tree_checked, ready_input_check) = check_input(proteome_folder, hogmap_in, species_tree, omamerdb, splice_folder)
    omamer_input_channel = proteomes
        .combine(omamerdb)
        .combine(hogmap_in)
        .combine(ready_input_check)
    hogmap = omamer_run(omamer_input_channel)
    nr_species = hogmap.count()

    (omamer_rhogs, gene_id_dic_xml, ready_infer_roothogs) = infer_roothogs(hogmap.collect(), proteome_folder, splice_folder)

    (rhogs_rest_batches, rhogs_big_batches) = batch_roothogs(omamer_rhogs)

    (pickle_big_rhog, msa_out_big, genetrees_out_rest) = hog_big(rhogs_big_batches.flatten(), species_tree_checked, nr_species)
    (pickle_rest_rhog,  msas_out_rest, genetrees_out_test) = hog_rest(rhogs_rest_batches.flatten(), species_tree_checked)
    channel.empty().concat(pickle_big_rhog, pickle_rest_rhog).set{ all_rhog_pickle }

    (orthoxml_file, OrthologousGroupsFasta, OrthologousGroups_tsv, rootHOGs_tsv)  = collect_subhogs(all_rhog_pickle.collect(), gene_id_dic_xml, omamer_rhogs, species_tree_checked, params.fasta_header_id_transformer)
    c = hogmap.count().branch{ n->
        TRUE: (n<=25 || params.force_pairwise_ortholog_generation)
        FALSE: n>25
    }
    extract_pairwise_ortholog_relations(orthoxml_file, c.TRUE)
    fastoma_report(notebook, orthoxml_file, proteome_folder, species_tree_checked)


  workflow.onComplete = {
    def String report = ( params.report ? "\nNextflow report : ${params.statsdir}" : "");
    println ""
    println "Completed at    : $workflow.complete"
    println "Duration        : $workflow.duration"
    println "Processes       : $workflow.workflowStats.succeedCount (success), $workflow.workflowStats.failedCount (failed)"
    println "Output in       : $params.output_folder" + report
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
  }
}


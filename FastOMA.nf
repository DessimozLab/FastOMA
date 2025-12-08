// Import nf-schema functions
include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'


//Set dynamic defaults for input/output paths before validation
params.input           = params.input ?: "${projectDir}/testdata/in_folder"
params.proteome_folder = params.proteome_folder ?: "${params.input}/proteome"
params.hogmap_in       = params.hogmap_in ?: "${params.input}/hogmap_in"
params.splice_folder   = params.splice_folder ?: "${params.input}/splice"
params.species_tree    = params.species_tree ?: "${params.input}/species_tree.nwk"

// Utility process to fetch remote datasets
process fetchRemoteData {
    // Cache in a dedicated cache directory
    storeDir "${params.test_data_cache ?: "${launchDir}/.test-datasets"}"
    tag "fetch data from ${url}"

    input:
    val url

    output:
    path "${dataset_name}", emit: testDataDir

    script:
    dataset_name=url.toString().tokenize('/').last().replaceAll(/\.(tar\.gz|tgz)(\?.*)?$/, '')
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


// Utility Process to extract a local archive file as input data
process extractLocalArchive {
    storeDir "${params.test_data_cache ?: "${launchDir}/.test-datasets"}"
    tag "extract ${archive.name}"

    input:
    path archive

    output:
    path "${archive.baseName}", emit: extractedDir

    script:
    """
    echo "Extracting local archive: ${archive}"
    
    mkdir -p ${archive.baseName}
    
    if [[ "${archive}" == *.zip ]]; then
        unzip -q ${archive} -d ${archive.baseName}
    else
        tar -xzf ${archive} -C ${archive.baseName}
    fi
    
    echo "Extraction completed."
    ls -la ${archive.baseName}/
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

// Helper function to detect input type
def detectInputType(input) {
    if (input.startsWith('http://') || input.startsWith('https://') || input.startsWith('ftp://')) {
        return 'url'
    } else if (input.endsWith('.tar.gz') || input.endsWith('.tgz') || input.endsWith('.zip')) {
        def archive=file(input)
        if (!archive.exists() || !archive.isFile()) {
            log.error "Input archive file does not exist: ${input}"
            exit 1
        }
        return 'archive'
    } else {
        def dir=file(input)
        if (!dir.exists() || !dir.isDirectory()) {
            log.error "Input directory does not exist: ${input}"
            exit 1
        }
        return 'directory'
    }
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

    // Detect input type 
    def inputType = detectInputType(params.input)
    log.info "Detected input type '${inputType}' for: ${params.input}"       
    if (inputType == "directory") {
        log.info "Using local input folder: ${params.input}"
        // Local/custom dataset - allow parameter overrides
        proteome_folder = Channel.value(params.proteome_folder)
        proteomes = Channel.fromPath("${params.proteome_folder}/*.{fa,fasta}", checkIfExists: true)
        species_tree = Channel.value(params.species_tree)
        splice_folder = Channel.value(params.splice_folder)
        hogmap_in = Channel.value(params.hogmap_in)
    } else {
        // Input is either a URL or an archive file - fetch and extract
        // Fetch test dataset from remote URL
        if (inputType == "url") {
            log.info "Fetching test dataset from URL: ${params.input}"
            input_path = fetchRemoteData(Channel.value(params.input))
        } else if (inputType == "archive") {
            log.info "Extracting test dataset from local archive: ${params.input}"

            input_path = extractLocalArchive(Channel.fromPath(params.input))
        }
        
        // Set up all channels based on the downloaded folder structure
        proteome_folder = input_path.map { "${it}/proteome" }
        proteomes = input_path.flatMap { dir ->
            file("${dir}/proteome").listFiles().findAll {
                it.name.endsWith('.fa') || it.name.endsWith('.fasta') || it.name.endsWith('.faa')
            }
        }
        species_tree = input_path.map { "${it}/species_tree.nwk" }
        splice_folder = input_path.map { "${it}/splice" }
        hogmap_in = input_path.map { "${it}/hogmap_in" }
    } 

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


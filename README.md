FastOMA
======
FastOMA is a scalable software package to infer orthology relationship.

Want to learn more about FastOMA and try it online, check out [FastOMA academy](https://omabrowser.org/oma/academy/module/fastOMA_2023) and FastOMA talk at ISMB 2023 on [YouTube](https://youtu.be/KGetTUMDvlA?si=efeqKKarwpIFgXyN)!

# Input and Output: 

### Input: 
1. Sets of protein sequences in FASTA format (with `.fa` extension) in the folder `proteome`.
The name of each fasta file is the name of species. Please make sure that the name of fasta records do not contain special characters including `||`.

2. Rooted Species tree in [newick format](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-newick-trees).
A rough species tree is enough, and it does not need to be binary (fully resolved). Besides, we do not need branch lengths. You could use the NCBI tree via ete3 package via `cat list_taxanomic_id.txt | ete3 ncbiquery --tree > species_tree.nwk` (see [this](http://etetoolkit.org/documentation/ete-ncbiquery/)). 
Note that the name of leaves of the tree (species name) should be the same as the file name of FASTAs (without `.fa` extension) (item 1). 
And there shouldn't be any repeated names in leaves names and internal node names. The tree should not be with quotation.  

3. The omamer database which is available for download from the [OMA browser](https://omabrowser.org/oma/current/).
The FastOMA workflow will automatically download the omamer database for LUCA (7.7 GB) if the argument `--omamer_db` is not
provided on the command line. The argument can be a local file (e.g. a previously downloaded omamer database file) or 
a URL to an alternative omamer database, e.g. a subset of the LUCA database which is smaller, like Primates with this [link](https://omabrowser.org/All/Primates-v2.0.0.h5) which is ~100MB. However, to have a broader reference gene families, we recommend to use the LUCA database if possible. 


You can see an example in the [testdata](https://github.com/DessimozLab/FastOMA/tree/main/testdata/in_folder) folder.
```
$ ls proteome
AQUAE.fa  CHLTR.fa  MYCGE.fa
$ cat species_tree.nwk
((AQUAE,CHLTR)inter1,MYCGE)inter2;
```

Besides, the internal node should not contain any special character (e.g. `\`  `/` or `space`).
The reason is that FastOMA write some files whose names contain the internal node's name. 
If the species tree does not have label for some/all internal nodes, FastOMA labels them sequentially. 
The updated tree will be stored in the output folder named as `species_tree_checked.nwk`.



### Main output:
Orthology information as HOG structure in [OrthoXML](https://orthoxml.org/) format
which can be used with [PyHAM](https://github.com/DessimozLab/pyham). Learn more about HOG [here](https://youtu.be/5p5x5gxzhZA?si=YxP-1VgKSH5e_wMS) and [here](https://omabrowser.org/oma/homologs/).
The details of output are described [below](#expected-output-structure-for-test-data).

Additionally, FastOMA generates TSV files for rootlevel HOGs (deepest level) and 
marker genes groups (one gene per species maximum) together with dumps of fasta 
files (one per rootlevel HOG / marker gene). 


# How to run FastOMA

FastOMA is implemented as a [nextflow-workflow](https://www.nextflow.io/). As such, FastOMA can be run without 
any installation steps given the system supports running either docker containers, singularity containers or has conda 
installed.

```bash
nextflow run dessimozlab/FastOMA -profile docker  --input_folder /path/to/in_folder --output_folder /path/to/out_folder 
```
You could also add specific version to be used by adding `-r v0.3.5` to the command line. Without any `-r` argument, 
always the latest available release will be used. With `-r dev` the latest development release can be used.

> [!WARNING]
> Nextflow caches pulled workflows. Git branches such as `dev` are therefore not automatically updated. You might need
> to first run `nextflow drop dessimozlab/FastOMA` before pulling again.

Nextflow will automatically fetch the [dessimozlab/FastOMA](https://github.com/dessimozlab/FastOMA) repository and starts 
the `FastOMA.nf` workflow. The `-profile` argument must be used to specify the profile to use. We support `docker`, 
`singularity` and `conda` which then automatically set up the necessary tools by downloading the required containers or creating 
a conda environment with the necessary dependencies.

See also [How to install FastOMA](#how-to-install-FastOMA) for additional ways how to install and run FastOMA. Note also the 
section on the different [profiles](#using-different-nextflow-profiles).

For more informaton on how Nextflow and Docker work together, see [here](https://www.nextflow.io/blog/2016/docker-and-nextflow.html).  

### More details on how to run
We provide for every commit of the repository a docker image for FastOMA on dockerhub. You can specify the container as 
part of the nextflow command with the parameter `container_version`. If you want to use the container of the current 
git checkout version, you can specify this in the following way:

```bash
nextflow run FastOMA.nf -profile docker \
    --container_version "sha-$(git rev-list --max-count=1 --abbrev-commit HEAD)" \
    --input_folder testdata/in_folder \
    --output_folder myresult/
```


# How to install FastOMA

There are four ways to run/install FastOMA detailed below:

### 1. Running workflow directly

The FastOMA workflow can be run directly without any installation using nextflow's ability to fetch a workflow from github. A specific version can be selected by specifying the `-r` option to nextflow to select a specific version of FastOMA:

```bash
nextflow run dessimozlab/FastOMA -r v0.3.5 -profile conda 
```

This will fetch version v0.3.5 from github and run the FastOMA workflow using the conda profile. See section [How to run fastOMA](#how-to-run-fastoma). 

### 2. Cloning the FastOMA repo and running from there

```bash
git clone https://github.com/DessimozLab/FastOMA.git
cd FastOMA
nextflow run FastOMA.nf -profile docker --container_version "sha-$(git rev-list --max-count=1 --abbrev-commit HEAD)" ...
```

### 3. Manual installation (for development) in python virtual environment

- install [mafft](https://mafft.cbrc.jp/alignment/software) and [FastTree](http://www.microbesonline.org/fasttree/) and ensure the software is accessible on the PATH.
- install python >= 3.9
- create virtual environment, activate it and install FastOMA with additional extras inside it:
  ```bash
  python3 -m venv .venv
  source .venv/bin/activate
  pip install FastOMA[report,nextflow] 
  ```
  You can also install FastOMA from a clone of the repository in editable mode with `pip install -e .[report,nextflow]`.

- run pipeline including with some testdata (For more details, see the section [How to run FastOMA on the test data](https://github.com/DessimozLab/fastoma?tab=readme-ov-file#how-to-run-fastoma-on-the-test-data) )
  ```bash
  nextflow run FastOMA.nf -profile standard --input_folder testdata/in_folder --output_folder output -with-report
  ```


### 4. Manual installation in conda/mamba environment
In the FastOMA repository, we provide a conda environment file that can be used to generate a conda / mamba 
environment. 

For Conda installation you need to first download the Miniconda installer from [this link](https://docs.anaconda.com/free/miniconda/). 

For MacOS:
```
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3.sh
bash Miniconda3.sh
```

For Linux
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o Miniconda3.sh
bashMiniconda3.sh
```

Then follow the instruction on the terminal. Finally, close and re-open the terminal and run
```
conda create -n fastoma python=3.9 --file environment-conda.yml
conda activate fastoma
```
Then, clone and install fastOMA using
```
git clone https://github.com/DessimozLab/FastOMA.git

```



Alternatively, you could use Mamba instead of Conda (which needs its own [installation](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)). 
Note that it is not encouraged to have both Mamba and conda on one system, [more info](https://stackoverflow.com/questions/76760906/installing-mamba-on-a-machine-with-conda).

```
mamba env create -n FastOMA -f environment-conda.yml
mamba activate FastOMA
```

Afterwards, you can run the workflow using nextflow (which is installed as part of the conda environment)

```
nextflow run FastOMA.nf -profile standard|slurm --input_folder /path/to/input_folder --output_folder /path/to/output
```
Note that you should use either the profile `standard` or `slurm` such the nextflow executor will use the activated environment.


## Using different nextflow profiles

Nextflow provides support to run a workflow on different infrastructures. Selection of this is done using the `-profile` argument. 
For FastOMA, we've implemented the following profiles below. Additional ones can also be created by specifying them in the `nextflow.config` file.

### Docker
With `-profile docker` one can use docker as an execution platform. It requires docker to be installed on the system (see [here](https://docs.docker.com/engine/install/)). The pipeline 
will automatically fetch missing containers from dockerhub (e.g. dessimozlab/fastoma) if not found locally. By default, the version
`latest` is used by the pipeline, however we provide images for any branch and release as well; even for every recent commit.
One can select the desired container via the `--container_version` argument

```
nextflow run FastOMA.nf -profile docker \
    --container_version "sha-$(git rev-list --max-count=1 --abbrev-commit HEAD)" \
    --input_folder testdata/in_folder \
    --output_folder myresult/
```
This will use the container that is tagged with the current commit id. Similarly, one could also use 
`--container_version "0.3.5"` to use the container with version `dessimozlab/fastoma:0.3.5` from dockerhub. Check the latest version on the [DockerHub](https://hub.docker.com/r/dessimozlab/fastoma/tags).

### Singularity
Since Docker needs administrator privileges (root access), [Singluarity](https://apptainer.org/index.html) (a.k.a Apptainer) is a good alternative. This can be installed using [Conda](https://anaconda.org/conda-forge/singularity) with `conda install conda-forge::singularity`. However, in most of the academic HPC cluster, singluarity is already installed and can be called with `module load`.
With `-profile singularity` singularity containers will be used to run the workflow. It requires singularity to 
be installed on your system. The containers are automatically pulled from dockerhub and converted to singularity 
containers. The same options as for [Docker](#docker) will be available.

### Conda
with `-profile conda`, the FastOMA workflow will create a conda environment which contains the necessary 
dependencies and use this environment to run the workflow steps. Note that this environment does not need 
to be activated manually. If you prefer to install the dependencies inside a conda or mamba environment 
yourself, this can be achieved as described in [](#manual-installation-for-development-in-python-virtual-environment).

### Slurm (with singularity/conda)
On a HPC system you typically run processes using a scheduler system such as slurm or LSF. We provide 
profiles `-profile slurm`, `-profile slurm_singularity` and `-profile slurm_conda` to run FastOMA with 
the respective engine using [slurm](https://slurm.schedmd.com/overview.html) as a scheduler system. 
If you need a different scheduler, it is quite straight forward to 
set it up in `nextflow.config` based on the existing profiles and the documentation of 
[nextflow executors](https://www.nextflow.io/docs/latest/executor.html).


# How to run FastOMA on the test data
Note : If you are using FastOMA with Docker or other profiles, check out the difference [here](#using-different-nextflow-profiles).   

First, cd to the `testdata` folder and download the omamer database (optional) and change its name to `omamerdb.h5`.
```
cd FastOMA/testdata
wget https://omabrowser.org/All/Primates-v2.0.0.h5     # 105MB
mv Primates-v2.0.0.h5    in_folder/omamerdb.h5 
```
(This is for the test however, I would suggest downloading the `LUCA-v2.0.0.h5` instead of `Primates-v2.0.0.h5` for your real analysis.).
Check the item 2 in the [input section](https://github.com/sinamajidian/FastOMA#input) for details.

Now we have such a structure in our testdata folder.
``` 
$ tree ../testdata/in_folder
   ├── omamerdb.h5
   ├── proteome
   │   ├── AQUAE.fa
   │   ├── CHLTR.fa
   │   └── MYCGE.fa
   └── species_tree.nwk
```

Finally, run the package using nextflow as below:
```
# cd FastOMA/testdata
nextflow run ../FastOMA.nf  \
         --input_folder in_folder  \
         --omamer_db in_folder/omamerdb.h5 \
         --output_folder out_folder \
         --report \
         -profile standard
```



Note that to have a comprehensive test, we set the default value of needed cpus as 10. 
If you face `.command.sh: line 2: papermill: command not found`, note that the orthology inference is finished and you have them in output folder and you may want to install `pip install -e .[report]` to have `papermill` generating the report and run the last step.


## Expected log for test data
After few minutes, the run for test data finishes. 
```
[] process > check_input ()     [100%] 1 of 1 ✔
[] process > omamer_run ()      [100%] 3 of 3 ✔
[] process > infer_roothogs ()  [100%] 1 of 1 ✔
[] process > batch_roothogs ()  [100%] 1 of 1 ✔
[] process > hog_big ()         [100%] 1 of 1 ✔
[] process > hog_rest ()        [100%] 2 of 2 ✔
[] process > collect_subhogs () [100%] 1 of 1 ✔
```

The first step is to run [OMAmer](https://github.com/DessimozLab/omamer) for finding the putative gene families (putative rootHOG) based on  kmer similarity.
Next, we write them in FASTA files, which could be used to run next steps in parrallel on each FASTA gene family.
Then, to have similar size jobs, we batch these FASTA files either as one big roothog (per job `hog_big`) or a few hundreds together as one job `hog_rest`.
These are decided based on the FASTA file size. Finally, once all jobs of `hog_big` and `hog_rest` are done, we `collect_subhog` and save all outputs.  

If the run interrupted, by adding `-resume` to the nextflow commond line, you might be able to continue your previous nextflow job.

Pro-tip. Nextflow creat a folder named `work` for storing its temprorary files. The characters in the bracket of the nextflow log (not shown here) are the short form of the folder address in `work/`
where the last task of such job were done.
e.g `[3f/2efg] process > check_input (1)` you can `cd work/3f/2efg` then use tab to complete the folder name, then you can see the temporary files of `check_input` task. In such folder there are some hidden files `.command.log/sh/run`.f


## Expected output structure for test data

The output of FastOMA includes several output files regarding orthology inference
(`OrthologousGroups.tsv`, `RootHOGs.tsv`, `FastOMA_HOGs.orthoxml`, `orthologs.tsv.gz` and `species_tree_checked.nwk`),
a jupyter notebook based report about the dataset (`report.ipynb` and `report.html`) and four folders
(`hogmap`, `OrthologousGroupsFasta`, `RootHOGsFasta` and `stats`).
  
The `hogmap` folder includes the output of [OMAmer](https://github.com/DessimozLab/omamer); each file corresponds to an input proteome.
The folder `OrthologousGroupsFasta` includes FASTA files, and all proteins inside each FASTA file are orthologous to each other. 
These could be used as gene markers for species tree inference with refined resolution, [more info](https://f1000research.com/articles/9-511).
Note that OrthologousGroups are groups of strict orthologs, with at most 1 representative per species.
Hierarchical Orthologous Groups are groups of orthologs and paralogs, defined at each taxonomic level. The file 
`FastOMA_HOGs.orthoxml` contains all the nested groups in orthoxml format. The `RootHOGs.tsv` and `RootHOGsFasta/` files contains
the groups at the deepest level.

So, following files and folders should appear in the folder `out_folder` which was the argument.
```
$tree out_folder
├── FastOMA_HOGs.orthoxml
├── hogmap
│   ├── AQUAE.fa.hogmap
│   ├── CHLTR.fa.hogmap
│   └── MYCGE.fa.hogmap
├── OrthologousGroupsFasta
│   ├── OG_0000001.fa.gz
│   ├── OG_0000002.fa.gz
│   ├── OG_0000003.fa.gz
│         ├ ...
├── OrthologousGroups.tsv
├── orthologs.tsv.gz
├── phylostratigraphy.html
├── report.html
├── report.ipynb
├── RootHOGsFasta
│   ├── HOG:0000001.fa.gz
│   ├── HOG:0000002.fa.gz
│   ├── HOG:0000003.fa.gz
│   ├ ...
├── RootHOGs.tsv
├── species_tree_checked.nwk
└── stats
    ├── pipeline_dag_<date>.html
    ├── report_<date>.html
    ├── timeline_<date>.html
    └── trace_<date>.txt
```
among which `FastOMA_HOGs.orthoxml` is the final output in [orthoXML format](https://orthoxml.org/0.4/orthoxml_doc_v0.4.html). Its content looks like this

```
<?xml version='1.0' encoding='utf-8'?>
<orthoXML xmlns="http://orthoXML.org/2011/" origin="FastOMA 0.1.6" originVersion="2024-01-10 17:36:45" version="0.5">
  <species name="MYCGE" taxonId="5" NCBITaxId="0">
    <database name="database" version="2023">
      <genes>
        <gene id="1000000001" protId="sp|P47500|RF1_MYCGE" />
        <gene id="1000000002" protId="sp|P13927|EFTU_MYCGE" />
        <gene id="1000000003" protId="sp|P47639|ATPB_MYCGE" />
            
 ...
    <orthologGroup id="HOG:0000001_1" taxonId="1">
      <score id="CompletenessScore" value="1.0" />
      <property name="OMAmerRootHOG" value="HOG:D0900115" />
      <property name="TaxRange" value="inter2" />
      <geneRef id="1000000005" />
      <orthologGroup id="HOG:0000001_2" taxonId="2">
        <score id="CompletenessScore" value="1.0" />
        <property name="TaxRange" value="inter1" />
        <geneRef id="1002000010" />
        <geneRef id="1001000009" />
      </orthologGroup>
    </orthologGroup>
  </groups>
</orthoXML>
```

If you are interested in specific gene in specific species, and wants to know 
proteins that are in the gene family, you can find its protein ID in the file `RootHOGs.tsv` using grep. 
The first column of this file `RootHOGs.tsv` shows the rootHOG ID which could be searched on the [OMA browser](https://omabrowser.org/). 
Note that some of the input genes might not appear in this file. 

To find list of genes that are orthologous to your gene of interest, you can search in the file `OrthologousGroups.tsv` 
where each line is an orthologous group. Each line corresponds to a FASTA file in the folder ` OrthologousGroupsFasta`. 


Note that some of the output files are symlink (a.k.a a symbolic link), linked to files in the folder `work` created by nextflow pipeline. 
This means that if you remove or rename the `work` and its parents folder, you will not have access to the output files anymore. 

If you are working on a large scale project, you may need to change the limitation on the number of files opened in linux using `ulimit -n 271072`. 

You can learn about OMA and FastOMA on [OMA Academy](https://omabrowser.org/oma/academy/).  


Regarding temp folders:
The folder `temp_output` includes `gene_id_dic_xml.pickle` storing mapping between gene name and gene integer ID used for orthoxml format,
`temp_omamer_rhogs` a folder that includes the fasta files of omamer-based gene families (described [here](https://github.com/DessimozLab/FastOMA#under-the-hood-what-are-fastoma-gene-families)).  

The folder `temp_pickles` includes the pickle file of orthoxml object which are final product of FastOMA for each gene family stored in `temp_omamer_rhogs`. 
These file can be empty when the gene family doesn't end up as a group (usually with size of 5 Byte). Gene trees and MSAs will be stored in `temp_pickles` 
if activated (in `_config.py` and fastOMA installed with `pip -e` ). 



### Using omamer's output
The first step of the FastOMA pipele is to run [OMAmer](https://github.com/DessimozLab/omamer). If you already have the hogmap files, you can put them in the `in_folder/hogmap_in`.
Then your structure of files will be 
```
$ tree ../testdata/
├── in_folder
│   ├── hogmap_in
│   │   ├── CHLTR.fa.hogmap
│   │   ├── MYCGE.fa.hogmap
│   ├── omamerdb.h5
│   ├── proteome
│   │   ├── AQUAE.fa
│   │   ├── CHLTR.fa
│   │   └── MYCGE.fa
│   └── species_tree.nwk
└── README.md
```
In this case, FastOMA uses two hogmap files and only run omamer for the `AQUAE.fa`. Then continue the rest of pipeline. 
Let's save the planet together with 
[green computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009324). 


### Run on a cluster 
For running on a SLURM cluster you can add `-c ../nextflow_slurm.config`  to the commond line.

```
# cd FastOMA/testdata
# rm -r out_folder work          # You may remove stuff from previous run
# ls ../FastOMA.nf 

nextflow ../FastOMA.nf  -c ../nextflow_slurm.config   --input_folder in_folder   --output_folder out_folder
```

You may need to re-run nextflow command line by adding `-resume`, if the allocated time is not enough for your dataset.

You may need to increase the number of opoened files in your system with `ulimit -n 131072` or higher as nextflow generates hundreds of files depending on the size of your input dataset.


## Handle splice files
You can put the splice files in the folder `in_folder/splice`. They should be named as `species_name.splice` for each species.
For each row of different isforoms of a preotien, FastOMA selects the best one (based on omamer family score and isoform length). 
We also use those proteins that are not in splice file but present in the FASTA proteome file. 
```
$ head HUMAN.splice 
HUMAN00001;HUMAN00002;HUMAN00003;HUMAN00004;HUMAN00005;HUMAN00006
HUMAN00007;HUMAN00008;HUMAN00009;HUMAN00010;HUMAN00011;HUMAN00012;
HUMAN00022;HUMAN00023;HUMAN00024;
HUMAN00027;HUMAN00028;HUMAN00029;HUMAN00030;HUMAN00031;HUMAN00032;HUMAN00033
HUMAN00034;HUMAN00035
HUMAN00036
HUMAN00037
```

The selected isforoms will be added as a new column to the input splice files stored as tsv at `out_folder/temp_output/selected_isoforms/`

## Under the hood: what are fastOMA gene families?
Firstly, those proteins that are mapped to the same OMAdb rootHOG (e.g. HOG:D0066142 for HOG:D0066142.1a.1a) by OMAmer are 
grouped together to create query rootHOGs (no protein from OMAdb is stored), from now on called rootHOG.
Then, as OMAmer provide us with alternative mapping, we try to merge those rootHOGs (high chance of split HOGs) that have 
many shared mappings. The query proteins of these rootHOGs will be stored in only one rootHOG. 
These will be saved as fasta files in `out_folder/temp_output/temp_omamer_rhogs` with file names format `HOG_LXXXXX.fa`. `L` is the release ID of OMADB. 
Replacing `_` with ':' gives the HOG ID which could be investigated in the [OMA Browser](https://omabrowser.org/oma/hog/HOG:D0114562/Sar/iham/).

There are some cases that only one protein is mapped to one rootHOG, called singleton (which is not good, we are hoping for orthologous groups/pairs).
Using alternative OMAmer mapping, FastOMA tries to put these to other rootHOGs. Still some will be left. 

FastOMA uses the [linclust](https://github.com/soedinglab/MMseqs2#cluster) software to find new gene families on set of unmapped proteins and singletons.
These will be saved as fasta files in `out_folder/temp_output/temp_omamer_rhogs` with file names format `HOG_clustXXXXX.fa`.
These are initial gene families that are used in `infer_subhogs` step, which could be split into a few smaller gene families. 

## Cite us

Majidian, Sina, Yannis Nevers, Ali Yazdizadeh Kharrazi, Alex Warwick Vesztrocy, Stefano Pascarelli, David Moi, Natasha Glover, Adrian M. Altenhoff, and Christophe Dessimoz. "Orthology inference at scale with FastOMA." bioRxiv (2024): 2024-01. https://www.biorxiv.org/content/10.1101/2024.01.29.577392v1.full


## Change log
- Update  v0.3.5:
  - Fixes an issue with reaching the maximum recursion limit. (#31)
  - Fixes a problem with parallel execution for big families. (#44)
- Update  v0.3.4:
  - Fixing a bug in marker gene groups extraction. Before, more than one gene per species were possible
- Update  v0.3.3: improvements for nextflow (selection of alternative versions) and updates on readme
- Update  v0.3.1: spliting HOG and sampling
- Update  v0.1.6: adding dynamic resources, additional and improved output
- Update  v0.1.5: docker, add help, clean nextflow 
- Update  v0.1.4: new gene families with linclust if mmseqs is installed, using quoted protein name to handle species chars, check input first 
- Update  v0.1.3: merge rootHOGs and handle singleton using omamer multi-hits
- Update  v0.1.2: improve rootHOG inference, splice, OMAmerv2 with multi-hits
- Release v0.1.0: improve nextflow pipeline and outputs. 
- prelease v.0.0.6: use `--fragment-detection` for `infer-subhogs` and `--low-so-detection --fragment-detection`
- prelease v.0.0.6: using input hogmpa
- prelease v.0.0.5: adding pip setup.py 
- prelease v.0.0.4: simple nextflow
- prelease v.0.0.3: with dask

FastOMA
======
FastOMA is a scalable software package to infer orthology relationship.

# Input and Output: 

### Input: 
1- Sets of protein sequences in FASTA format (with `.fa` extension) in the folder `proteome`.
The name of each fasta file is the name of species. Please make sure that the name of fasta records do not contain `||`. 


2- The omamer database which you can download [this](https://omabrowser.org/All/LUCA-v2.0.0.h5) 
which is from [OMA browser](https://omabrowser.org/oma/current/). 
This file is `13 Gb` containing all the gene families of the Tree of Life or you can download it for a subset of them, e.g. Primates (352MB). 

3- Rooted Species tree in [newick format](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-newick-trees).
A rough species tree is enough and it does not need to be binary. Besides, we do not need branch lengths. 
Note that the name of leaves of the tree (species name) should be the same as the file name of FASTAs (without `.fa` extension) (item 1). 
And there shouldn't be any repeated names in leaves names and internal node names. The tree should not be with quotation.  



You can see an example in the [testdata](https://github.com/sinamajidian/FastOMA/tree/master/testdata/in_folder) folder.
```
$ ls proteome
AQUAE.fa  CHLTR.fa  MYCGE.fa
$ cat species_tree.nwk
((AQUAE,CHLTR)inter1,MYCGE)inter2;
```

Besides, the internal node should not contain any special character (e.g. `\`  `/` or space).
The reason is that FastOMA write some files whose names contain the internal node's name. 
If the species tree does not have label for some/all internal nodes, FastOMA labels them sequentially.  


### Input check:

After installing FastOMA, you can have a initial check for your input dataset by running the following in the folder `in_folder`:

```
cd in_folder
check-fastoma-input
```



### Main output:
Orthology information as HOG strcutre in [OrthoXML](https://orthoxml.org/) format
which can be used with [PyHAM](https://github.com/DessimozLab/pyham).
The details of output are described [below](https://github.com/DessimozLab/FastOMA#expected-output-structure-for-test-data).


# How to run FastOMA
In summary, you need to 1) install FastOMA and its prerequisites (below), and  2) put the input files in the folder `in_folder` 
and 3) run FastOMA using the nextflow recipe `FastOMA_light.nf`. 
```
nextflow  FastOMA_light.nf  --input_folder /path/to/in_folder   --output_folder /path/to/out_folder 
```
The script `FastOMA_light.nf` is tailored for a few species. To run FastOMA with hundreds of species, please use `FastOMA.nf`. 


# How to install FastOMA

## prerequisites

First, we create a fresh [conda](https://docs.conda.io/en/latest/miniconda.html) environment.
```
conda create --name FastOMA python=3.9
conda activate FastOMA
python -m pip install --upgrade pip
```
You may use conda to install [fasttree](http://www.microbesonline.org/fasttree/), [mafft](http://mafft.cbrc.jp/alignment/software/). and [openjdk](https://jdk.java.net/java-se-ri/17) (the alternative for Java 11< version <17 which is needed for nextflow). 
```  
conda install -c bioconda mafft fasttree
conda install -c conda-forge openjdk=16
```

## How to install FastOMA 
First, download the FastOMA package:
```
wget https://github.com/DessimozLab/FastOMA/archive/refs/heads/main.zip
unzip main.zip
mv FastOMA-main FastOMA
```
Then install it
```
ls FastOMA/setup.py
python -m pip install -e FastOMA 
```

The output would be 
```
...
Running setup.py develop for FastOMA
Successfully installed Cython-3.0.1 DendroPy-4.6.1  biopython-1.81 blosc2-2.0.0 ete3-3.1.3 future-0.18.3 humanfriendly-10.0 llvmlite-0.40.1 lxml-4.9.3 msgpack-1.0.5 nextflow-23.4.3 numba-0.57.1 numexpr-2.8.5 numpy-1.24.4 omamer-0.2.6 packaging-23.1 pandas-2.0.3 property-manager-3.0 py-cpuinfo-9.0.0 pyparsing-3.1.1 pysais-1.1.0 python-dateutil-2.8.2 pytz-2023.3 scipy-1.11.2 six-1.16.0 tables-3.8.0 tqdm-4.66.1 tzdata-2023.3 verboselogs-1.7
FastOMA-0.0.6
```

You can check your installation with running one of submodules of FastOMa
``` 
infer-roothogs --version
```

You can make sure that omamer and nextflow is installed with running  
``` 
omamer -h
nextflow -h
```

If it doesn't work, you may need to have the following for nextflow to work.
```
JAVA_HOME="/path/to/jdk-17"
NXF_JAVA_HOME="/path/to/jdk-17"
export PATH="/path/to/jdk-17/bin:$PATH"
```
You can always make sure whether you are using the python that you intended to use with `which python`  and `which python3`.
If you face any difficulty during installation, feel free to create a [github issue](https://github.com/DessimozLab/FastOMA/issues), we'll try to solve it toghter.



# How to run FastOMA on the test data
Then, cd to the `testdata` folder and download the omamer database and change its name to `omamerdb.h5`.
```
cd FastOMA/testdata
wget https://omabrowser.org/All/Primates-v2.0.0.h5     # 105MB
mv Primates-v2.0.0.h5    in_folder/omamerdb.h5 
```
(This is for the test however, I would suggest downloading the `LUCA-v2.0.0.h5` instead of `Primates-v2.0.0.h5` for your real analysis.). Check the item 2 in the [input section](https://github.com/sinamajidian/FastOMA#input) for details.

Now we have such a structure in our  testdata folder.
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
nextflow ../FastOMA_light.nf  --input_folder in_folder   --output_folder out_folder  -with-report
```
The script `FastOMA_light.nf` is tailored for a few species. In real case scenario, please use `FastOMA.nf`.  
The only difference between these two scripts is the amount of CPU and memory assigned to each job. 


Note that to have a comprehensive test, we set the default value of needed cpus as 10.

## expected log for test data
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
These are decided based on the FASTA file size. Finally once all jobs of `hog_big` and `hog_rest` are done, we `collect_subhog` and save all outputs.  


If the run interrupted, by adding `-resume` to the nextflow commond line, you might be able to continue your previous nextflow job.


## expected output structure for test data

The output of FastOMA includes four files 
(`OrthologousGroupsFasta.tsv`, `rootHOGs.tsv`, `output_hog.orthoxml` and `species_tree_checked.nwk`) and four folders
(`hogmap`, `OrthologousGroupsFasta`, `temp_pickles` and `temp_output`).
  
The `hogmap` folder includes the output of [OMAmer](https://github.com/DessimozLab/omamer); each file corresponds to an input proteome.
The folder `OrthologousGroupsFasta` includes FASTA files, and all proteins inside each FASTA file are orthologous to each other. 
These could be used as gene markers for species tree inference with refined resolution, [more info](https://f1000research.com/articles/9-511).
Note that Orthologous Groups are groups of strict orthologs, with at most 1 representative per species.
Hierarchical Orthologous Groups are groups of orthologs and paralogs, defined at each taxonomic level.

So, following files and folders should appear in the folder `out_folder` which was the argument.
```
$ls out_folder
hogmap  OrthologousGroupsFasta  OrthologousGroups.tsv  output_hog.orthoxml  rootHOGs.tsv  species_tree_checked.nwk  temp_output  temp_pickles
```
among which `output_hog.orthoxml` is the final output in [orthoXML format](https://orthoxml.org/0.4/orthoxml_doc_v0.4.html). Its content looks like this

```
<?xml version="1.0" ?>
<orthoXML xmlns="http://orthoXML.org/2011/" origin="OMA" originVersion="Nov 2021" version="0.3">
   <species name="MYCGE" NCBITaxId="1">
      <database name="QFO database " version="2020">
         <genes>
            <gene id="1000000000" protId="sp|P47500|RF1_MYCGE"/>
            <gene id="1000000001" protId="sp|P13927|EFTU_MYCGE"/>
            <gene id="1000000002" protId="sp|P47639|ATPB_MYCGE"/>
            
 ...
      <orthologGroup id="HOG:B0885011_sub10003">
         <property name="TaxRange" value="inter1"/>
         <geneRef id="1002000004"/>
         <geneRef id="1001000004"/>
      </orthologGroup>
   </groups>
</orthoXML>
```

If you are interested in specific gene in specific species, and wants to know 
proteins that are in the gene family, you can find its protein ID in the file `rootHOGs.tsv` using grep. 
The first column of this file `rootHOGs.tsv` shows the rootHOG ID which could be searched on the [OMA browser](https://omabrowser.org/). 
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



### using omamer's output
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


## Run on a cluster 
For running on a SLURM cluster you can add `-c ../nextflow_slurm.config`  to the commond line.

```
# cd FastOMA/testdata
# rm -r out_folder work          # You may remove stuff from previous run
# ls ../FastOMA.nf 

nextflow ../FastOMA.nf  -c ../nextflow_slurm.config   --input_folder in_folder   --output_folder out_folder
```

You may need to re-run nextflow command line by adding `-resume`, if the allocated time is not enough for your dataset.

You may need to increase the number of opoened files in your system with `ulimit -n 131072` or higher.


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



# Downstream analysis

- High resolution tree inference

- Phylostragraphy with pyham 


## Change log
- Update  v0.1.4: new gene families with linclust if mmseqs is installed, using quoted protein name to handle species chars, check input first 
- Update  v0.1.3: merge rootHOGs and handle singleton using omamer multi-hits
- Update  v0.1.2: improve rootHOG inference, splice, OMAmerv2 with multi-hits
- Release v0.1.0: improve nextflow pipeline and outputs. 
- prelease v.0.0.6: use `--fragment-detection` for `infer-subhogs` and `--low-so-detection --fragment-detection`
- prelease v.0.0.6: using input hogmpa
- prelease v.0.0.5: adding pip setup.py 
- prelease v.0.0.4: simple nextflow
- prelease v.0.0.3: with dask

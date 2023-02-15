GetHOG v3. (Under development)
======



## prerequisites

GETHOG3 needs following software packages:  [omamer](https://github.com/DessimozLab/omamer),  [Nextflow](https://nextflow.io/),
[Biopython](https://github.com/biopython/biopython),  [ete3](http://etetoolkit.org), [fasttree](http://www.microbesonline.org/fasttree/)
and [mafft](http://mafft.cbrc.jp/alignment/software/) (multiple sequence aligner).

For installing omamer check its [github page]( [omamer](https://github.com/DessimozLab/omamer). 
For the rest, you can install them using [conda](https://docs.conda.io/en/latest/miniconda.html).
```
conda install -c conda-forge biopython ete3 
conda install -c bioconda mafft iqtree fasttree nextflow
```



# Input and Output: 

### Input: 
1- Sets of protein sequences in FASTA format (with `.fa` extension) in the folder `proteome`. The name of each fasta file is the name of species.

2- The omamer database which you can download from [here](https://omabrowser.org/oma/current/) which is this [link](https://omabrowser.org/All/LUCA.h5). 
This file is `13 Gb` containing all the gene families of the Tree of Life. 

3- Sepecies tree in nwk or phyloxml format. Note that the internal node should not contain any special character (e.g. `\`  `/` or space). 
The reason is that gethog3 write some files whose names contains the internal node's name. 

### Output:
Orthology information  as HOG structre in [OrthoXML](https://orthoxml.org/) format.


# How to config and run GETHOG3

First, download the GETHOG3 package:
```
wget https://github.com/sinamajidian/gethog3/archive/refs/heads/master.zip
unzip master.zip
mv gethog3-master gethog3
```
or clone it 
```
git clone git@github.com:sinamajidian/gethog3.git
```


GETHOG3 is based on nextflow. We consider a working folder which contains the omamer database `Primates.h5`,
the proteome folder of the species of interst `proteome` (of fa files inside),
and the speceis tree `species_tree.phyloxml`.
```
$ ls working_folder
Primates.h5  proteome  species_tree.phyloxml
```
After running the package, the outputs will appear in this working folder.  
Please set the working_folder in two places `gethog3/gethog3/_config.py` and `gethog3/gethog3/nextflow.config` :

1- The variable `working_folder` in the file `_config.py`

2- The variable `params.working_folder` in the file `nextflow.config` (the same as item 1).

Then, provide the address of the gethog3 code as `params.gethog3` in the file `nextflow.config`.

Finally you can run:
```
nextflow gethog3_script.nf
```

# TO DO
- unit test
- 


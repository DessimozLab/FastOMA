GetHOG v3. (Under development)
======



## prerequisites

GETHOG3 needs following softwares:  [omamer](https://github.com/DessimozLab/omamer),  [Nextflow](https://nextflow.io/)
[Biopython](https://github.com/biopython/biopython),  [ete3](http://etetoolkit.org), [fasttree](http://www.microbesonline.org/fasttree/)
 [mafft](http://mafft.cbrc.jp/alignment/software/) (multiple sequence aligner),
[iqtree](http://www.iqtree.org/)  (not necessary). 

For installing omamer check its [github page]( [omamer](https://github.com/DessimozLab/omamer). For the rest, you can install them using [conda](https://docs.conda.io/en/latest/miniconda.html).

```
conda install -c conda-forge biopython ete3 
conda install -c bioconda mafft iqtree fasttree nextflow
```



# Input and Output: 

### Input: 
1- Sets of protein sequences in FASTA format (with `.fa` extension) in the folder `proteome`. The name of each fasta file is the name of species.

2- The omamer database which you can download from [here](https://omabrowser.org/oma/current/) which is this [link](https://omabrowser.org/All/LUCA.h5). 
This file is `13 Gb` containing all the gene families of the Tree of Life. 

### Output:
Orthology information  as HOG structre in [OrthoXML](https://orthoxml.org/) format.


# Run

GETHOG3 is based on nextflow. We consider a working folder which contains the omamer db `Primates.h5`,
the proteome folder of the species of interst `proteome` (of fa files inside),
and the speceis tree `species_tree.phyloxml`.
```
$ ls working_folder
Primates.h5  proteome  species_tree.phyloxml
``` 
Please set the working_folder in two places:
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


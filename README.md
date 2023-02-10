GetHOG v3. (Under development)
======



## Input: 

1- Sets of protein sequences in FASTA format (with `.fa` extension) in the folder `proteome`. The name of each fasta file is the name of species.

2- The omamer database which you can download from [here](https://omabrowser.org/oma/current/) which is this [link](https://omabrowser.org/All/LUCA.h5). 
This file is `13 Gb` containing all the gene families of the Tree of Life. 

## Output:
Orthology information  as HOG structre in [OrthoXML](https://orthoxml.org/) format.



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


# Run

GETHOG3 is based on nextflow. You can run it through

```
nextflow gethog3_script.nf
```

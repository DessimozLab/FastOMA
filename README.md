GetHOG v3. (Under development)
======



## prerequisites

GETHOG3 needs following software packages:  [omamer](https://github.com/DessimozLab/omamer),  [Nextflow](https://nextflow.io/),
[Biopython](https://github.com/biopython/biopython),  [ete3](http://etetoolkit.org), [fasttree](http://www.microbesonline.org/fasttree/)
and [mafft](http://mafft.cbrc.jp/alignment/software/) (multiple sequence aligner).

You can start with [conda](https://docs.conda.io/en/latest/miniconda.html).
```
conda create --name gethog3 python=3.9
conda activate gethog3
```
For installing omamer, please check its page [github page](https://github.com/DessimozLab/omamer). (You may need to install omamer with `scipy==1.4.1 numpy==1.20.0 pytables==3.6.1`)

You can install the rest using [conda](https://docs.conda.io/en/latest/miniconda.html).
```
conda install -c conda-forge biopython ete3 
conda install -c bioconda mafft iqtree fasttree nextflow
```


# Input and Output: 

### Input: 
1- Sets of protein sequences in FASTA format (with `.fa` extension) in the folder `proteome`. The name of each fasta file is the name of species.

2- The omamer database which you can download from [here](https://omabrowser.org/oma/current/) which is this [link](https://omabrowser.org/All/LUCA.h5). 
This file is `13 Gb` containing all the gene families of the Tree of Life or a subset of them, e.g. Primates (352MB). 

3- Sepecies tree in nwk or phyloxml format. Note that the internal node should not contain any special character (e.g. `\`  `/` or space). 
The reason is that gethog3 write some files whose names contains the internal node's name. 

### Output:
Orthology information as HOG strcutre in [OrthoXML](https://orthoxml.org/) format.




# How to run GETHOG3 the test data
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
Then, cd to the `testdata` folder and download the omamer database.
```
cd gethog3/testdata
wget https://omabrowser.org/All/Primates.h5    # 352MB
mv Primates.h5  working_folder 
```
If you are using omamer database of different name, please change `params.omamer_db` in `gethog3/gethog3/nextflow.config`. 


Next, set the path to working_folder (as a global path) in two places `gethog3/gethog3/_config.py` and `gethog3/gethog3/nextflow.config` :

1- The variable `working_folder` in the file `_config.py`
```
params.working_folder= "/work/folder/gethog3/testdata/working_folder"+ "/"
params.gethog3= "/work/folder/gethog3/gethog3/"
```
2- The variable `params.working_folder` in the file `nextflow.config` (the same as item 1).
```
```

Finally run the package using nextflow as below:
```
cd gethog3/testdata
nextflow ../gethog3/gethog3_script.nf
```

Now following files and folders should appear in the `gethog3/testdata`.
```
gene_id_dic_xml.pickle  hogmap  output_hog_.orthoxml  pickle_rhogs 
 Primates.h5  proteome  rhogs_all  rhogs_big  rhogs_rest  species_tree.nwk

```
among which `output_hog_.orthoxml` is the output and its content looks like this

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



# How to config and run GETHOG3
Please first try the test data. Now you should have the GETHOG3 package.

GETHOG3 is based on nextflow. We consider a working folder which contains the omamer database `Primates.h5`,
the proteome folder of the species of interst `proteome` (of fa files inside),
and the speceis tree `species_tree.phyloxml` (or nwk).
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
nextflow gethog3/gethog3/gethog3_script.nf
```


## log changes


prelease v.0.0.1  (Feb 15 2022)

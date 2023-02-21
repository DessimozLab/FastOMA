GetHOG v3.0
======
GetHOG is a scalable software package for inferring orthology relationship.  

# Input and Output: 

### Input: 
1- Sets of protein sequences in FASTA format (with `.fa` extension) in the folder `proteome`.
The name of each fasta file is the name of species.

2- The omamer database which you can download [this](https://omabrowser.org/All/LUCA.h5) 
which is from [OMA browser](https://omabrowser.org/oma/current/). 
This file is `13 Gb` containing all the gene families of the Tree of Life or you can download it for a subset of them, e.g. Primates (352MB). 

3- Sepecies tree in [newick format](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-newick-trees).
Note that the leaves should corrospond to the name of fasta files (item 1). 
Besides, the internal node should not contain any special character (e.g. `\`  `/` or space). 
The reason is that gethog3 write some files whose names contains the internal node's name. 

### Output:
Orthology information as HOG strcutre in [OrthoXML](https://orthoxml.org/) format
which can be used with [PyHAM](https://github.com/DessimozLab/pyham)


# How to run GETHOG3
In summary, you need to 1) install GETHOG3 using pip after its prerequisites (below), and  2) put the input files in the folder `in_folder` 
and 3) run GETHOG3 using the nextflow recipe `gethog3_script.nf`. 
```
python -m pip install ./gethog3 
nextflow  gethog3_script.nf  --input_folder /path/to/in_folder   --output_folder /path/to/out_folder
```
For a detailed instruction start from prerequisites and continue to the section [How to run GETHOG3 on the test data (details)](https://github.com/sinamajidian/gethog3#how-to-run-gethog3-on-the-test-data-details).


## prerequisites

GETHOG3 needs following software packages:  [omamer](https://github.com/DessimozLab/omamer),  [Nextflow](https://nextflow.io/),
[Biopython](https://github.com/biopython/biopython), [dendropy](https://dendropy.org/),
[pyparsing](https://github.com/pyparsing/pyparsing/) , [ete3](http://etetoolkit.org), [fasttree](http://www.microbesonline.org/fasttree/)
and [mafft](http://mafft.cbrc.jp/alignment/software/).

To do so, you can start with a fresh [conda](https://docs.conda.io/en/latest/miniconda.html) environment.
```
conda create --name gethog3 python=3.9
conda activate gethog3
```
You can install the packages using [pip](https://pypi.org/).
You can always make sure whether you are using the python that you intended to use with `which python`  and `which python3`.
``` 
python -m pip install --upgrade pip
python -m pip install biopython
python -m pip install omamer
# python -m pip install pytables==3.6.1 # if you had trouble with pytables for omamer.
python -m pip install ete3  
python -m pip install nextflow
python -m pip install pyparsing
python -m pip install DendroPy 
```
You may want to install mafft and fasttree using conda as well. 
```  
conda install -c bioconda mafft fasttree
```

Note that for Nextflow you need Java (11< version <17), and you can download [openjdk](https://jdk.java.net/java-se-ri/17) instead.
You may also need to have the following for nextflow to work.
```
JAVA_HOME="/path/to/jdk-17"
NXF_JAVA_HOME="/path/to/jdk-17"
export PATH="/path/to/jdk-17/bin:$PATH"
```

You can make sure that omamer is installed with running  
``` 
omamer
```


## How to run GETHOG3 on the test data (details)
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
Then install it
```
ls gethog3/setup.py
python -m pip install -e gethog3 
```

The output would be 
```
  Preparing metadata (setup.py) ... done
Building wheels for collected packages: gethog3
  Building wheel for gethog3 (setup.py) ... done
  Created wheel for gethog3: filename=gethog3-0.0.5-py3-none-any.whl size=29386 sh
Successfully built gethog3
Installing collected packages: gethog3
Successfully installed gethog3-0.0.5
```

You can check your installation with 
``` 
infer-roothogs --version
```



Then, cd to the `testdata` folder and download the omamer database and change its name to `omamerdb.h5`.
```
cd gethog3/testdata
wget https://omabrowser.org/All/Primates.h5    # 352MB
mv Primates.h5  in_folder/omamerdb.h5 
```
(This is for the test however, I would suggest downloading the `LUCA.h5` instead of `Primates.h5` for your real analysis.)


Now we have such a structre in our  testdata folder.
``` 
$ tree ../testdata/
├── in_folder
│   ├── omamerdb.h5
│   ├── proteome
│   │   ├── AQUAE.fa
│   │   ├── CHLTR.fa
│   │   └── MYCGE.fa
│   └── species_tree.nwk
└── README.md
```


Finally, run the package using nextflow as below:
```
# cd gethog3/testdata
nextflow ../gethog3_script.nf  --input_folder in_folder   --output_folder out_folder
```

After few minutes, the run for test data finishes. 
```
[] process > omamer_run (2)      [100%] 3 of 3 ✔
[] process > infer_roothogs (1)  [100%] 1 of 1 ✔
[] process > batch_roothogs (1)  [100%] 1 of 1 ✔
[] process > hog_big (1)         [100%] 1 of 1 ✔
[] process > hog_rest (2)        [100%] 2 of 2 ✔
[] process > collect_subhogs (1) [100%] 1 of 1 ✔
```
If the run interrupted, by adding `-resume` to the nextflow commond line, you can continue your previous nextflow job. 

Then, following files and folders should appear in the folder `out_folder` which was the argument.
```
$ tree -L 1  out_folder
├── hogmap
├── output_hog.orthoxml
├── pickle_rhogs
├── rhogs_all
├── rhogs_big
└── rhogs_rest
```
among which `output_hog_.orthoxml` is the final output. Its content looks like this

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


## Run on a cluster 
For running on a SLURM cluster you can add `-c ../nextflow_slurm.config`  to the commond line.

```
# cd gethog3/testdata
# rm -r out_folder work          # You may remove stuff from previous run
# ls ../gethog3_script.nf 

nextflow ../gethog3_script.nf  -c ../nextflow_slurm.config   --input_folder in_folder   --output_folder out_folder
```


You may need to re-run nextflow command line by adding `-resume`, if the allocated time is not enough for your dataset.


## Change log

- prelease v.0.0.5: adding pip setup.py 
- prelease v.0.0.4: simple nextflow
- prelease v.0.0.3: with dask

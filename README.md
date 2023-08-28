FastOMA
======
FastOMA is a scalable software package to infer orthology relationship.



# Input and Output: 

### Input: 
1- Sets of protein sequences in FASTA format (with `.fa` extension) in the folder `proteome`.
The name of each fasta file is the name of species. Please make sure that the name of fasta records do not contain `||`. 


2- The omamer database which you can download [this](https://omabrowser.org/All/LUCA.h5) 
which is from [OMA browser](https://omabrowser.org/oma/current/). 
This file is `13 Gb` containing all the gene families of the Tree of Life or you can download it for a subset of them, e.g. Primates (352MB). 

3- Sepecies tree in [newick format](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-newick-trees).
Note that the name of leaves of the tree (species name) should be the same as the file name of fastas (without `.fa` extension) (item 1). 
You can see an example in the [testdata](https://github.com/sinamajidian/FastOMA/tree/master/testdata/in_folder) folder.
```
$ ls proteome
AQUAE.fa  CHLTR.fa  MYCGE.fa
$ cat species_tree.nwk
((AQUAE,CHLTR)inter1,MYCGE)inter2;
```

Besides, the internal node should not contain any special character (e.g. `\`  `/` or space).
The reason is that FastOMA write some files whose names contains the internal node's name. 
If the species tree does not have label for some/all internal nodes, FastOMA labels them sequentially.  

### Output:
Orthology information as HOG strcutre in [OrthoXML](https://orthoxml.org/) format
which can be used with [PyHAM](https://github.com/DessimozLab/pyham).


# How to run FastOMA
In summary, you need to 1) install FastOMA and its prerequisites (below), and  2) put the input files in the folder `in_folder` 
and 3) run FastOMA using the nextflow recipe `FastOMA_light.nf`. 
```
nextflow  FastOMA_light.nf  --input_folder /path/to/in_folder   --output_folder /path/to/out_folder 
```


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
conda install -c conda-forge openjdk
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

## How to install FastOMA 
First, download the FastOMA package:
```
wget https://github.com/DessimozLab/FastOMA/archive/refs/heads/master.zip
unzip master.zip
mv FastOMA-master FastOMA
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


# How to run FastOMA on the test data
Then, cd to the `testdata` folder and download the omamer database and change its name to `omamerdb.h5`.
```
cd FastOMA/testdata
wget https://omabrowser.org/All/Primates.h5    # 352MB
mv Primates.h5  in_folder/omamerdb.h5 
```
(This is for the test however, I would suggest downloading the `LUCA.h5` instead of `Primates.h5` for your real analysis.). Check the item 2 in the [input section](https://github.com/sinamajidian/FastOMA#input) for details.

Now we have such a structre in our  testdata folder.
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
Note that to have a comprehensive test, we set the default value of needed cpus as 10.
## expected log for test data
After few minutes, the run for test data finishes. 
```
[] process > omamer_run ()      [100%] 3 of 3 ✔
[] process > infer_roothogs ()  [100%] 1 of 1 ✔
[] process > batch_roothogs ()  [100%] 1 of 1 ✔
[] process > hog_big ()         [100%] 1 of 1 ✔
[] process > hog_rest ()        [100%] 2 of 2 ✔
[] process > collect_subhogs () [100%] 1 of 1 ✔
```
If the run interrupted, by adding `-resume` to the nextflow commond line, you can continue your previous nextflow job. 
## expected output structure for test data
Then, following files and folders should appear in the folder `out_folder` which was the argument.
```
$ls out_folder
hogmap  OrthologousGroupsFasta  OrthologousGroups.tsv  output_hog.orthoxml
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

### using omamer's output
The first step of the FastOMA pipele is to run [OMAmer](https://github.com/DessimozLab/omamer). If you already have the hogmap files, you can put them in the `in_folder/hogmap_input_folder`.
Then your structure of files will be 
```
$ tree ../testdata/
├── in_folder
│   ├── hogmap_input_folder
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
# ls ../FastOMA_light.nf 

nextflow ../FastOMA_light.nf  -c ../nextflow_slurm.config   --input_folder in_folder   --output_folder out_folder
```

You may need to re-run nextflow command line by adding `-resume`, if the allocated time is not enough for your dataset.

You may need to increase the number of opoened files in your system with `ulimit -n 131072` or higher.




## Change log
- prelease
- prelease v.0.0.6: use `--fragment-detection` for `infer-subhogs` and `--low-so-detection --fragment-detection`
- prelease v.0.0.6: using input hogmpa
- prelease v.0.0.5: adding pip setup.py 
- prelease v.0.0.4: simple nextflow
- prelease v.0.0.3: with dask

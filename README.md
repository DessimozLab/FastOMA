GetHOG v3. (Under development)
======



## Input: 

Sets of protein sequence in FASTA format. 
The name of each fasta file is the name of species.
Please put them in `omamer_search/proteome` folder.



## Steps

### Step 1: Run [omamer](https://github.com/DessimozLab/omamer) in bash:

First download the [species tree](https://omabrowser.org/All/speciestree.nwk) and [OMA database](https://omabrowser.org/All/OmaServer.h5) including Hierarchical Orthologous Groups (HOG), i.e. gene subfamily and put them in `omamer_database/oma_path`.

Then, create the HOG kmer table `omamerdb.h5`.
```
num_threads=1
cd omamer_database
omamer mkdb --db omamerdb.h5 --oma_path ./oma_path/  --nthreads ${num_threads} --root_taxon LUCA
```

Map query protein onto the HOG database.
```
num_threads=1
mkdir omamer_search/hogmap

for prot in $(ls omamer_search/proteome/| sed "s/.fa//"); do 
omamer search --db  omamer_database/omamerdb.h5 --query omamer_search/proteome/${prot}.fa \
 --nthreads ${num_threads}  --out omamer_search/hogmap/${prot}.hogmap
done
```


### Step 2:
Run rootHOG inference. 
```
python main.py rhogs
```

### Step 3:
Run HOG inference. 
```
python main.py hogs
```


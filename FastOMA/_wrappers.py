import subprocess
import logging

from ete3 import Tree

from .zoo.wrappers.aligners import mafft
from .zoo.wrappers.treebuilders import fasttree
from .zoo.wrappers.trimmers.trimal import TrimAl

logger = logging.getLogger(__name__)

seed_random = 1234  # Also in _hog_class.py
tree_tool = "fasttree"   #  "fasttree"  "iqtree"   # todo iqtree is very slow and not tested properly
rooting_mad_executable_path = "mad"  # it could be also a full address ends with mad like  /user/myfolder/mad
# mmseqs_executable_path ="mmseqs" # todo move run_linclust to _wrapper.py


def merge_msa(list_msas):
    """
    merge orthoxml_to_newick.py list of MSAs (multiple sequnce aligmnet)
    by run mafft on them.
    Each element of msa should be orthoxml_to_newick.py MultipleSeqAlignment object.

    output: merged (msa)
    """
    logger.debug("Number of items in list_msas " + str(len(list_msas)))
    logger.debug(str(list_msas[:4]) + "...")
    logger.debug("max length is " + str(max([len(i[0]) for i in list_msas])) + " .")

    #logger.debug("we are mergin subhogs"+len(list_msas))
    # logger.debug(str(list_msas[0][0].id ) + "\n")
    # SeqIO.write(list_msas ?? , genetree_msa_file_addr + ".unaligned.fa", "fasta")

    # todo using more cpus ?  (now a bit better using --thread -1)
    # sometimes better not to merge and remove gapps and do from scratch!
    wrapper_mafft_merge = mafft.Mafft(list_msas, datatype="PROTEIN")
    # todo add merge
    # if len(list_msas) < _config.num_msas_merge_mafft:
    #     wrapper_mafft_merge.options['--merge'].active = True
    # else:
    #     wrapper_mafft_merge.options['--merge'].active = False
    
    #mafft --auto Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
    wrapper_mafft_merge.options['--auto'].set_value(False) # todo we can make it as an argument in fastoma-infer-subhogs. 
    # wrapper_mafft_merge.options['--anysymbol'].active = True
    wrapper_mafft_merge.options['--anysymbol'].set_value(True)
    wrapper_mafft_merge.options['--thread'].set_value(-1) # -1 uses a largely appropriate number of threads in each step, after automatically counting the number of physical cores the computer has.
    # --randomseed
    wrapper_mafft_merge.options['--randomseed'].set_value(seed_random)
    merged = wrapper_mafft_merge()
    logger.debug("running mafft took " + str(wrapper_mafft_merge.elapsed_time))
    return merged


def infer_gene_tree(msa):
    """
    infere gene tree using fastTree for the input msa
    and write it as orthoxml_to_newick.py file

    output: gene tree in nwk format
    """
    # prot_ids = [i.id for i in msa]
    # assert len(set(prot_ids)) == len(prot_ids), "non uniq fasta record in msa "+str(prot_ids)
    # since quote is activated, the description of prots will be in the gene tree leaves. so we need to remove that.
    msa_edited = []
    for rec in msa:
        rec.description = ""
        msa_edited.append(rec)
    if tree_tool == "fasttree":
        wrapper_tree = fasttree.Fasttree(msa_edited, datatype="PROTEIN")
        wrapper_tree.options.options['-fastest'].active = True  # speed up the neighbor joining phase in fasttree & reduce memory usage  (recommended for >50,000 sequences)
        #  we don't really need fastest for small dataset and making this False didn't make qfo result better
        wrapper_tree.options.options['-quote'].active = True   # .set_value(True)  doesnt work.
        wrapper_tree.options.options['-seed'].set_value(seed_random)
        #wrapper_tree.options.options['-nt'].active = True
        # todo using more cpus ?
        # elif _config.tree_tool == "iqtree": # very slow not recommanded
        #     wrapper_tree = iqtree.Iqtree(msa, datatype="PROTEIN")
        #     wrapper_tree.options.options['-m'].set_value("LG+I+G")
        #     wrapper_tree.options.options['-nt'].set_value(1)

    result_tree1 = wrapper_tree()
    if wrapper_tree.stderr:
        logger.debug("tree inference stderr " + str(wrapper_tree.stderr))
    time_taken_tree = wrapper_tree.elapsed_time
    result_tree2 = wrapper_tree.result
    tree_nwk = result_tree2["tree"].as_string(schema='newick') #str(result_tree2["tree"])
    return tree_nwk

# to avoid confusion: right now we are using _utils_roothog.run_linclust in infer_roothogs.py
# def run_linclust(fasta_to_cluster="singleton_unmapped.fa"):

#     num_threads = 5
#     command_clust= "mmseqs easy-linclust --threads" +str(num_threads) + " " +fasta_to_cluster+"singleton_unmapped tmp_linclust"

#     logger.debug("linclust rooting started" + command_clust)
#     process = subprocess.Popen(command_clust.split(), stdout=subprocess.PIPE)
#     output, error = process.communicate()
#     #if "Error analyzing file" in str(output) or error:
#     #    try:

#     return "done"


def trim_msa(msa):
    try:
        trimal = TrimAl(msa)
        trimal.options['-automated1'].set_value(True)
        msa_out = trimal()

    except:
        # when the input msa is too much gappy
        # output trimmed msa is empty and trimal do not generate output file and raise FileNotFoundError  and  WrapperError
        # return AlignIO.read(output, 'fasta')
        msa_out = []

    return msa_out


def mad_rooting(input_tree_file_path: str):  # , mad_executable_path: str = "./mad"
    """
    Rooting a tree using MAD algorithm.
    :param input_tree_file_path: path to the input tree file.
    :param mad_excutable_path: path to the MAD executable.
    :return: rooted tree as PhyloTree object.

    https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip
        # Example usage:

    input_tree_file_path = "input_tree.nwk"
    "./mad"

    # Root the tree using MAD.
    rooted_tree = mad_rooting(input_tree_file_path, mad_executable_path)

    """
    # Create the command to run MAD.
    mad_command = f"{rooting_mad_executable_path} -i {input_tree_file_path}"
    # Run MAD and ignore the output.
    #os.system(mad_command)


    logger.debug("MAD rooting started")
    #bashCommand = f"mad {gene_tree_address}"
    process = subprocess.Popen(mad_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    if "Error analyzing file" in str(output) or error:
        raise RuntimeError("Error running MAD rooting: \n{}\n{}".format(output, error))

    elif "rooted trees written" in str(output): # 3 rooted trees written to tree_664187_Eukaryota.nwk.rooted Warning: Trees with repeating branch lengths are suspicious (3 repeating values).
        with open(input_tree_file_path + ".rooted", 'rt') as f_in:
            trees = []
            for line_tree1 in f_in:
                if line_tree1.strip():
                    trees.append(line_tree1.strip())
        # todo which one to choose ?
        rooted_tree = Tree(trees[0])
    else:  # Rooted tree written to 'tree_672375_Theria.nwk.rooted'
        rooted_tree = Tree(input_tree_file_path + ".rooted")
    logger.debug("MAD rooting finished")
    return rooted_tree


def infer_gene_tree(members_list_lowerLevel_ready, gene_tree_file_addr):
    """
    infere gene tree using fastTree for the input msa
    and write it as orthoxml_to_newick.py file


    output: gene tree in nwk format
    """
    #
    # if _config.tree_tool == "fasttree":
    #     wrapper_tree = fasttree.Fasttree(msa, datatype="PROTEIN")
    #     wrapper_tree.options.options['-fastest'].active = True   # .set_value(True)  is wrong.
    #     #wrapper_tree.options.options['-quote'].active = True
    #     #wrapper_tree.options.options['-nt'].active = True
    #
    # # todo using more cpus ?
    # # elif _config.tree_tool == "iqtree": # very slow not recommanded
    # #     wrapper_tree = iqtree.Iqtree(msa, datatype="PROTEIN")
    # #     wrapper_tree.options.options['-m'].set_value("LG+I+G")
    # #     wrapper_tree.options.options['-nt'].set_value(1)
    #
    # result_tree1 = wrapper_tree()
    # if wrapper_tree.stderr:
    #     logger_hog.debug("tree inference stderr " + str(wrapper_tree.stderr))
    # time_taken_tree = wrapper_tree.elapsed_time
    # result_tree2 = wrapper_tree.result
    # tree_nwk = str(result_tree2["tree"])
    # print(time_taken_tree)


    # tree_nwk=1


    import os
    import glob
    import json
    import datetime
    import re
    import numpy as np
    import pandas as pd
    import toytree
    import shutil
    import subprocess
    import shutil

    #from Bio.PDB import *
    from FastOMA.fold_tree import AFDB_tools
    from FastOMA.fold_tree import foldseek2tree
    #from fold_tree.src import AFDB_tools
    #from fold_tree.src import foldseek2tree
    #
    # def filter_plddt(pdb_path, thresh=.6, minthresh=.5):
    #     '''
    #     Extracts the plddt (in the beta factor column) of the first atom of each residue in a PDB file and returns bool if the pdb is accepted or not.
    #
    #     Parameters:
    #         pdb_path (str): The path to the PDB file.'''
    #
    #     thresh =  40
    #     minthresh =  0
    #     lddt = []
    #     parser = PDBParser()
    #     struc = parser.get_structure("a", pdb_path)
    #     for res in struc.get_residues():
    #         for at in res.get_atoms():
    #             lddt.append(at.get_bfactor())
    #             break
    #     if np.mean(lddt) < thresh or np.amin(lddt) < minthresh:
    #         return False
    #     else:
    #         return True
    #     return 1

    def struct_f(ids, infolder):


        structfolder = infolder + 'structs/'
        rejectedfolder = infolder + 'rejected/'
        try:
            os.mkdir(structfolder)
        except:
            print("folder exist " + structfolder)
        try:
            os.mkdir(rejectedfolder)
        except:
            print("folder exist " + rejectedfolder)

        resdf = AFDB_tools.grab_entries(ids, verbose=False) # download structures
        # as part of grab_entries : if not os.path.isfile(structfolder + uniID +'.pdb'):

        missing = [AFDB_tools.grab_struct(i, structfolder, rejectedfolder) for i in ids]
        found = glob.glob(structfolder + '*.pdb') + glob.glob(rejectedfolder + '*.pdb')
        found = {i.split('/')[-1].replace('.pdb', ''): i for i in found}
        #missing_structs = set(ids) - set(found.keys())

        return  len(set(found.keys()))


    def struct_f2(ids, infolder):
        structfolder = infolder + 'structs/'
        #rejectedfolder = infolder + 'rejected/'
        try:
            os.mkdir(structfolder)
        except:
            print("folder exist " + structfolder)
        #try:
        #    os.mkdir(rejectedfolder)
        #except:
        #    print("folder exist " + rejectedfolder)

        not_found=[]
        for id1 in ids:
            try :
                structs_folders="downloaded_structures/"
                logger_hog.debug(" *1* we are copying this file struct pdb  " + structs_folders +id1+".pdb")
                shutil.copyfile(structs_folders+id1+".pdb", structfolder+id1+".pdb")
            except:
                not_found.append(id1)
                print("struct pdb file not found"+id1)
                logger_hog.debug(" *2*  struct pdb file not found"+id1)

        #resdf = AFDB_tools.grab_entries(ids, verbose=False) # download structures
        # as part of grab_entries : if not os.path.isfile(structfolder + uniID +'.pdb'):

        #missing = [AFDB_tools.grab_struct(i, structfolder, rejectedfolder) for i in ids]
        #found = glob.glob(structfolder + '*.pdb') + glob.glob(rejectedfolder + '*.pdb')
        #found = {i.split('/')[-1].replace('.pdb', ''): i for i in found}
        #missing_structs = set(ids) - set(found.keys())

        return  len(ids)-len(not_found)

    def foldseek_dist(infolder):
        # fold_tree/foldseek/foldseek
        #
        logger_hog.debug("foldseek started")
        #command = "foldseek easy-search " + infolder + "structs/  " + infolder + "structs/ " + infolder + "allvall_1.csv " + infolder + "tmp --format-output query,target,fident,evalue,bits --exhaustive-search --alignment-type 2 -e inf"
        command = "foldseek easy-search " + infolder + "structs/  " + infolder + "structs/ " + infolder + "allvall_1.csv " + infolder + "tmp --format-output query,target,fident,evalue,bits --exhaustive-search --alignment-type 2 -e inf"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

        # output allvall_1.csv
        logger_hog.debug("foldseek finished")
        res = pd.read_table(infolder + "allvall_1.csv", header=None)
        print(res.head())
        # get the folder of the input file
        # infolder = snakemake.input[0].split('/')[:-1]
        # infolder = ''.join( [i + '/' for i in infolder])+'/'
        res[0] = res[0].map(lambda x: x.replace('.pdb', ''))
        res[1] = res[1].map(lambda x: x.replace('.pdb', ''))
        res.columns = 'query,target,fident,evalue,bits'.split(',')
        ids = list(set(list(res['query'].unique()) + list(res['target'].unique())))
        pos = {protid: i for i, protid in enumerate(ids)}
        kernels = ['fident']

        # set kernel columns to float
        for k in kernels:
            res[k] = res[k].astype(float)
        # change nan to 0
        res = res.fillna(0)
        matrices = {k: np.zeros((len(pos), len(pos))) for k in kernels}
        print(res)

        # calc kernel for tm, aln score, lddt
        for idx, row in res.iterrows():
            for k in matrices:
                matrices[k][pos[row['query']], pos[row['target']]] += row[k]
                matrices[k][pos[row['target']], pos[row['query']]] += row[k]

        output = ["foldtree_fastmemat.txt"]
        for i, k in enumerate(matrices):
            matrices[k] /= 2
            matrices[k] = 1 - matrices[k]
            print(matrices[k], np.amax(matrices[k]), np.amin(matrices[k]))
            #np.save(infolder + k + '_distmat.npy', matrices[k])
            distmat_txt = foldseek2tree.distmat_to_txt(ids, matrices[k], infolder + output[i])

        return 1

    def quicktree_f(infolder):
        delta=0
        logger_hog.debug("quicktree started")
        #command = "quicktree -i m " + infolder + "foldtree_fastmemat.txt "  # > foldtree_struct_tree.nwk
        command = "quicktree -i m " + infolder + "foldtree_fastmemat.txt "  # > foldtree_struct_tree.nwk
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        logger_hog.debug("quicktree finished")
        # treefile = foldseek2tree.postprocess("foldtree_struct_tree.nwk" , )

        outree = infolder + "foldtree_struct_tree.PP.nwk"
        t = str(output)[2:-1]
        treein = t.split("\\n")
        treestr = ' '.join([i.strip() for i in treein])
        tre = toytree.tree(treestr, format=0)
        # print(tre)
        for n in tre.treenode.traverse():
            if n.dist < 0:
                n.dist = delta
        tre.write(outree, tree_format=0)
        tree1 = tre.write(tree_format=0)

        return tree1

    from Bio.Align import MultipleSeqAlignment
    ids_dic = {}
    ids = [i.split("|")[1] for i in members_list_lowerLevel_ready]

    for prot_i in members_list_lowerLevel_ready:
        ids_dic[prot_i.split("|")[1]] = prot_i

    #
    # for in_msa in msa:
    #     # msa could be a list of MSAs. [for structure tree we don't do msa]
    #     if isinstance(in_msa, MultipleSeqAlignment):
    #         msa_i=in_msa
    #         for prot_ii in msa_i:
    #             ids_dic[prot_ii.id.split("|")[1]] = prot_ii.id
    #             ids.append(prot_ii.id.split("|")[1])
    #     else:
    #         prot_i = in_msa
    #         ids_dic[prot_i.id.split("|")[1]]=prot_i.id
    #         ids.append(prot_i.id.split("|")[1])
    #ids = [i.id.split("|")[1] for i in msa]



    time_date_raw =str(datetime.datetime.now())
    infolder ="fold_tmp/" +re.sub('[^A-Za-z0-9]+', '', time_date_raw)+"/"
    try:
        os.makedirs(infolder)
    except:
        print("folder exist "+infolder )

    num_prot_struct = struct_f2(ids,infolder)
    if num_prot_struct >1:
        foldseek_dist(infolder)

        tree_nwk_raw = quicktree_f(infolder)
        tree1= Tree(tree_nwk_raw)
        for node in tree1.traverse():
            if node.is_leaf():
                node_name_old = node.name  #
                node.name = ids_dic[node_name_old]
        tree_nwk = tree1.write()

        print("\n \n " + tree_nwk)
        a=2
        # current_time = datetime.now().strftime("%H:%M:%S")
        # for development we write the gene tree, the name of file should be limit in size in linux.
        # danger of overwriting
        # instead -> hash thing
        # ??? hashlib.md5(original_name).hexdig..it()

        if _config.gene_trees_write_all or _config.rooting_method == "mad":
            file_gene_tree = open(gene_tree_file_addr, "w")
            file_gene_tree.write(tree_nwk)
            #file_gene_tree.write(";\n")
            file_gene_tree.close()
    else:
        # not enough structures  downlaoded  to use
        tree_nwk="("+str(ids[0])+");"

    return tree_nwk



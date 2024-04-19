from Bio import SeqIO
from .zoo.wrappers.aligners import mafft
from .zoo.wrappers.treebuilders import fasttree
from .zoo.wrappers.trimmers.trimal import TrimAl
from ete3 import Tree

import subprocess

import logging

seed_random=1234 # Also in _hog_class.py
tree_tool = "fasttree"  #  "fasttree"  "iqtree"  # todo iqtree is very slow and not tested properly
rooting_mad_executable_path = "mad"  # it could be also a full address ends with mad like  /user/myfolder/mad
# mmseqs_executable_path ="mmseqs" # todo move run_linclust to _wrapper.py

logger_level = "DEBUG"            # DEBUG INFO  # TRACE  DEBUG INFO  WARN  ERROR  FATAL
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger("hog")
if logger_level == "INFO":
    logger.setLevel(logging.INFO)
if logger_level == "DEBUG":
    logger.setLevel(logging.DEBUG)



def merge_msa(list_msas, genetree_msa_file_addr, conf_infer_subhhogs):
    """
    merge orthoxml_to_newick.py list of MSAs (multiple sequnce aligmnet)
    by run mafft on them.
    Each element of msa should be orthoxml_to_newick.py MultipleSeqAlignment object.

    output: merged (msa)
    """
    logger.debug("Number of items in list_msas "+str(len(list_msas)))
    logger.debug(str(list_msas[:4])+"...")
    logger.debug("max length is "+ str(max([len(i[0]) for i in list_msas]))+" .")

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

    # wrapper_mafft_merge.options['--anysymbol'].active = True
    wrapper_mafft_merge.options['--anysymbol'].set_value(True)
    wrapper_mafft_merge.options['--thread'].set_value(-1) # -1 uses a largely appropriate number of threads in each step, after automatically counting the number of physical cores the computer has.
    # --randomseed
    wrapper_mafft_merge.options['--randomseed'].set_value(seed_random)

    merged = wrapper_mafft_merge()
    # time_duration = wrapper_mafft_merge.elapsed_time
    # print(time_duration)
    # logger.info(str(len(list_msas)) + " msas are merged with length of "+ str(len(merged)) + "  " + str (len(merged[0])))
    if conf_infer_subhhogs.msa_write:
        SeqIO.write(merged, genetree_msa_file_addr + ".fa", "fasta")

    return merged


def infer_gene_tree(msa, genetree_msa_file_addr, conf_infer_subhhogs):
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

    # current_time = datetime.now().strftime("%H:%M:%S")
    # for development we write the gene tree, the name of file should be limit in size in linux.
    # danger of overwriting
    # instead -> hash thing
    # ??? hashlib.md5(original_name).hexdig..it()

    if conf_infer_subhhogs.gene_trees_write or conf_infer_subhhogs.gene_rooting_method == "mad":
        file_gene_tree = open(genetree_msa_file_addr+".nwk", "w")
        file_gene_tree.write(tree_nwk) #file_gene_tree.write(";\n")
        file_gene_tree.close()

    return tree_nwk


def run_linclust(fasta_to_cluster="singleton_unmapped.fa"):

    num_threads = 5
    command_clust= "mmseqs easy-linclust --threads" +str(num_threads) + " " +fasta_to_cluster+"singleton_unmapped tmp_linclust"

    logger.debug("linclust rooting started" + command_clust)
    process = subprocess.Popen(command_clust.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    #if "Error analyzing file" in str(output) or error:
    #    try:

    return "done"


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
    #if verbose:
    #    print("output:\n", output)
    #    print("error:\n", error)

    if "Error analyzing file" in str(output) or error:
        raise RuntimeError("Error running MAD rooting: \n{}\n{}".format(output, error))

    elif "rooted trees written" in str(output): # 3 rooted trees written to tree_664187_Eukaryota.nwk.rooted Warning: Trees with repeating branch lengths are suspicious (3 repeating values).
        file_multiple_tree_handle= open(input_tree_file_path + ".rooted",'r')
        trees=[]
        for line_tree1 in file_multiple_tree_handle:
            if line_tree1.strip():
                trees.append(line_tree1.strip())
        # todo which one to choose ?
        rooted_tree = Tree(trees[0])
    else:#Rooted tree written to 'tree_672375_Theria.nwk.rooted'
        rooted_tree = Tree(input_tree_file_path + ".rooted")
    logger.debug("MAD rooting finished")
    return rooted_tree

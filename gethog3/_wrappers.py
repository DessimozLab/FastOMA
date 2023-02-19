
from Bio import SeqIO
from gethog3.zoo.wrappers.aligners import mafft
from gethog3.zoo.wrappers.treebuilders import fasttree
# from trimmers.trimal import TrimAl
from ete3 import Tree


from ._utils_subhog import logger_hog
from . import _config


def merge_msa(list_msas, gene_tree_file_addr):
    """
    merge orthoxml_to_newick.py list of MSAs (multiple sequnce aligmnet)
    by run mafft on them.
    Each element of msa should be orthoxml_to_newick.py MultipleSeqAlignment object.

    output: merged (msa)
    """
    logger_hog.debug(list_msas)
    # logger_hog.debug(str(list_msas[0][0].id ) + "\n")
    # SeqIO.write(list_msas ?? , gene_tree_file_addr + ".unaligned.fa", "fasta")

    wrapper_mafft_merge = mafft.Mafft(list_msas, datatype="PROTEIN")
    wrapper_mafft_merge.options['--merge'].active = True
    wrapper_mafft_merge.options['--anysymbol'].active = True
    merged = wrapper_mafft_merge()
    time_duration = wrapper_mafft_merge.elapsed_time
    # print(time_duration)
    #logger_hog.info \
    #    (str(len(list_msas)) + " msas are merged with length of "+ str(len(merged)) + "  " + str (len(merged[0])))
    if _config.gene_trees_write:
        SeqIO.write(merged, gene_tree_file_addr + "msa.fa", "fasta")

    return merged


def infer_gene_tree(msa, gene_tree_file_addr):
    """
    infere gene tree using fastTree for the input msa
    and write it as orthoxml_to_newick.py file


    output: gene tree in nwk format
    """

    if len(msa) <= 2:
        wrapper_tree = fasttree.Fasttree(msa, datatype="PROTEIN")
        wrapper_tree.options.options['-fastest'].active = True

    elif _config.tree_tool == "fasttree":
        wrapper_tree = fasttree.Fasttree(msa, datatype="PROTEIN")
        wrapper_tree.options.options['-fastest'].active = True

    elif _config.tree_tool == "iqtree":
        wrapper_tree = iqtree.Iqtree(msa, datatype="PROTEIN")
        wrapper_tree.options.options['-m'].set_value("LG+I+G")
        wrapper_tree.options.options['-nt'].set_value(1)

    result_tree1 = wrapper_tree()
    logger_hog.debug("iqtree stderr " + str(wrapper_tree.stderr))
    time_taken_tree = wrapper_tree.elapsed_time
    result_tree2 = wrapper_tree.result
    tree_nwk = str(result_tree2["tree"])
    # print(time_taken_tree)

    # current_time = datetime.now().strftime("%H:%M:%S")
    # for development we write the gene tree, the name of file should be limit in size in linux.
    # danger of overwriting
    # instead -> hash thing
    # ??? hashlib.md5(original_name).hexdig..it()

    if _config.gene_trees_write or _config.rooting_method == "mad":
        file_gene_tree = open(gene_tree_file_addr, "w")
        file_gene_tree.write(tree_nwk)
        file_gene_tree.write(";\n")
        file_gene_tree.close()

    return tree_nwk





def trim_msa(msa):


    try :
        trimal = TrimAl(msa)
        trimal.options['-automated1'].set_value(True)
        msa_out = trimal()

    except:
        # when the input msa is too much gappy
        # output trimmed msa is empty and trimal do not generate output file and raise FileNotFoundError  and  WrapperError
        # return AlignIO.read(output, 'fasta')
        msa_out = []

    return msa_out


def mad_rooting(input_tree_file_path: str, mad_executable_path: str = "./mad"):
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
    mad_command = f"{_config.rooting_mad_executable_path} -i {input_tree_file_path}"
    # Run MAD and ignore the output.
    #os.system(mad_command)

    import subprocess

    #bashCommand = f"mad {gene_tree_address}"
    process = subprocess.Popen(mad_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    #if verbose:
    #    print("output:\n", output)
    #    print("error:\n", error)

    if "Error analyzing file" in str(output) or error:
        rooted_tree = Tree(input_tree_file_path)
    else:
        rooted_tree = Tree(input_tree_file_path + ".rooted")

    return rooted_tree

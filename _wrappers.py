
from datetime import datetime
import zoo.wrappers.aligners.mafft as mafft
import zoo.wrappers.treebuilders.fasttree as fasttree

from _utils import logger_hog


def merge_msa(list_msas):
    """
    merge a list of MSAs (multiple sequnce aligmnet)
    by run mafft on them.
    Each element of msa should be a MultipleSeqAlignment object.

    output: merged (msa)
    """
    logger_hog.debug(list_msas)
    logger_hog.debug(str(list_msas[0][0].id ) + "\n")
    wrapper_mafft_merge = mafft.Mafft(list_msas, datatype="PROTEIN")
    wrapper_mafft_merge.options['--merge'].active = True
    merged = wrapper_mafft_merge()
    logger_hog.info \
        (str(len(list_msas)) + " msas are merged into one with the length of  " + str(len(merged)) + "  " + str
        (len(merged[0])))
    return merged


def infer_gene_tree(msa, gene_tree_file_addr):
    """
    infere gene tree using fastTree for the input msa
    and write it as a file


    output: gene tree in nwk format
    """
    wrapper_tree = fasttree.Fasttree(msa, datatype="PROTEIN")
    wrapper_tree.options.options['-fastest']
    result_tree1 = wrapper_tree()

    time_taken_tree = wrapper_tree.elapsed_time
    result_tree2 = wrapper_tree.result
    tree_nwk = str(result_tree2["tree"])
    # current_time = datetime.now().strftime("%H:%M:%S")
    # for development we write the gene tree, the name of file should be limit in size in linux.
    # danger of overwriting
    # instead -> hash thing
    # ??? hashlib.md5(original_name).hexdig..it()

    # as the debug==True
    # if len(gene_tree_file_addr) > 255: gene_tree_file_addr = gene_tree_file_addr[:250] + ".nwk"
    # file_gene_tree = open(gene_tree_file_addr, "w")
    # file_gene_tree.write(tree_nwk)
    # file_gene_tree.write(";\n")
    # file_gene_tree.close()

    return tree_nwk

__author__ = 'admin'

import lxml.etree as etree
from familyanalyzer.orthoxmlquery import ElementError, OrthoXMLQuery
import familyanalyzer as fa

# To parse the orthoxml and get the related etree object: etree.parse(fn_orthoxml) denote here by XML
# To get the root of your etree object: etree_object.getroot() denote here by XML_root
# The orthoxml is using two ids:
#   - The internal id denotes by the tag "id" which is found in the mapping in the beginning of the orthoxml and inside the orthologous groups
#   - The external id denotes by the tag "protId"(sometimes with the tag "geneId" or BOTH) which is found only in the begining of the orthoxml

def get_mapping_orthoxml(XML_root):

    map_ext_2_int = {}
    map_int_2_ext = {}

    for el in OrthoXMLQuery.getInputGenes(XML_root, species=None):
        int_id = None
        ext_id = []
        for _idtag,_id in el.attrib.items():
            if _idtag == "id":
                int_id = _id
            else:
                ext_id.append(str(_id))

        if int_id != None and ext_id != []:
            map_ext_2_int[frozenset(ext_id)] = str(int_id)
            map_int_2_ext[str(int_id)] = ext_id

    return map_ext_2_int, map_int_2_ext

def get_all_hogs(XML_root):
    '''
    Return a list of HOGs (etree object)
    '''
    return OrthoXMLQuery.getToplevelOrthologGroups(XML_root)

def get_hog_by_id(XML_root, hog_id):
    '''
    Return the HOG of interest (etree object)
    '''
    for OG_etree in OrthoXMLQuery.getToplevelOrthologGroups(XML_root):
        if int(OG_etree.get("id")) == int(hog_id):
            return OG_etree

def get_hogs_by_id(XML_root, list_hog_id):
    '''
    Return HOGs of interest (etree object)
    '''
    list_hog = []
    for OG_etree in OrthoXMLQuery.getToplevelOrthologGroups(XML_root):
        if int(OG_etree.get("id")) in list_hog_id:
            list_hog.append(OG_etree)
    return list_hog


def get_hogs_with_query_genes(XML_root, list_genes_you_want_to_look, use_internal_id = False, return_hog_etree = False):
    '''
    XML_root: root of the etree element
    list_genes_you_want_to_look: list of the genes external id
    use_internal_id: if the list of genes given is composed of internal id
    return_hog_etree: return the etree element of the hod instead of the id
    ext_id: external id use to parse the orthoxml

    return: each gene with its related HOG id (dictionnary)

    Protips: list(set(gene_with_its_HOG.values())) will return the list of hogs return

    '''
    map_ext_2_int, map_int_2_ext = get_mapping_orthoxml(XML_root)
    gene_with_its_HOG = {} # keys:gene & values:hog

    if use_internal_id == False:
        list_genes_int = []
        for ext_id in list_genes_you_want_to_look:
            for ext_id_keys in map_ext_2_int.keys():
                if str(ext_id) in ext_id_keys:
                    list_genes_int.append(str(map_ext_2_int[ext_id_keys]))
    else:
        list_genes_int = []
        for intid in list_genes_you_want_to_look:
            list_genes_int.append(str(intid))
    for OG_etree in get_all_hogs(XML_root):
        OG_genes = [ gene_etree for gene_etree in OG_etree.getiterator() if gene_etree.tag == "{http://orthoXML.org/2011/}geneRef" ]
        for gene_etree in OG_genes:
            if str(gene_etree.get("id")) in list_genes_int:
                if return_hog_etree:
                    gene_with_its_HOG[str(map_int_2_ext[gene_etree.get("id")])]=OG_etree
                else:
                    gene_with_its_HOG[str(map_int_2_ext[gene_etree.get("id")])]=OG_etree.get("id")
    return gene_with_its_HOG

def get_size_all_HOGs(XML_root):
    '''
    XML_root: root of the etree element
    return: dictionary key:hog_id and value:number of genes
    '''
    OGs = {}
    for OG_etree in get_all_hogs(XML_root):
        OG_genes = [ gene_etree for gene_etree in OG_etree.getiterator() if gene_etree.tag == "{http://orthoXML.org/2011/}geneRef" ]
        OGs[OG_etree.get("id")]= len(OG_genes)
    return OGs

def get_speciesTree(fn_orthoxml):
    '''
    Extract from a orthoxml the species tree based on the internal topology
    :param XML: orthoxml file
    :return: newick string
    '''
    op = fa.OrthoXMLParser(fn_orthoxml)
    tax = fa.TaxonomyFactory.newTaxonomy(op)
    op.augmentTaxonomyInfo(tax)
    return tax.newick()


def get_list_species(XML_root):
    '''
    return a list of species contains into an orthoxml.

    :XML_root: root of the etree element
    :return: list of species
    '''
    species = [ e.get("name") for e in XML_root.getiterator() if e.tag == "{http://orthoXML.org/2011/}species" ]
    return species
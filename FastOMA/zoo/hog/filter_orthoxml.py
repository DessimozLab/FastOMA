from ..utils import auto_open
import collections
from time import time
from lxml import etree as ET
import Bio.Phylo
from typing import Iterable
from pathlib import Path
import logging
logger = logging.getLogger(__name__)
class HOGFilter:
    def __init__(self, score:str, value:float):
        self.score = score
        self.value = value

    def remove(self, score_id, value):
        return score_id == self.score and self.value > float(value)


class OrthoXMLFilterProcesser:

    def __init__(self, filters:Iterable[HOGFilter]=None):
        self.filters = list(filters)

    def add_filter(self, filter:HOGFilter):
        self.filters.append(filter)

    def process(self, fh):
        NS = "http://orthoXML.org/2011/"
        self.doc = ET.parse(fh)
        root = self.doc.getroot()
        to_rem = []
        for hog in root.iterfind('.//{{{0}}}orthologGroup'.format(NS)):
            score = hog.find('./{{{0}}}score'.format(NS))
            if score is None:
                continue
                #return
            for filt in self.filters:
                if filt.remove(score.get('id'), score.get('value')):
                    to_rem.append(hog)
                    break
        logger.info(f"will remove {len(to_rem)} hogs")
        for h in to_rem:
            h.getparent().remove(h)


    def write(self, fh):
        self.doc.write(fh, xml_declaration=True, encoding="UTF-8")



def filter_orthoxml_file(source_orthoxml, out, filter: HOGFilter):
    processor = OrthoXMLFilterProcesser([filter])
    if isinstance(source_orthoxml, (str, bytes, Path)):
        with auto_open(source_orthoxml, 'rt') as fh:
            processor.process(fh)
    else:
        processor.process(source_orthoxml)
    processor.write(out)









# __author__ = 'admin'
#
#
# import FastOMA.zoo.familyanalyzer as fa
# # import lxml.etree as etree
# # from familyanalyzer.orthoxmlquery import ElementError, OrthoXMLQuery
#
#
#
#
# # To parse the orthoxml and get the related etree object: etree.parse(fn_orthoxml) denote here by XML
# # To get the root of your etree object: etree_object.getroot() denote here by XML_root
# # The orthoxml is using two ids:
# #   - The internal id denotes by the tag "id" which is found in the mapping in the beginning of the orthoxml and inside the orthologous groups
# #   - The external id denotes by the tag "protId"(sometimes with the tag "geneId" or BOTH) which is found only in the begining of the orthoxml
#
#
# def remove_useless_information(etree_element):
#     '''
#     This function cleaned the orthoxml file by removing all information related to genes/species that are not
#     belonging to any OGs
#     :param etree_element: etree.parse(fn_orthoxml)
#     :return: cleaned etree_element
#     '''
#
#     root_XML = etree_element.getroot()
#
#     useful_genes = []
#
#     # Look at all genes present in OGs
#     for gene_etree in root_XML.getiterator():
#         if gene_etree.tag == "{http://orthoXML.org/2011/}geneRef" :
#             useful_genes.append(str(gene_etree.get("id")))
#
#     # Get all genes present in the orthoxml
#     input_genes = fa.OrthoXMLQuery.getInputGenes(root_XML)
#
#     # Remove all genes that are not present in an OG
#     for gene_el in input_genes:
#         if str(gene_el.get("id")) not in useful_genes:
#             gene_el.getparent().remove(gene_el)
#
#
#     # Remove empty species (if all genes have been deleted)
#     genes_map = [ e for e in root_XML.getiterator() if e.tag == "{http://orthoXML.org/2011/}genes" ]
#     for gene in genes_map:
#         if gene.getchildren() == []:
#             db = gene.getparent()
#             db.getparent().remove(db)
#     sp_map = [ e for e in root_XML.getiterator() if e.tag == "{http://orthoXML.org/2011/}species" ]
#     for sp in sp_map:
#         if sp.getchildren() == []:
#             sp.getparent().remove(sp)
#
#     return etree_element
#
#
#
# def filter_hogs_by_id(XML, list_id, output):
#     '''
#     Filter out all toplevel OG of an orthoxml file that are not in the list_id and write it into a new orthoxml file.
#     :param XML: a orthoxml etree obj
#     :param list_id: the list of top level OG id we want to keep
#     :param output: fn for the filtered orthoxml
#     '''
#
#     root_XML = XML.getroot()
#     topLevel_OG = fa.OrthoXMLQuery.getToplevelOrthologGroups(root_XML)
#     for OG in topLevel_OG:
#         hog_id = OG.get("id")
#         if str(hog_id) not in list_id:
#             OG.getparent().remove(OG)
#     XML_cleaned = remove_useless_information(XML)
#     XML_cleaned.write(output)

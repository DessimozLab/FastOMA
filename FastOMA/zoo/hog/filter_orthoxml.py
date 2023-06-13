
from ..utils import auto_open
# import collections
# from time import time
from lxml import etree as ET
# import Bio.Phylo
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
            for filt in self.filters:
                if filt.remove(score.get('id'), score.get('value')):
                    to_rem.append(hog)
                    break
        logger.info(f"will remove {len(to_rem)} hogs")
        for h in to_rem:
            parent = h.getparent()
            if 'id' in h.attrib:
                logger.info("removing hog " + str(h) + " line " + str(h.sourceline) + " " +str(h.attrib['id']))
            else:
                logger.info("removing hog " + str(h) + " line " + str(h.sourceline))
            if parent:
                parent.remove(h)
                if sum(c.tag == "{{{0}}}orthologGroup".format(NS) for c in parent) == 0:
                    if 'id' in parent.attrib:
                        logger.info("consider deleting the empty parent hog "+str(parent)+" line "+str(parent.sourceline)+" "+str(parent.attrib['id']))
                    else:
                        logger.info("consider deleting the empty parent hog " + str(parent) + " line "+str(parent.sourceline))
                    to_rem.append(parent)

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





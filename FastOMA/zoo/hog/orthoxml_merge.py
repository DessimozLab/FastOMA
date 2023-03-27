from xml.etree import ElementTree as ET
from typing import List, Iterable
from random import randint



class GeneRefManager:
    def __init__(self):
        self.xrefs = {}
        self.ids = set([])

    def _random_unused_id(self):
        while True:
            cand = randint(100000, 1000000000)
            if str(cand) not in self.ids:
                return str(cand)

    def register_and_reassign(self, gene_nodes:Iterable[ET.Element]):
        update_ids = {}
        to_rem = []
        for gene in gene_nodes:
            if gene.attrib['id'] in self.ids:
                if gene.attrib['protId'] in self.xrefs:
                    # protId already in set. is it unique? if yes, no action, otherwise error
                    if self.xrefs[gene.attrib['protId']] != gene.attrib['id']:
                        raise ValueError("protId '{}' is used several times with different gene id :'{},'{}'"
                                         .format(gene.attrib['protId'], self.xrefs[gene.attrib['protId']], gene.attrib['id']))
                    else:
                        to_rem.append(gene.attrib['id'])
                        continue
                else:
                    # reassign internal gene id.
                    new_id = self._random_unused_id()
                    update_ids[gene.attrib['id']] = new_id
                    gene.attrib['id'] = new_id

            self.xrefs[gene.attrib['protId']] = gene.attrib['id']
            self.ids.add(gene.attrib['id'])
        return update_ids, to_rem


class Merger:
    def __init__(self, first):
        self.NS = "http://orthoXML.org/2011/"
        ET.register_namespace("", self.NS)
        self.doc = ET.parse(first)
        self.root = self.doc.getroot()

        self.all_species = set(z.attrib['name'] for z in self.doc.findall('./{{{}}}species'.format(self.NS)))
        self.all_genes = GeneRefManager()
        self.all_genes.register_and_reassign(
            self.doc.findall("./{{{0}}}species/{{{0}}}database/{{{0}}}genes/{{{0}}}gene".format(self.NS))
        )

    def merge_file(self, other):
        gene_id_updates, to_rem = self.all_genes.register_and_reassign(
            other.findall("./{{{0}}}species/{{{0}}}database/{{{0}}}genes/{{{0}}}gene".format(self.NS)))
        self._remove_unnecessary_genes(other, to_rem)
        self._update_geneRef_ids(other.find('./{{{}}}groups'.format(self.NS)), gene_id_updates)

        for sp in other.findall("./{{{}}}species".format(self.NS)):
            if sp.attrib['name'] not in self.all_species:
                species_seen = False
                for i, el in enumerate(self.root):
                    if el.tag == "{{{}}}species".format(self.NS):
                        species_seen = True
                    elif species_seen:
                        break
                self.root.insert(i, sp)
                self.all_species.add(sp.attrib['name'])
            else:
                db = self.root.find("./{{{0}}}species[@name='{1}']/{{{0}}}database/{{{0}}}genes".format(self.NS, sp.attrib['name']))
                for g in sp.iterfind(".//{{{}}}gene".format(self.NS)):
                    db.append(g)
        grps = self.root.find("./{{{}}}groups".format(self.NS))
        for g in other.find("./{{{}}}groups".format(self.NS)):
            grps.append(g)

    def _update_geneRef_ids(self, root, gene_id_updates):
        for old_id, new_id in gene_id_updates.items():
            for g in root.iterfind(".//{{{0}}}geneRef[@id='{1}']".format(self.NS, old_id)):
                g.attrib['id'] = new_id

    def _remove_unnecessary_genes(self, root, to_rem):
        for e in to_rem:
            parent = root.find("./{{{0}}}species/{{{0}}}database/{{{0}}}genes/{{{0}}}gene[@id='{1}']/.."
                               .format(self.NS, e))
            child = parent.find("./{{{0}}}gene[@id='{1}']".format(self.NS, e))
            parent.remove(child)




    def write(self, fh):
        self.doc.write(fh, xml_declaration=True, encoding="UTF-8", default_namespace=None)


def merge_orthoxml_files(out, files):
    """function to merge several orthoxml files into a single orthoxml file that contains all groups.

    This function combines several orthoxml files into a single orthoxml file that
    contains all the groups and maintains a valid definition block of the species
    and their genes. The protId attributes among all the orthoxml files need to be
    either unique or being at least assigned to the same internal gene id; in that
    case it is assumed that it is the same gene across the different files and it
    can be merged.
    if the gene id attribute is the same two or more orthoxml files, but their
    protId value is different, a new gene id value is generated and the geneRef
    values are updated accordingly.

    :param out: a path or a filehandle object where the combined orthoxml data should
                be written to.

    :param files: a list of paths or filehandle objects (of valid orthoxml format) that
                should be merged.

    """

    first = files.pop()
    merger = Merger(first)
    for f in files:
        merger.merge_file(ET.parse(f).getroot())

    return merger.write(out)

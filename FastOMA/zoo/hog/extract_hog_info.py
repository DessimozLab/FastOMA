from ..utils import auto_open
import collections
from time import time
import xml.etree.ElementTree as etree
from pathlib import Path
import logging
logger = logging.getLogger(__name__)

Gene = collections.namedtuple("Gene", "xref species internal_id")


class SpeciesAnalyser:
    def __init__(self, gene_attr="protId"):
        self.gene_attr = gene_attr
        self.genes = {}
        self.nr_genes_per_species = collections.defaultdict(int)

    def add_genome_genes(self, genome_node):
        genome_name = genome_node.get('name', None)
        if genome_name is None:
            genome_name = genome_node.get("NCBITaxId")

        generef_2_xref = {}
        for gene in genome_node.findall('.//{http://orthoXML.org/2011/}gene'):
            gene_id = gene.get('id')
            gene_prot_id = gene.get(self.gene_attr)
            generef_2_xref[gene_id] = Gene(gene_prot_id, genome_name, gene_id)
            self.nr_genes_per_species[genome_name] += 1
        self.genes.update(generef_2_xref)

    def gene_in_group(self, gene_id):
        self.genes.pop(gene_id)

    def get_singletons(self):
        return self.genes

    def summary(self):
        single = collections.defaultdict(int)
        for g in self.genes.values():
            single[g.species] += 1
        return [{'species': g, 'genes': self.nr_genes_per_species[g], 'not_in_group': single[g]}
                for g in self.nr_genes_per_species]


def parse_orthoxml(fh, genome_watcher: SpeciesAnalyser):
    taxonomy = {}
    og_level = 0

    def collect_genes(elem):
        genes = 0
        for child in elem.iter():
            if child == elem:
                continue
            if child.tag == "{http://orthoXML.org/2011/}geneRef":
                genes += 1
                if genome_watcher is not None:
                    genome_watcher.gene_in_group(child.get('id'))
            elif child.tag == "{http://orthoXML.org/2011/}orthologGroup":
                genes += child.text
        elem.clear()
        elem.text = genes
        return genes

    logger.info("start mapping of orthoxml formatted input file")
    for event, elem in etree.iterparse(fh, events=('start', 'end')):
        if event == "start":
            if elem.tag == "{http://orthoXML.org/2011/}orthoXML":
                if elem.get('version') != "0.5":
                    raise RuntimeError(f"Expecting orthoXML version 0.5, but is {elem.get('version')}")
            elif elem.tag == '{http://orthoXML.org/2011/}orthologGroup':
                og_level += 1
        elif event == 'end':
            if elem.tag == "{http://orthoXML.org/2011/}orthologGroup":
                og_level -= 1
                data = {'id': elem.get('id'), 'level': taxonomy[elem.get('taxonId')]}
                for child in elem.findall('./{http://orthoXML.org/2011/}score'):
                    data[child.get('id')] = float(child.get('value'))
                data['nr_members'] = collect_genes(elem)
                data['is_roothog'] = og_level == 0
                yield data
                if og_level == 0:
                    elem.clear()
            elif elem.tag == "{http://orthoXML.org/2011/}species":
                if genome_watcher is not None:
                    genome_watcher.add_genome_genes(elem)
                elem.clear()
            elif elem.tag == "{http://orthoXML.org/2011/}taxon":
                taxonomy[elem.get('id')] = elem.get('name')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthoxml", required=True)
    conf = parser.parse_args()
    genome_coverage_stats = SpeciesAnalyser()
    with open(conf.orthoxml, 'rt') as xml:
        for group in parse_orthoxml(xml, genome_coverage_stats):
            print(group)
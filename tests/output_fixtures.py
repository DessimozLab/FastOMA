"""Hand-made inputs for fastoma-collect-subhogs (the output stage).

Species tree (((A,B)AB,C)ABC,D)root. Proteins are named a1..a7, b1..b7, c1..c7, d1..d7;
the numeric gene id of e.g. b3 is 2003, its protId 'sp|P02003|B3_B'.

rootHOG E0000001 (nested duplications, like _hog_class.HOG.to_orthoxml produces them):

    og root
      paralogGroup
        og ABC
          paralogGroup
            og AB [a1 b1]
            og AB [a2 b2]
          c1
        og ABC
          paralogGroup
            og AB [a3 b3]
            og AB [a4 b4]
          c2
      d1

rootHOG E0000002 (root-level duplication, wrapped as in _infer_subhog.py:119-128, plus a
fragment paralogGroup with a Type property as in _hog_class.py:520-527):

    og AB
      paralogGroup
        og AB [a5 b5]
        og AB [paralogGroup(Type=Dubiousfragment)[a6 a7] b6]

Proteins a7 can optionally be marked as a minor isoform of a1. b7, c3.. etc. are in the
header but in no group.
"""
import pickle
import xml.etree.ElementTree as ET
from pathlib import Path

from FastOMA._utils_roothog import Gene

SPECIES_TREE = "(((A,B)AB,C)ABC,D)root;"
SPECIES_IDX = {"A": 1, "B": 2, "C": 3, "D": 4}
PROTEIN_NAMES = [f"{s}{i}" for s in "abcd" for i in range(1, 8)]
RHOG1_MEMBERS = ["a1", "a2", "a3", "a4", "b1", "b2", "b3", "b4", "c1", "c2", "d1"]
RHOG2_MEMBERS = ["a5", "a6", "a7", "b5", "b6"]


def species(name):
    return name[0].upper()


def gene_id(name):
    return SPECIES_IDX[species(name)] * 1000 + int(name[1:])


def prot_id(name):
    return f"sp|P{gene_id(name):05d}|{name.upper()}_{species(name)}"


def accession(name):
    return f"P{gene_id(name):05d}"


def sequence(name):
    return "MKV" + name.upper()


def _og(tax, *children, id_="x", score=True):
    e = ET.Element("orthologGroup", {"id": id_})
    if score:
        ET.SubElement(e, "score", {"id": "CompletenessScore", "value": "1.0"})
    ET.SubElement(e, "property", {"name": "TaxRange", "value": tax})
    for c in children:
        e.append(c)
    return e


def _pg(*children, prop=None):
    e = ET.Element("paralogGroup")
    if prop:
        ET.SubElement(e, "property", {"name": "Type", "value": prop})
    for c in children:
        e.append(c)
    return e


def _g(name):
    return ET.Element("geneRef", {"id": str(gene_id(name))})


def rhog1_hog(score=True):
    def og(*a, **kw):
        return _og(*a, score=score, **kw)
    return og("root",
              _pg(og("ABC", _pg(og("AB", _g("a1"), _g("b1")), og("AB", _g("a2"), _g("b2"))), _g("c1")),
                  og("ABC", _pg(og("AB", _g("a3"), _g("b3")), og("AB", _g("a4"), _g("b4"))), _g("c2"))),
              _g("d1"), id_="HOG_E0000001_sub10001")


def rhog2_hog(score=True):
    def og(*a, **kw):
        return _og(*a, score=score, **kw)
    return og("AB",
              _pg(og("AB", _g("a5"), _g("b5")),
                  og("AB", _pg(_g("a6"), _g("a7"), prop="Dubiousfragment"), _g("b6"))),
              id_="HOG_E0000002_sub10002")


def build_collect_inputs(folder, scores=True, isoforms=False):
    """Write the inputs of fastoma-collect-subhogs into `folder` and return their paths.

    pickles/rhog_<n>/file_HOG_E000000<n>.pickle, gene_id_dic_xml.pickle, species_tree.nwk,
    rhogs/HOG_E000000<n>.fa
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    tree = folder / "species_tree.nwk"
    tree.write_text(SPECIES_TREE)

    main_of = {"a7": "a1"} if isoforms else {}
    gene_map = {sp: [Gene(gene_id(n), prot_id(n), gene_id(main_of[n]) if n in main_of else None)
                     for n in PROTEIN_NAMES if species(n) == sp]
                for sp in "ABCD"}
    if isoforms:
        # the main isoform references itself (_utils_roothog.save_gene_id_mapping)
        gene_map["A"][0] = Gene(gene_id("a1"), prot_id("a1"), gene_id("a1"))
    gene_pickle = folder / "gene_id_dic_xml.pickle"
    with open(gene_pickle, "wb") as fh:
        pickle.dump(gene_map, fh)

    pickles = folder / "pickles"
    for nr, hog in ((1, rhog1_hog(scores)), (2, rhog2_hog(scores))):
        sub = pickles / f"rhog_{nr}"
        sub.mkdir(parents=True, exist_ok=True)
        with open(sub / f"file_HOG_E000000{nr}.pickle", "wb") as fh:
            pickle.dump([hog], fh)

    rhogs = folder / "rhogs"
    rhogs.mkdir(exist_ok=True)
    for fn, members in (("HOG_E0000001.fa", RHOG1_MEMBERS), ("HOG_E0000002.fa", RHOG2_MEMBERS)):
        with open(rhogs / fn, "w") as fh:
            for n in members:
                fh.write(f">{prot_id(n)}||{species(n)}||{gene_id(n)} desc of {n}\n{sequence(n)}\n")
    return {"folder": folder, "species_tree": tree, "gene_id_pickle": gene_pickle,
            "pickle_folder": pickles, "rhogs": rhogs}

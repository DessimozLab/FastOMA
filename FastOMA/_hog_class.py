
import xml.etree.ElementTree as ET
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from random import sample
from typing import Optional, List, Union
from ete3 import Tree, TreeNode
import random
from ._utils_subhog import MSAFilter

# todo some of these could also come under conf_infer_subhhogs
seed_random=1234 # Also in _wrappers.py
random.seed(seed_random)
fragment_detection = True  # this also need to be consistent in _infer_subhog.py
fragment_detection_msa_merge = True
subsampling_hogclass = True
hogclass_max_num_seq = 10              # 40 subsampling in msa # ver very 2
hogclass_min_cols_msa_to_filter = hogclass_max_num_seq * 50
#hogclass_tresh_ratio_gap_col = 0.05     # 0.8 for very very big


import itertools
from . import _utils_subhog, logger


def _member_species(members) -> set:
    return set(m.split("||")[1] for m in members)


def _species_name_index(tree: TreeNode) -> dict:
    """species name -> TreeNode lookup for `tree`, built once and cached on its root.

    ete3's search_nodes()/get_common_ancestor(names) each do a full linear scan of the tree;
    at scale (e.g. a rootHOG with 100k+ genes, each needing a per-species lookup) calling those
    once per gene dominates runtime. Building this index once per species tree and doing O(1)
    dict lookups instead is ~1000x faster in practice and is safe to cache on the tree: a fresh
    species tree object is built per rootHOG (see prepare_species_tree), so the cache can't go
    stale across rootHOGs, and it's built lazily on first use so callers never pass it explicitly."""
    root = tree.get_tree_root()
    index = getattr(root, "_species_name_index", None)
    if index is None:
        index = {n.name: n for n in root.traverse()}
        root.add_feature("_species_name_index", index)
    return index


def _species_mrca(tree: TreeNode, species_of_members: set) -> TreeNode:
    """Common ancestor of `species_of_members` within `tree`.

    `tree` should be the real, unpruned species tree whenever one is available: the working
    species (sub)tree used during inference is pruned to the species present in the rootHOG
    (see `_utils_subhog.prepare_species_tree`), so a species entirely absent from the rootHOG
    no longer exists in it, and computing losses against that tree would silently ignore it.
    A node that genuinely branches for `species_of_members` is never affected by that pruning
    (pruning only ever removes nodes that become unary), so this mrca has the same identity/name
    whether computed against the pruned or the unpruned tree -- only the topology below it, which
    ImpliedLosses/TCSScore need, differs."""
    index = _species_name_index(tree)
    return tree.get_common_ancestor(*[index[x] for x in species_of_members])


def _count_implied_losses(mrca: TreeNode, species_with_members: set) -> int:
    """Dollo-parsimony style count of minimal loss events between `mrca` and the species
    that actually have a member gene: a whole clade lacking any member counts as one loss,
    regardless of how many species it spans."""
    def recurse(node):
        leaves_under = set(n.name for n in node.iter_leaves())
        if leaves_under.isdisjoint(species_with_members):
            return 1
        if node.is_leaf():
            return 0
        return sum(recurse(c) for c in node.children)
    return sum(recurse(c) for c in mrca.children)


def _count_untouched_clades(tree_node: TreeNode, reached_nodes: list) -> int:
    """Counts, below `tree_node`, the maximal clades that contain none of `reached_nodes` --
    i.e. clades no sub-hog lineage touches at all -- as one loss each. Does not descend into a
    clade that a lineage already reaches: that clade's internal losses (if any) are counted
    separately, per lineage, by _hog_implied_losses -- re-scanning them here would double-count."""
    if tree_node in reached_nodes:
        return 0
    if not any(tree_node in r.get_ancestors() for r in reached_nodes):
        return 1
    return sum(_count_untouched_clades(c, reached_nodes) for c in tree_node.children)


def _hog_implied_losses(hog: "HOG", tree_node: TreeNode) -> int:
    """Recursively counts implied gene losses for `hog` within `tree_node`'s subtree.

    Mirrors _tax_overlap's _tax_now-based grouping of _subhogs, but sums losses per group member
    instead of intersecting/folding lineages: a group of sub-hogs sharing the same _tax_now
    represents a duplication, and each such paralogous copy's losses are counted independently
    so that one copy's presence in a species can't mask a sibling copy's loss there."""
    subhogs = getattr(hog, "_subhogs", None)
    if not subhogs:
        return 0

    losses = 0
    clade_groups = {}
    for sh in subhogs:
        clade_groups.setdefault(sh._tax_now.name, []).append(sh)
    for subhogs_of_clade in clade_groups.values():
        clade_node = subhogs_of_clade[0]._tax_now
        for sh in subhogs_of_clade:
            losses += _hog_implied_losses(sh, clade_node)

    reached_nodes = [subhogs_of_clade[0]._tax_now for subhogs_of_clade in clade_groups.values()]
    losses += _count_untouched_clades(tree_node, reached_nodes)
    return losses


def _species_lineage_index(tree: TreeNode) -> dict:
    """species name -> frozenset of ancestor node names (inclusive of the species itself and the
    tree's root), cached on the tree's root like _species_name_index.

    Unlike _species_mrca (which needs a lineage relative to some already-known mrca), the TCS
    taxonomy-overlap score below wants each species' *absolute* lineage: it discovers the HOG's
    own mrca as a side effect of intersecting members' lineages, rather than needing it supplied
    up front (see attach_scores)."""
    root = tree.get_tree_root()
    index = getattr(root, "_species_lineage_index", None)
    if index is None:
        index = {}
        for leaf in root.iter_leaves():
            ancestors = set()
            n = leaf
            while True:
                ancestors.add(n.name)
                if n is root:
                    break
                n = n.up
            index[leaf.name] = frozenset(ancestors)
        root.add_feature("_species_lineage_index", index)
    return index


def _combine_tax_overlap(parts):
    """Folds sibling (nset, leaf_size, leaf_acc, tax_score) tuples into their parent's, per
    Moi/Kim's taxonomy-overlap algorithm: a "match" (nonempty lineage shared by every
    contributing part) earns len(nset) points per gene (leaf_size) below this node, on top of
    whatever each part already scored deeper down. A part with an empty nset (no informative
    match anywhere in its own subtree) is excluded from the intersection rather than zeroing it
    out -- one already-incongruent branch shouldn't poison a genuinely-matching sibling."""
    leaf_size = sum(p[1] for p in parts)
    leaf_acc = sum(p[2] for p in parts) + leaf_size
    tax_score = sum(p[3] for p in parts)
    nonempty = [p[0] for p in parts if p[0]]
    nset = frozenset.intersection(*nonempty) if nonempty else frozenset()
    tax_score += len(nset) * leaf_size
    return nset, leaf_size, leaf_acc, tax_score


def _tax_overlap(hog: "HOG", species_lineage_index: dict):
    """Recursively computes (nset, leaf_size, leaf_acc, tax_score) for `hog`, mirroring the same
    _tax_now-based grouping of _subhogs used by HOG.to_orthoxml() (a group of sub-hogs sharing a
    _tax_now is a duplication -- its members are folded together as siblings, same as a group of
    literal gene-tree children would be). `leaf_size` is read off len(hog).
    See attach_scores for how the four returned values become TCSScore."""
    if not hog._subhogs:
        species = next(iter(_member_species(hog.get_members())))
        return species_lineage_index[species], len(hog), 0, 0

    groups = []
    for _, subhogs_of_clade in itertools.groupby(
            sorted(hog._subhogs, key=lambda h: h._tax_now.name), key=lambda h: h._tax_now.name):
        subhogs_of_clade = list(subhogs_of_clade)
        if len(subhogs_of_clade) == 1:
            groups.append(_tax_overlap(subhogs_of_clade[0], species_lineage_index))
        else:
            groups.append(_combine_tax_overlap([_tax_overlap(sh, species_lineage_index) for sh in subhogs_of_clade]))
    return groups[0] if len(groups) == 1 else _combine_tax_overlap(groups)


def attach_scores(hog_element: ET.Element, hog: "HOG", mrca: TreeNode, species_of_members: set) -> None:
    """Computes and attaches CompletenessScore, TCSScore and ImpliedLosses as <score>
    sub-elements of `hog_element`.

    TCSScore follows Moi et al. 2025 / Kim et al. 2026's taxonomy-overlap score: _tax_overlap()
    computes it using each species' absolute lineage (not one relative to `mrca`), and the HOG's
    own mrca-relative "ideal" contribution is discovered and subtracted algebraically at the end
    (tax_score - leaf_acc * len(nset)) rather than needing `mrca` supplied up front, then
    normalized by leaf_size -- the gene count of `hog` itself."""
    completeness_score = round(len(species_of_members) / mrca.size, 4)
    ET.SubElement(hog_element, "score", attrib={"id": "CompletenessScore", "value": str(completeness_score)})

    if getattr(hog, "_subhogs", None):
        implied_losses = _hog_implied_losses(hog, mrca)
    else:
        implied_losses = _count_implied_losses(mrca, species_of_members)
    ET.SubElement(hog_element, "score", attrib={"id": "ImpliedLosses", "value": str(implied_losses)})

    nset, leaf_size, leaf_acc, tax_score = _tax_overlap(hog, _species_lineage_index(mrca))
    tcs_score = (tax_score - leaf_acc * len(nset)) / leaf_size
    ET.SubElement(hog_element, "score", attrib={"id": "TCSScore", "value": str(round(tcs_score, 4))})


# from .infer_subhogs import conf_infer_subhhogs #fastoma_infer_subhogs #


class Representative:
    def __init__(self, the_one, elements:Optional[List]=None):
        self.enabled = True
        if isinstance(the_one, Representative):
            self._seq = the_one.get_record()
            self._subelements = the_one.get_subelements()
            self._species = the_one.get_species()
        elif isinstance(the_one, SeqRecord):
            self._seq = the_one
            self._subelements = set([the_one.id])
            # TODO: this should be done more genenerally. copied from singleton_hog function
            self._species = {the_one.id.split('||')[1]}
        if elements is not None:
            for e in elements:
                self._subelements |= e.get_subelements()
                self._species |= e.get_species()

    def get_id(self):
        return self._seq.id

    def get_record(self):
        return self._seq

    def get_subelements(self):
        return self._subelements

    def get_species(self):
        return self._species

    def __repr__(self):
        return '<Representative: {} [{} elements; {} species]>'.format(self.get_id(), len(self.get_subelements()), len(self.get_species()))

class HOG:
    _hogid_iter = 10000

    def __init__(self, input_instantiate, taxnomic_range:TreeNode, rhogid:str, msa:Optional[MultipleSeqAlignment] = None,
                 representatives=Optional[List[Representative]], conf_infer_subhhogs=None):
        """HOG class constructor.

        A HOG can be initialized with either a single SeqRecord or from a list of HOGs
        (which will become its sub-hog elements). Furthermore, the taxonomic range
        (a node from the species tree) needs to be provided.
        Parameters:
        :param input_instantiate"""
        # fragment_list list of  sets , each set contains protein ID of fragments
        # the input_instantiate could be either
        #     1) orthoxml_to_newick.py protein as the biopython seq record  SeqRecord(seq=Seq('MAPSSRSPSPRT. ]
        # or  2) orthoxml_to_newick.py set of intances of class HOG   wit orthoxml_to_newick.py big msa
        # those variable starting with _ are local to the class, should not access directly  (although it is possbile)
        self._rhogid = rhogid
        self.__class__._hogid_iter += 1
        # 0070124

        self._hogid = "HOG_" + self._rhogid+ "_sub" + str(self.__class__._hogid_iter)
        self._tax_now = taxnomic_range    # the taxnomic level that we are considering now, checking for duplication
        self.active = True

        if isinstance(input_instantiate, SeqRecord):  # if len(sub_hogs)==1:
            only_protein = input_instantiate  # only one seq, only on child, leaf
            self._members = set([only_protein.id])
            self._representatives = [Representative(only_protein)]
            self._msa = MultipleSeqAlignment([only_protein])
            self._subhogs = []
            self._dubious_members = set()
            # <<class 'Bio.Align.MultipleSeqAlignment'> instance (1 records of length 314) at 7f0e86c713d0>

        elif msa and all(isinstance(x, HOG) for x in input_instantiate):
            # here we want to merge few subhogs and create orthoxml_to_newick.py new HOG.   #the n

            sub_hogs = input_instantiate
            hog_members = set()
            for sub_hog in sub_hogs:
                hog_members |= sub_hog.get_members()
            self._members = hog_members  # set.union(*tup)    # this also include dubious_members
            self._subhogs = list(input_instantiate)  # full members of subhog, children
            self._representatives = [r for r in representatives]
            representatives_id = set(r.get_id() for r in representatives)
            assert len(representatives_id) == len(self._representatives), "Non-unique ids among representatives"
            assert representatives_id.issubset(self._members), "Representatives must be subset of HOG members"

            dubious_members = set()
            for sub_hog in sub_hogs:
                dubious_members |= sub_hog.get_dubious_members()  # union
            self._dubious_members = dubious_members
            assert self._dubious_members.isdisjoint(representatives_id), \
                   'Dubious members has non-empty intersection with representatives:' + str(self._dubious_members.intersection(representatives_id))

            records = [record for record in msa if (record.id in representatives_id)]

            if len(records[0]) > hogclass_min_cols_msa_to_filter and conf_infer_subhhogs is not None:
                filter = MSAFilter(levelprocessor=None, conf=conf_infer_subhhogs)
                records = filter.msa_filter_col(records)
                # the challange is that one of the sequences might be complete gap
            self._msa = MultipleSeqAlignment(records)
            # without replacement sampling ,  # self._children = sub_hogs # as legacy  ?
        else:
            logger.error("Error 142769,  check the input format to instantiate HOG class")
            assert False

    def __repr__(self):
        return f"<HOG:{self.hogid},size={len(self._members)},tax={self._tax_now.name}>"

    @property
    def hogid(self):
        return self._hogid

    @property
    def taxname(self):
        return self._tax_now.name

    @property
    def taxlevel(self):
        return self._tax_now

    @property
    def rhogid(self):
        return self._rhogid

    def __len__(self):
        return len(self._members)

    def get_members(self):
        return set(self._members)

    def get_dubious_members(self):
        return self._dubious_members

    def get_representatives(self):
        return self._representatives

    def get_msa(self):
        return self._msa

    def find_representative(self, seq_id):
        """Find the representative of the given sequence using its ID."""
        for r in self._representatives:
            if r.get_id() == seq_id:
                return r
        for r in self._representatives:
            if seq_id in r.get_members():
                return r
        raise KeyError("No representative found for " + seq_id)

    def remove_prot_from_hog(self, prot_to_remove):
        prot_members_hog_old = self._members
        assert prot_members_hog_old
        prot_members_hog_edited = prot_members_hog_old - set([prot_to_remove])
        self._members = prot_members_hog_edited # self._members.remove(prot_to_remove)    # discard
        msa_old = self._msa     #  we may want to edit the msa of children level to be consistent
        msa_edited = MultipleSeqAlignment([i for i in msa_old if i.id != prot_to_remove])
        self._msa = msa_edited
        if len(prot_members_hog_edited) == 0:  # hog should be removed, no members is left
            return 0
        return 1

    def insert_dubious_prots(self, fragment_host, fragments_list_nothost):

        self._dubious_members |= set(fragments_list_nothost)
        self._dubious_members |= {fragment_host}
        if fragment_host not in self._members:
            logger.error("Error 1252769, fragment_host not in hog"+str(fragment_host))
        self._members |= set(fragments_list_nothost)

        msa_old = self._msa
        msa_edited = MultipleSeqAlignment([i for i in msa_old if i.id not in self._dubious_members])
        self._msa = msa_edited
        # species_name = fragment_host.split("||")[1]   # if self._tax_now == species_name:
        return 1

    def __contains__(self, item):
        if isinstance(item, Representative):
            return item.get_id() in self._members
        elif isinstance(item, HOG):
            return item in self._subhogs
        else:
            return item in self._members

    def get_subhog_path(self, seq_id, max_depth=-1):
        if seq_id not in self._members:
            raise KeyError(f"{seq_id} is not part of this HOG")
        path = []
        h = self
        while max_depth != 0 and len(h._subhogs) > 0:
            for s in h._subhogs:
                if seq_id in s:
                    path.append(s)
                    h = s
                    break
            max_depth -= 1
        return path

    def merge_prots_name_hog(self, fragment_name_host, merged_fragment_name):
        prot_members_hog = list(self._members)
        assert fragment_name_host in prot_members_hog, str(fragment_name_host)+str("  prot_members_hog:")+str(prot_members_hog)
        fragment_host_idx = prot_members_hog.index(fragment_name_host)
        # merged_fragment_name = fragment_name_host + "_|_" + fragment_name_remove1 + ... # there might be more than one fragment
        prot_members_hog[fragment_host_idx] = merged_fragment_name
        # 'BUPERY_R03529||BUPERY||1105002086_|_BUPERY_R10933||BUPERY||1105008975']
        self._members = set(prot_members_hog)

        return 1

    def merge_prots_msa(self, merged_fragment_name, merged_msa_new):   # merged_fragment_name 'BUPERY_R03529||BUPERY||1105002086_|_BUPERY_R10933||BUPERY||1105008975']
        prot_members_hog = list(self._members)
        assert merged_fragment_name in prot_members_hog  # merged_fragment_name in self._members   -> True . already the name of host prot is changed using merge_prots_name_hierarchy_toleaves
        msa_new = []
        for seq_record in merged_msa_new:
            if seq_record.id in prot_members_hog:
                msa_new.append(seq_record)
            # else:  # some of the sequences in the msa is for other subhogs, so we only care about those in this subhog.
        # msa_old = self._msa
        # host_fragment_name = merged_fragment_name.split("_|_")[0]
        # merged_msa_new_ids = [i.id for i in merged_msa_new]         # for seq_record in msa_old:
        self._msa = MultipleSeqAlignment(msa_new)

        return 1


    # def merge_prots_msa(self, fragment_name_host, fragment_name_remove, merged_sequence, merged_msa_new):
    #     merged_fragment_name = fragment_name_host + "_|_" + fragment_name_remove
    #     prot_members_hog = list(self._members)
    #     assert merged_fragment_name in prot_members_hog
    #     # merged_fragment_name in self._members   -> True . already the name of host prot is changed using merge_prots_name_hierarchy_toleaves
    #     msa_new = []
    #     for prot_name in self._members:
    #         if seq_rec.id == fragment_name_host:
    #             seq_rec_edited = SeqRecord(Seq(merged_sequence), id=merged_fragment_name, name=merged_fragment_name)
    #         msa_new.append(seq_rec_edited)
    #     self._msa = MultipleSeqAlignment(msa_new)
    #     return 1

    def to_orthoxml(self, full_species_tree: Optional[TreeNode] = None):
        if len(self._subhogs) == 0:
            list_member = list(self._members)
            if len(list_member) == 1:
                list_member_first = list(self._members)[0]  # ['tr|A0A3Q2UIK0|A0A3Q2UIK0_CHICK||CHICK_||1053007703']
                if fragment_detection and fragment_detection_msa_merge and "_|_" in list_member_first:
                    paralog_element = ET.Element('paralogGroup')
                    #property_element = ET.SubElement(paralog_element, "property", attrib={"name": "TaxRange", "value": str(self._tax_now)})
                    property_element = ET.SubElement(paralog_element, "property", attrib={"name": "Type", "value": "DubiousMergedfragment"})
                    #todo is it happening in output?
                    # property_element = ET.SubElement(paralog_element, "property", attrib={"TaxRange": str(self._tax_now), "Type": "DubiousMergedfragment"})
                    list_member_first_fragments = list_member_first.split("_|_")
                    for fragment in list_member_first_fragments:
                        prot_name_integer = fragment.split("||")[2].strip()
                        geneRef_elemnt = ET.Element('geneRef', attrib={'id': str(prot_name_integer)})
                        paralog_element.append(geneRef_elemnt)
                    return paralog_element
                else:
                    prot_name_integer = list_member_first.split("||")[2].strip()
                    geneRef_elemnt = ET.Element('geneRef', attrib={'id': str(prot_name_integer)})
                    #'id': str(gene_id_name[list_member_first])})  # # gene_id_name[query_prot_record.id]   # hog_elemnt.append(geneRef_elemnt)
                    # to do could be improved when the rhog contains only one protein
                    return geneRef_elemnt

            elif len(list_member) > 1 and fragment_detection and (not fragment_detection_msa_merge):
                paralog_element = ET.Element('paralogGroup')
                property_element = ET.SubElement(paralog_element, "property", attrib={"name": "Type", "value":"Dubiousfragment"})
                # property_element = ET.SubElement(paralog_element, "property", attrib={"TaxRange": str(self._tax_now),"Type": "DubiousMergedfragment"})
                # removed proteins which are dubious are reported only in the log files.
                for member in list_member:
                    prot_name_integer = member.split("||")[2].strip()
                    geneRef_elemnt = ET.Element('geneRef', attrib={'id': str(prot_name_integer)})
                    paralog_element.append(geneRef_elemnt)
                return paralog_element

        # We continue this function as an Implicit else :  len(self._subhogs) >=1

        def _sorter_key(sh):
            return sh._tax_now.name
        self._subhogs.sort(key=_sorter_key)

        element_list = []
        for sub_clade, sub_hogs in itertools.groupby(self._subhogs, key=_sorter_key):  # sub_clade is the taxrange
            list_of_subhogs_of_same_clade = list(sub_hogs)
            # list_of_subhogs_of_same_clade is  [object of HOGclass hogID=HOG_D0634402_sub10094,len=12, tax_least=Amniota, tax_now= Amniota, object of HOGclass hogID=HOG_D0634402_sub10096,len=20, tax_least=Amniota, tax_now= Amniota, object of HOGclass hogID=HOG_D0634402_sub10081,len=1, tax_least=MONDO, tax_now= Amniota, object of HOGclass hogID=HOG_D0634402_sub10093,len=12, tax_least=Amniota, tax_now= Amniota, object of HOGclass hogID=HOG_D0634402_sub10095,len=7, tax_least=Amniota, tax_now= Amniota]
            # following only for debugging, can be deleted later
            for subhog in list_of_subhogs_of_same_clade:
                if len(subhog._members) == 0:
                    logger.warning("issue 12314506" + str(subhog) + str(sub_clade))

            if len(list_of_subhogs_of_same_clade) > 1:
                paralog_element = ET.Element('paralogGroup')
                # todo add score how much overlap they have, report as score
                # the following line could be improved, instead of tax_now we can use the least common ancestor of all members
                # property_element = ET.SubElement(paralog_element, "property",attrib={"name": "TaxRange", "value": str(sub_clade)}) # self._tax_now
                for sh in list_of_subhogs_of_same_clade:
                    element_p = sh.to_orthoxml(full_species_tree)
                    if str(element_p):
                        paralog_element.append(element_p)  # ,**gene_id_name  indent+2
                    else:
                        logger.warning("issue 123434" + str(sh) + str(sub_clade))

                element_list.append(paralog_element)
            elif len(list_of_subhogs_of_same_clade) == 1:
                subhog = list_of_subhogs_of_same_clade[0]
                if len(subhog._members):
                    element = subhog.to_orthoxml(full_species_tree)
                    if str(element):  # element could be  <Element 'geneRef' at 0x7f7f9bacb450>
                        element_list.append(element)  # indent+2
                    else:
                        logger.warning("issue 12359 "+str(subhog))
                else:
                    logger.warning("issue 12369 " + str(subhog))

        if len(element_list) == 1:
            return element_list[0]

        elif len(element_list) > 1:
            #hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(self._hogid)})
            hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(self._hogid)}, )
            species_of_members = _member_species(self._members)  #  'tr|H2MU14|H2MU14_ORYLA||ORYLA||1056022282'
            mrca = _species_mrca(full_species_tree if full_species_tree is not None else self.taxlevel, species_of_members)
            if mrca.name != self.taxlevel.name:
                logger.info(f"mrca ({mrca.name}) != self.taxlevel ({self.taxlevel.name})")
                logger.info(f"<{hog_elemnt.tag} {hog_elemnt.attrib}>")

            attach_scores(hog_elemnt, self, mrca, species_of_members)
            property_element = ET.SubElement(hog_elemnt, "property", attrib={"name": "TaxRange", "value": str(mrca.name)})

            for element in element_list:
                hog_elemnt.append(element)

            return hog_elemnt

""" 
        Example of orthoxml format
        <orthologGroup id="HOG_0022401_sub10169">
           <property name="TaxRange" value="Alcidae"/>
           <paralogGroup>
              <geneRef id="1163015305"/>
              <orthologGroup id="HOG_0022401_sub10166">
                 <property name="TaxRange" value="Alca_torda-Uria"/>
                 <orthologGroup id="HOG_0022401_sub10164">
                    <property name="TaxRange" value="Uria"/>
                    <geneRef id="1164015054"/>
                    <geneRef id="1180013314"/>
                 </orthologGroup>
              </orthologGroup>
           </paralogGroup>
        </orthologGroup>
"""


def split_hog(hog:HOG, level_name:str, *partitions):
    """splits a hog into parts based to the partitions provided.

    The partitions need to be a lists of Representatives (or simply members)
    of the given HOG.

    :param partitions:  a list or a tuple of Representatives of the current HOG
    :type partitions:   List[Representatives]
    :returns List[HOG]: a list of HOG objects splitting the original HOG in the
                        desired partitions.
    """
    max_depth = 4
    if len(partitions) < 2:
        raise ValueError(f"Must provide at least two partitions to split {hog}")
    partitions = [set(hog.find_representative(r) for r in part) for part in partitions]
    if not all(p1.isdisjoint(p2) for p1, p2 in itertools.combinations(partitions, 2)):
        raise ValueError(f"Not all partitions are non-overlapping {partitions}")
    #if sum(len(p) for p in partitions) != len(set(hog.get_representatives())):
    #    raise ValueError(f"Partitions to split {hog} must cover all {len(hog.get_representatives())} representatives, but cover only {sum(len(p) for p in partitions)}.")
    rep_to_subhog_list = {rep: hog.get_subhog_path(rep.get_id(), max_depth) for part in partitions for rep in part}

    logger.debug(f"Subhog paths of represenatatives for {hog}")
    for i, part in enumerate(partitions):
        logger.debug(f"---Partition {i} -----")
        for rep in part:
            path = " --> ".join(str(h) for h in rep_to_subhog_list[rep])
            logger.debug(f"{rep.get_id()}: {path}")

    first_level_subhogs = [set(rep_to_subhog_list[rep][0] for rep in part) for part in partitions]
    if all(s1.isdisjoint(s2) for s1, s2 in itertools.combinations(first_level_subhogs, 2)):
        logger.info(f"Splitting {hog} into {len(first_level_subhogs)} subhogs: {first_level_subhogs}")
        hogs = []
        for p, subhogs in enumerate(first_level_subhogs):
            h = HOG(subhogs, taxnomic_range=hog.taxlevel, rhogid=hog.rhogid, msa=hog.get_msa(), representatives=partitions[p])
            hogs.append(h)
        return hogs
    else:
        rep_to_subhog_list = {rep: hog.get_subhog_path(rep.get_id(), max_depth=-1) for part in partitions for rep in part}
        for depth in range(max(len(sh) for sh in rep_to_subhog_list.values())):
            rep_sets = [set(rep_to_subhog_list[rep][depth] for rep in part if len(rep_to_subhog_list[rep]) > depth) for part in partitions]
            if all(s1.isdisjoint(s2) for s1, s2 in itertools.combinations(rep_sets, 2)):
                break
        depth -= 1
        rep_sets = [set(rep_to_subhog_list[rep][depth] for rep in part if len(rep_to_subhog_list[rep]) > depth) for part in partitions]
        merged_in = set.union(*(s1.intersection(s2) for s1, s2 in itertools.combinations(rep_sets, 2)))
        involved_reps = [(r, h[depth].hogid, h[depth].taxname) for r, h in rep_to_subhog_list.items() if len(h)>depth and h[depth] in merged_in]
        logger.warning(f"{hog} should be split in {len(partitions)} partitions at {level_name}. Not implemented yet.")
        for rep in involved_reps:
            logger.warning(f"  - Rep {rep[0]} merged into {rep[1]} at level {rep[2]}")
    #raise RuntimeError("this part of the code needs more thinking")


import xml.etree.ElementTree as ET
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from random import sample
import itertools

from . import _utils_subhog
from ._utils_subhog import logger_hog
from . import _config


class HOG:
    _hogid_iter = 10000

    def __init__(self, input_instantiate, taxnomic_range, rhogid, msa=None, num_species_tax_speciestree=None):
        # fragment_list list of  sets , each set contains protein ID of fragments
        # the input_instantiate could be either
        #     1) orthoxml_to_newick.py protein as the biopython seq record  SeqRecord(seq=Seq('MAPSSRSPSPRT. ]
        # or  2) orthoxml_to_newick.py set of intances of class HOG   wit orthoxml_to_newick.py big msa
        # those variable starting with _ are local to the class, should not access directly  (although it is possbile)
        self._rhogid = rhogid
        self.__class__._hogid_iter += 1
        # 0070124
        #todo add release id,
        self._hogid = "HOG_" + self._rhogid+ "_sub" + str(self.__class__._hogid_iter)
        self._tax_least = taxnomic_range  #  least taxnomic level
        self._tax_now = taxnomic_range    # the taxnomic level that we are considering now, checking for duplication


        if isinstance(input_instantiate, SeqRecord):  # if len(sub_hogs)==1:
            only_protein = input_instantiate  # only one seq, only on child, leaf
            self._members = set([only_protein.id])
            self._msa = MultipleSeqAlignment([only_protein])
            self._subhogs = []
            self._dubious_members = set()
            self._num_species_tax_speciestree = 1  # at the leaf level, there is only one species in the species tree
            # <<class 'Bio.Align.MultipleSeqAlignment'> instance (1 records of length 314) at 7f0e86c713d0>

        elif msa and all(isinstance(x, HOG) for x in input_instantiate):
            # here we want to merge few subhogs and create orthoxml_to_newick.py new HOG.   #the n

            sub_hogs = input_instantiate
            hog_members = set()
            for sub_hog in sub_hogs:
                hog_members |= sub_hog.get_members()  # union
            self._members = hog_members  # set.union(*tup)    # this also include dubious_members
            self._subhogs = list(input_instantiate)  # full members of subhog, children
            self._num_species_tax_speciestree = num_species_tax_speciestree

            dubious_members = set()
            for sub_hog in sub_hogs:
                dubious_members |= sub_hog.get_dubious_members()  # union
            self._dubious_members = dubious_members

            records_full = [record for record in msa if (record.id in self._members) and (record.id not in self._dubious_members) ]

            if len(records_full[0]) > _config.hogclass_min_cols_msa_to_filter:
                records_sub_filt = _utils_subhog.msa_filter_col(records_full, _config.hogclass_tresh_ratio_gap_col)
                # the challange is that one of the sequences might be complete gap
            else:
                records_sub_filt = records_full  # or even for rows # msa_filt_row_col = _utils.msa_filter_row(msa_filt_row, tresh_ratio_gap_row)

            if _config.subsampling_hogclass and len(records_sub_filt) > _config.hogclass_max_num_seq:
                # to do in future:  select best seq, not easy to defin, keep diversity,
                records_sub_sampled_raw = sample(list(records_sub_filt), _config.hogclass_max_num_seq)  # without replacement.
                records_sub_sampled = _utils_subhog.msa_filter_col(records_sub_sampled_raw, 0.01) # to make sure no empty column
                logger_hog.info("we are doing subsamping in hog class from " + str(len(records_full)) + " to " + str(_config.hogclass_max_num_seq) + " seqs.")
            else:
                records_sub_sampled = records_sub_filt
                # removing some columns completely gap - (not x   )
            self._msa = MultipleSeqAlignment(records_sub_sampled)
            # without replacement sampling ,  # self._children = sub_hogs # as legacy  ?
        else:
            logger_hog.error("Error 142769,  check the input format to instantiate HOG class")
            assert False

    def __repr__(self):
        return "HOGobj:" + self._hogid + ",size=" + str(
            len(self._members))+", taxLeast=" + str(self._tax_least) + ", taxNow= " + str(self._tax_now)

    def get_members(self):
        return set(self._members)

    def get_dubious_members(self):
        return self._dubious_members

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
            logger_hog.error("Error 1252769, fragment_host not in hog"+str(fragment_host))
        self._members |= set(fragments_list_nothost)

        msa_old = self._msa
        msa_edited = MultipleSeqAlignment([i for i in msa_old if i.id not in self._dubious_members])
        self._msa = msa_edited
        # species_name = fragment_host.split("||")[1]   # if self._tax_now == species_name:
        return 1


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



    def to_orthoxml(self):

        if len(self._subhogs) == 0:
            list_member = list(self._members)
            if len(list_member) == 1:
                list_member_first = list(self._members)[0]  # ['tr|A0A3Q2UIK0|A0A3Q2UIK0_CHICK||CHICK_||1053007703']
                if _config.fragment_detection and _config.fragment_detection_msa_merge and "_|_" in list_member_first:
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

            elif len(list_member) > 1 and _config.fragment_detection and (not _config.fragment_detection_msa_merge):
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
            return sh._tax_now
        self._subhogs.sort(key=_sorter_key)

        element_list = []
        for sub_clade, sub_hogs in itertools.groupby(self._subhogs, key=_sorter_key):  # sub_clade is the taxrange
            list_of_subhogs_of_same_clade = list(sub_hogs)
            # list_of_subhogs_of_same_clade is  [object of HOGclass hogID=HOG_D0634402_sub10094,len=12, tax_least=Amniota, tax_now= Amniota, object of HOGclass hogID=HOG_D0634402_sub10096,len=20, tax_least=Amniota, tax_now= Amniota, object of HOGclass hogID=HOG_D0634402_sub10081,len=1, tax_least=MONDO, tax_now= Amniota, object of HOGclass hogID=HOG_D0634402_sub10093,len=12, tax_least=Amniota, tax_now= Amniota, object of HOGclass hogID=HOG_D0634402_sub10095,len=7, tax_least=Amniota, tax_now= Amniota]
            # following only for debugging, can be deleted later
            for subhog in list_of_subhogs_of_same_clade:
                if len(subhog._members) == 0:
                    logger_hog.warning("issue 12314506" + str(subhog) + str(sub_clade))

            if len(list_of_subhogs_of_same_clade) > 1:
                paralog_element = ET.Element('paralogGroup')
                # todo add score how much overlap they have, report as score
                # the following line could be improved, instead of tax_now we can use the least common ancestor of all members
                # property_element = ET.SubElement(paralog_element, "property",attrib={"name": "TaxRange", "value": str(sub_clade)}) # self._tax_now
                for sh in list_of_subhogs_of_same_clade:
                    element_p = sh.to_orthoxml()
                    if str(element_p):
                        paralog_element.append(element_p)  # ,**gene_id_name  indent+2
                    else:
                        logger_hog.warning("issue 123434" + str(sh) + str(sub_clade))

                element_list.append(paralog_element)
            elif len(list_of_subhogs_of_same_clade) == 1:
                subhog = list_of_subhogs_of_same_clade[0]
                if len(subhog._members):
                    element = subhog.to_orthoxml()
                    if str(element):  # element could be  <Element 'geneRef' at 0x7f7f9bacb450>
                        element_list.append(element)  # indent+2
                    else:
                        logger_hog.warning("issue 12359 "+str(subhog))
                else:
                    logger_hog.warning("issue 12369 " + str(subhog))

        if len(element_list) == 1:
            return element_list[0]

        elif len(element_list) > 1:
            #hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(self._hogid)})
            hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(self._hogid)}, )
            num_species_tax_hog = len(set([i.split("||")[1] for i in self._members]))  #  'tr|H2MU14|H2MU14_ORYLA||ORYLA||1056022282'
            completeness_score = round(num_species_tax_hog/self._num_species_tax_speciestree,4)
            property_element = ET.SubElement(hog_elemnt, "score", attrib={"id": "CompletenessScore", "value": str(completeness_score)})
            property_element = ET.SubElement(hog_elemnt, "property", attrib={"name": "TaxRange", "value": str(self._tax_now)})

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
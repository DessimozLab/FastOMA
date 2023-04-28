
import xml.etree.ElementTree as ET
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from random import sample
import itertools
from . import _utils_subhog
from ._utils_subhog import logger_hog
from . import _config
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class HOG:
    _hogid_iter = 10000

    def __init__(self, input_instantiate, taxnomic_range, rhogid_num, msa=None):
        # fragment_list list of  sets , each set contains protein ID of fragments
        # the input_instantiate could be either
        #     1) orthoxml_to_newick.py protein as the biopython seq record  SeqRecord(seq=Seq('MAPSSRSPSPRT. ]
        # or  2) orthoxml_to_newick.py set of intances of class HOG   wit orthoxml_to_newick.py big msa
        # those variable starting with _ are local to the class, should not access directly  (although it is possbile)
        self._rhogid_num = rhogid_num
        self.__class__._hogid_iter += 1
        # 0070124
        self._hogid = "HOG_" + str(self._rhogid_num).zfill(7) + "_sub" + str(self.__class__._hogid_iter)
        self._tax_least = taxnomic_range  #  least taxnomic level
        self._tax_now = taxnomic_range    # the taxnomic level that we are considering now, checking for duplication

        if isinstance(input_instantiate, SeqRecord):  # if len(sub_hogs)==1:
            only_protein = input_instantiate  # only one seq, only on child, leaf
            self._members = set([only_protein.id])
            self._msa = MultipleSeqAlignment([only_protein])
            self._subhogs = []
            self._dubious_members = set()
            # <<class 'Bio.Align.MultipleSeqAlignment'> instance (1 records of length 314) at 7f0e86c713d0>

        elif msa and all(isinstance(x, HOG) for x in input_instantiate):
            # here we want to merge few subhogs and create orthoxml_to_newick.py new HOG.   #the n

            sub_hogs = input_instantiate
            hog_members = set()
            for sub_hog in sub_hogs:
                hog_members |= sub_hog.get_members()  # union
            self._members = hog_members  # set.union(*tup)    # this also include dubious_members
            self._subhogs = list(input_instantiate)  # full members of subhog, children

            dubious_members = set()
            for sub_hog in sub_hogs:
                dubious_members |= sub_hog.get_dubious_members()  # union
            self._dubious_members = dubious_members

            records_full = [record for record in msa if (record.id in self._members) and (record.id not in self._dubious_members) ]
            if len(records_full) > _config.hogclass_max_num_seq:
                # to do in future:  select best seq, not easy to defin, keep diversity,
                records_sub_sampled_raw = sample(records_full, _config.hogclass_max_num_seq)  # without replacement.

                if _config.subsampling_hogclass:
                    if len(records_sub_sampled_raw[0]) > _config.hogclass_min_cols_msa_to_filter:
                        records_sub_sampled = _utils_subhog.msa_filter_col(records_sub_sampled_raw, _config.hogclass_tresh_ratio_gap_col)
                    else:
                        records_sub_sampled = records_sub_sampled_raw
                    # or even for rows # msa_filt_row_col = _utils.msa_filter_row(msa_filt_row, tresh_ratio_gap_row)
                    logger_hog.info( "we are doing subsamping in hig class from " + str(len(records_full)) + " to " + str(_config.hogclass_max_num_seq) + " seqs.")
                else:
                    records_sub_sampled = records_sub_sampled_raw

            else:
                records_sub_sampled = records_full
                # removing some columns completely gap - (not x   )
            self._msa = MultipleSeqAlignment(records_sub_sampled)
            # without replacement sampling ,  # self._children = sub_hogs # as legacy  ?
        else:
            logger_hog.error("Error 142769,  check the input format to instantiate HOG class")
            assert False

    def __repr__(self):
        return "object of HOGclass hogID=" + self._hogid + ",len=" + str(
            len(self._members))+", tax_least=" + str(self._tax_least) + ", tax_now= " + str(self._tax_now)

    def get_members(self):
        return set(self._members)

    def get_dubious_members(self):
        return self._dubious_members

    def remove_prots_from_hog(self, prots_to_remove):
        prot_members_hog_old = self._members
        if len(prot_members_hog_old) & len(prots_to_remove):
            prot_members_hog_edited = prot_members_hog_old - set(prots_to_remove)
            self._members = prot_members_hog_edited
            msa_old = self._msa
            # we may not to edit the msa of children level
            msa_edited = MultipleSeqAlignment([i for i in msa_old if i.id not in prots_to_remove])
            self._msa = msa_edited
            if len(prot_members_hog_edited) == 0:  # hog should be removed , no members is left
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

        # species_name = fragment_host.split("||")[1]
        # if self._tax_now == species_name:
        return 1


    def merge_prots_name_hog(self, fragment_name_host, fragment_name_remove):

        prot_members_hog = list(self._members)
        assert fragment_name_host in prot_members_hog, str(fragment_name_host)+str("  prot_members_hog:")+str(prot_members_hog)
        fragment_host_idx = prot_members_hog.index(fragment_name_host)
        merged_fragment_name = fragment_name_host + "_|_" + fragment_name_remove
        prot_members_hog[fragment_host_idx] = merged_fragment_name
        # 'BUPERY_R03529||BUPERY||1105002086_|_BUPERY_R10933||BUPERY||1105008975']
        self._members = prot_members_hog

        return 1



    # hog_host.merge_prots_msa(fragment_name_host, fragment_name_remove, merged_sequence, merged_msa_new)
    def merge_prots_msa(self,fragment_name_host, merged_fragment_name, merged_msa_new):
        # merged_fragment_name 'BUPERY_R03529||BUPERY||1105002086_|_BUPERY_R10933||BUPERY||1105008975']
        prot_members_hog = list(self._members)
        assert merged_fragment_name in prot_members_hog
        # merged_fragment_name in self._members   -> True . already the name of host prot is changed using merge_prots_name_hierarchy_toleaves
        # msa_old = self._msa
        msa_new = []
        for seq_record in merged_msa_new:
            if seq_record.id in  prot_members_hog:
                msa_new.append(seq_record)
            else:
                logger_hog.debug("issue 132851"+seq_record.id +"not in merged_msa_new ")

        self._msa = MultipleSeqAlignment(msa_new)

        return 1

    # # hog_host.merge_prots_msa(fragment_name_host, fragment_name_remove, merged_sequence, merged_msa_new)
    # def merge_prots_msa(self, fragment_name_host, fragment_name_remove, merged_sequence, merged_msa_new):
    #     merged_fragment_name = fragment_name_host + "_|_" + fragment_name_remove
    #     # 'BUPERY_R03529||BUPERY||1105002086_|_BUPERY_R10933||BUPERY||1105008975']
    #     prot_members_hog = list(self._members)
    #     assert merged_fragment_name in prot_members_hog
    #     # merged_fragment_name in self._members   -> True . already the name of host prot is changed using merge_prots_name_hierarchy_toleaves
    #     # msa_old = self._msa
    #     msa_new = []
    #     for prot_name in self._members:
    #
    #         if seq_rec.id == fragment_name_host:
    #             seq_rec_edited = SeqRecord(Seq(merged_sequence), id=merged_fragment_name, name=merged_fragment_name)
    #         msa_new.append(seq_rec_edited)
    #     self._msa = MultipleSeqAlignment(msa_new)
    #
    #     return 1


    def to_orthoxml(self):
        hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(self._hogid)})
        property_element = ET.SubElement(hog_elemnt, "property",
                                        attrib={"name": "TaxRange", "value": str(self._tax_now)})

        # todo double check we report the taxanomic level

        # to do the following could be improved ???   without this if it will be like, one property is enough
        # <orthologGroup>
        #    <property name="TaxRange" value="GORGO_HUMAN_PANTR"/>
        #    <property name="TaxRange" value="GORGO_HUMAN_PANTR"/>
        # if property_element not in hog_elemnt:
        #    hog_elemnt.append(property_element)
        #    print("*")
        # gene = ET.SubElement(species, "gene", attrib={"id":str(gene_counter), "protId":query_prot_record.id})
        # hog_elemnt = ET.SubElement(species,

        if len(self._subhogs) == 0:
            list_member = list(self._members)

            if len(list_member) == 1:
                list_member_first = list(self._members)[0]  # ['tr|A0A3Q2UIK0|A0A3Q2UIK0_CHICK||CHICK_||1053007703']
                if _config.merge_fragments_detected and  "_|_" in  list_member_first :
                    paralog_element = ET.Element('paralogGroup')
                    property_element = ET.SubElement(paralog_element, "property", attrib={"name": "Type", "value": "DubiousMergedfragment"})
                    list_member_first_fragments = list_member_first.split("_|_")
                    for fragment in list_member_first_fragments:
                        prot_name_integer = fragment.split("||")[2].strip()
                        geneRef_elemnt = ET.Element('geneRef', attrib={'id': str(prot_name_integer)})
                        paralog_element.append(geneRef_elemnt)
                    return paralog_element
                else:
                    prot_name_integer = list_member_first.split("||")[2].strip()
                    geneRef_elemnt = ET.Element('geneRef', attrib={'id': str(prot_name_integer)})
                    #'id': str(gene_id_name[list_member_first])})  # # gene_id_name[query_prot_record.id]
                # hog_elemnt.append(geneRef_elemnt)
                # to do could be improved when the rhog contains only one protein
                return geneRef_elemnt

            elif len(list_member) > 1 and not _config.merge_fragments_detected:
                # todo  a better flag is needed probably becuase of inserting dubious prots
                paralog_element = ET.Element('paralogGroup')
                property_element = ET.SubElement(paralog_element, "property", attrib={"name": "Type", "value":"Dubiousfragment"})
                # removed proteins which are dubious are reported only in the log files.
                for member in list_member:
                    prot_name_integer = member.split("||")[2].strip()
                    geneRef_elemnt = ET.Element('geneRef', attrib={'id': str(prot_name_integer)})
                    paralog_element.append(geneRef_elemnt)
            else:
                print("** issue 455922")

            return paralog_element  # hog_elemnt


        def _sorter_key(sh):
            return sh._tax_now
            # todo for checking whether it is paralogous group we check the level we are looking at. Not the tax_least (could be species level).

        self._subhogs.sort(key=_sorter_key)
        for sub_clade, sub_hogs in itertools.groupby(self._subhogs, key=_sorter_key):
            list_of_subhogs_of_same_clade = list(sub_hogs)
            # print(f'{" "*(indent+1)} clade: {sub_clade} with {str(len(list_of_subhogs_of_same_clade))}')
            if len(list_of_subhogs_of_same_clade) > 1:
                paralog_element = ET.Element('paralogGroup')
                for sh in list_of_subhogs_of_same_clade:
                    paralog_element.append(sh.to_orthoxml())  # ,**gene_id_name  indent+2
                hog_elemnt.append(paralog_element)
            else:
                hog_elemnt.append(list_of_subhogs_of_same_clade[0].to_orthoxml())  # indent+2

        return hog_elemnt


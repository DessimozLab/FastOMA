
import xml.etree.ElementTree as ET
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from random import sample
import itertools
# from . import _utils_subhog
from ._utils_subhog import logger_hog
from . import _config


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

                # todo not sure to do filtering for columns in hog class , may be problamtic with fragment detection
                records_sub_sampled = records_sub_sampled_raw
                # if len(records_sub_sampled_raw[0]) > _config.hogclass_min_cols_msa_to_filter:
                #     records_sub_sampled = _utils_subhog.msa_filter_col(records_sub_sampled_raw, _config.hogclass_tresh_ratio_gap_col)
                # else:
                #     records_sub_sampled = records_sub_sampled_raw
                ## or even for rows
                #         # msa_filt_row_col = _utils.msa_filter_row(msa_filt_row, tresh_ratio_gap_row)
                # logger_hog.info( "we are doing subsamping in hig class from " + str(len(records_full)) + " to " + str(max_num_seq) + " seqs.")
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




    def to_orthoxml(self):
        hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(self._hogid)})
        property_element = ET.SubElement(hog_elemnt, "property",
                                        attrib={"name": "TaxRange", "value": str(self._tax_now)})  # todo double check we report the taxanomic level

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
            list_member_first = list(self._members)[0]
            # 'tr|A0A3Q2UIK0|A0A3Q2UIK0_CHICK||CHICK_||1053007703'
            prot_name_integer = list_member_first.split("||")[2].strip()
            geneRef_elemnt = ET.Element('geneRef', attrib={'id': str(prot_name_integer)})
                #'id': str(gene_id_name[list_member_first])})  # # gene_id_name[query_prot_record.id]
            # hog_elemnt.append(geneRef_elemnt)
            # to do could be improved when the rhog contains only one protein
            return geneRef_elemnt  # hog_elemnt

        def _sorter_key(sh):
            return sh._tax_now  # for checking whether it is paralogous group we check the level we are looking at. Not the tax_least (could be species level).

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


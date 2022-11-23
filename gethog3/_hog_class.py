
import xml.etree.ElementTree as ET
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from random import sample
import itertools

import _utils
from _utils import logger_hog
import _config


class HOG:
    _hogid_iter = 10000

    def __init__(self, input_instantiate, taxnomic_range, rhogid_num, msa=None):  # _prot_names
        # the input_instantiate could be either
        #     1) orthoxml_to_newick.py protein as the biopython seq record  SeqRecord(seq=Seq('MAPSSRSPSPRT. ]
        # or  2) orthoxml_to_newick.py set of intances of class HOG   wit orthoxml_to_newick.py big msa
        # those variable starting with _ are local to the class, should not access directly  (although it is possbile)
        self._rhogid_num = rhogid_num
        self.__class__._hogid_iter += 1
        # 0070124
        self._hogid = "HOG:B" + str(self._rhogid_num).zfill(7) + "_sub" + str(self.__class__._hogid_iter)
        self._taxnomic_range = taxnomic_range  # print("**** orthoxml_to_newick.py new HOG is instantiated with id", self._hogid)

        if isinstance(input_instantiate, SeqRecord):  # if len(sub_hogs)==1:
            only_protein = input_instantiate  # only one seq, only on child, leaf
            self._members = set([only_protein.id])
            self._msa = MultipleSeqAlignment([only_protein])
            self._subhogs = []
            # <<class 'Bio.Align.MultipleSeqAlignment'> instance (1 records of length 314) at 7f0e86c713d0>

        elif msa and all(isinstance(x, HOG) for x in input_instantiate):
            # here we want to merge few subhogs and creat orthoxml_to_newick.py new HOG.   #the n
            sub_hogs = input_instantiate
            hog_members = set()
            for sub_hog in sub_hogs:
                hog_members |= sub_hog.get_members()  # union
            self._members = hog_members  # set.union(*tup)
            self._subhogs = list(input_instantiate)  # full members

            records_full = [record for record in msa if record.id in self._members]

            if len(records_full) > _config.hogclass_max_num_seq:
                records_sub_sampled_raw = sample(records_full, _config.hogclass_max_num_seq)  # without replacement.
                # to do in future:  select best seq

                if len(records_sub_sampled_raw[0]) > _config.hogclass_min_cols_msa_to_filter:
                    records_sub_sampled = _utils.msa_filter_col(records_sub_sampled_raw, _config.hogclass_tresh_ratio_gap_col)
                else:
                    records_sub_sampled = records_sub_sampled_raw
            #         # msa_filt_row_col = _utils.msa_filter_row(msa_filt_row, tresh_ratio_gap_row)
            #         # print(len(msa_filt_row_col), len(msa_filt_row_col[0]))
            # logger_hog.info( "we are doing subsamping now from " + str(len(records_full)) + " to " + str(max_num_seq) + " seqs.")
            else:
                records_sub_sampled = records_full
                # removing some columns completely gap -  (not x   )
            self._msa = MultipleSeqAlignment(records_sub_sampled)
            # without replacement sampling ,  # self._children = sub_hogs # as legacy  ?
        else:
            logger_hog.error("Error 169,  check the input format to instantiate orthoxml_to_newick.py HOG class")
            assert False

    def __repr__(self):
        return "an object of class HOG of hogID=" + self._hogid + ", length=" + str(
            len(self._members)) + ", taxonomy= " + str(self._taxnomic_range)

    def get_members(self):
        return set(self._members)

    def to_orthoxml(self):

        hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(self._hogid)})
        property_element = ET.SubElement(hog_elemnt, "property",
                                        attrib={"name": "TaxRange", "value": str(self._taxnomic_range)})
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
            return sh._taxnomic_range

        self._subhogs.sort(key=_sorter_key)
        for sub_clade, sub_hogs in itertools.groupby(self._subhogs, key=_sorter_key):
            list_of_subhogs_of_same_clade = list(sub_hogs)
            # print(f'{" "*(indent+1)} clade: {sub_clade} with {str(len(list_of_subhogs_of_same_clade))}')
            if len(list_of_subhogs_of_same_clade) > 1:
                paralog_element = ET.Element('paralogGroup')
                for sh in list_of_subhogs_of_same_clade:
                    # paralog_element.append(sh.to_orthoxml(**gene_id_name))
                    paralog_element.append(sh.to_orthoxml())  # ,**gene_id_name  indent+2
                hog_elemnt.append(paralog_element)
            else:
                hog_elemnt.append(list_of_subhogs_of_same_clade[0].to_orthoxml())  # indent+2
        return hog_elemnt
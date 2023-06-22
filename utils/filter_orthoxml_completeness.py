

import logging
logging.basicConfig(level=logging.DEBUG)


print("started 1 ")


#python  /work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/FastOMA/utils/filter_orthoxml_completeness.py  output_hog_.orthoxml 0.3  CompletenessScore

from FastOMA.zoo.hog import filter_orthoxml_file, HOGFilter

# folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/run_1june/"
# input_orthoxml_add = folder+"out_folder/test5"
#
# folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_qfo/run24may_xmlscore/"
# input_orthoxml_add = folder+"out_folder/output_hog_.orthoxml"
folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/run_1june/"

orthoxml_path=folder+"out_folder/output_hog_.orthoxml_scientific.orthoxml_edited"
input_orthoxml_add = orthoxml_path
threshold_filt = 0.3


#threshold_filt = 0.6
score_type = "CompletenessScore"


output_name = input_orthoxml_add + "_filt_"+str(threshold_filt)+".orthoxml"
with open(output_name, 'wb') as output_file:
    filt = HOGFilter(score_type, threshold_filt)
    filter_orthoxml_file(input_orthoxml_add, output_file, filt)

print("we wrote the output in "+output_name)


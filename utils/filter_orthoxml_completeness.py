
import sys
import io

import logging

logging.basicConfig(level=logging.DEBUG)

#python  /work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/FastOMA/utils/filter_orthoxml_completeness.py  output_hog_.orthoxml 0.3  CompletenessScore

from FastOMA.zoo.hog import filter_orthoxml_file, HOGFilter

input_orthoxml_add = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_qfo/run24may_xmlscore/out_folder/output_hog_.orthoxml"
#"output_hog_2.orthoxml" #sys.argv[1]
threshold_filt = 0.93
score_type = "CompletenessScore"

# if len(sys.argv) > 2:
#     threshold_filt = float(sys.argv[2])
# if len(sys.argv) > 3:
#     score_type = sys.argv[3]

# "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/FastOMA/archive//HOG_file_652949.pickle_12.orthoxml"


# res = io.BytesIO()

with open(input_orthoxml_add + "_filt2.orthoxml", 'wb') as output_file:
    filt = HOGFilter(score_type, threshold_filt)
    filter_orthoxml_file(input_orthoxml_add, output_file, filt)
    # out_orthoxml_str = res.getvalue()

print("we wrote the out in "+input_orthoxml_add + "_filt2.orthoxml")
#     for i in str(out_orthoxml_str).split("\\n"):
#         if i.startswith("b'<?"):
#             output_file.write(i[2:] + "\n")
#         else:
#             output_file.write(i + "\n")
#
# output_file.close()


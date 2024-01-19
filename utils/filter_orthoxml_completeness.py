


###  How to use:  python  filter_orthoxml_completeness.py  FastOMA_HOGs.orthoxml  0.3  

import sys
import logging
logging.basicConfig(level=logging.DEBUG)
from FastOMA.zoo.hog import filter_orthoxml_file, HOGFilter

print("started ")

input_orthoxml_add = sys.argv[1]
threshold_filt = float(sys.argv[2])

score_type = "CompletenessScore"


output_name = input_orthoxml_add + "_filt_"+str(threshold_filt)+".orthoxml"
with open(output_name, 'wb') as output_file:
    filt = HOGFilter(score_type, threshold_filt)
    filter_orthoxml_file(input_orthoxml_add, output_file, filt)

print("we wrote the output in "+output_name)



"""
I accidanetly comment <property na in hog class

so with this code I can edit the file


add       <orthologGroup id="HOG:B0811125_sub1201">
         <property name="TaxRange" value="CHLTR_MYCGE"/>


"""


file_in = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/xml_output/out_27aug_6pm.xml_no_property"
file_out = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/xml_output/out_27aug_6pm_property.xml"

file_in_handle = open(file_in, 'r')
file_out_handle = open(file_out, 'w')
property_str ="<property name=\"TaxRange\" value=\"test\"/>"
print("started")
for line in file_in_handle:
    if not "<orthologGroup" in line:
        file_out_handle.write(line)
    else:
        file_out_handle.write(line)

        for i, st in enumerate(line):
            if st == "<":
                needed_num = i

        file_out_handle.write(" "*(needed_num+2) + property_str+"\n")

print("finished")





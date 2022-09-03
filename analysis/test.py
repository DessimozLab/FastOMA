


from orthoxml_to_newick import orthoxml_to_newick


test_xml = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/xml_output/small_output.xml"

tree_dic = orthoxml_to_newick(test_xml)

for keys, values in tree_dic.items():
    print("\n")
    print(keys)
    print(values)
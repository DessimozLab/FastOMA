


# import OrthoXMLSplitter


folder="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/run_1june/out_folder/"
hog_file = folder + "/output_hog_.orthoxml"
outdir=folder+"/perrhog_folder"

from OrthoXMLSplitter import OrthoXMLSplitter

splitter = OrthoXMLSplitter(hog_file, outdir)

splitter()



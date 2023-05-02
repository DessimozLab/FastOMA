




import os
from os import listdir


folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_bird_/run_1may/out_folder/"



project_files = listdir(folder + "/rhogs_all/")
rhogs = []
for file in project_files:
    file_name_split = file.split(".")
    if file_name_split[-1] == "fa":
        rhog_id = int(file_name_split[0].split("_")[1])
        rhogs.append(rhog_id)

print("number of rhogs is ", len(rhogs))

folder_pickle = folder + "/pickle_rhogs/"
project_files = listdir(folder_pickle)
pickles = []
for file in project_files:
    if os.path.getsize(folder_pickle + file) > 2:
        file_name_split = file.split(".")
        if file_name_split[-1] == "pickle":
            rhog_id = int(file_name_split[0].split("_")[1])
            pickles.append(rhog_id)
    else:
        print("this file is empty", file)

print("number of pickles is ", len(pickles))

no_pickle_list = set(rhogs) - set(pickles)

print("number of rhogs not finished is ", len(no_pickle_list))

print("\n \n ", no_pickle_list)

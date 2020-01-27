# Python program to convert monthly model output to multi-year mean climat.
# The climate consistes of ANN, DJF, MAM, JJA, SON, and 01-12 months.
# NCO is used to do the average. 

import os
import numpy as np

#==========================================
# Python function to convert a list 
# to string using join() function 
def listToString(s):  
# initialize an empty string 
    str1 = ""  
# traverse in the string   
    for ele in s:  
        str1 += ele.rstrip()+" "   
# return string   
    return str1  

#==========================================
# The main program starts here.

# Input 
#caseid="E3SM_DECKv1b_H1.ne30"
#monthly_data_path="./E3SM_DECKv1b_H1.ne30/remap_180x360"
#caseid="E3SM_coupled_restart_20TR_Yr2000-emis_Scat.Year2000_2014"
caseid="E3SMv2_offline_noEmis"
print(caseid)
#monthly_data_path="/raid00/xianwen/E3SM_output/"+caseid+"/remap_180x360"
monthly_data_path="/global/cscratch1/sd/xianwen/acme_scratch/cori-knl/"+caseid+"/remap_180x360/"
#years=np.arange(0001,0004)
years=np.array(["2000","2001","2002"])
#years=np.array(["2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014"])
# Output
months_to_do=["01","02","03","04","05","06","07","08","09","10","11","12"]
seasons_to_do=["ANN","DJF","MAM","JJA","SON"]

out_path=monthly_data_path+"../climo"
if not os.path.exists(out_path):
     os.makedirs(out_path)

# create list of all input files
# Monthly climo
for mon in months_to_do:
    print("-- doing "+mon+"--")
    list_file="list_"+mon+".txt"
    if os.path.exists(list_file):
        os.system("rm "+list_file)
    for yr in years:
        os.system("ls "+monthly_data_path+"/*cam.h0*"+str(yr)+"-"+mon+".nc|cat >>"+list_file)
    #os.system("ls "+monthly_data_path+"/*"+mon+".nc >"+list_file)
    with open(list_file) as f_obj:
        lines=f_obj.readlines()
    lists=listToString(lines)
    
    climo_file=caseid+"_climo_"+mon+".nc"
    cmd="ncra "+lists+" "+out_path+"/"+climo_file
    os.system(cmd)
    os.system("mv "+list_file+" "+out_path+"/")
# Seasonal and Annual climo
for seasn in seasons_to_do:
    print("-- doing "+seasn+"--")
    list_file="list_"+seasn+".txt"
    if seasn == "ANN":
        mons_for_seasn=["01","02","03","04","05","06","07","08","09","10","11","12"]
        if os.path.exists(list_file):
            os.system("rm "+list_file)
        for mons in mons_for_seasn:
            os.system("ls "+out_path+"/*climo_"+mons+".nc|cat >>"+list_file)
    elif seasn == "DJF":
        mons_for_seasn=["12","01","02"]
        if os.path.exists(list_file):
            os.system("rm "+list_file)
        for mons in mons_for_seasn:
            os.system("ls "+out_path+"/*climo_"+mons+".nc|cat >>"+list_file)
    elif seasn == "MAM":
        mons_for_seasn=["03","04","05"]
        if os.path.exists(list_file):
            os.system("rm "+list_file)
        for mons in mons_for_seasn:
            os.system("ls "+out_path+"/*climo_"+mons+".nc|cat >>"+list_file)
    elif seasn == "JJA":
        mons_for_seasn=["06","07","08"]
        if os.path.exists(list_file):
            os.system("rm "+list_file)
        for mons in mons_for_seasn:
            os.system("ls "+out_path+"/*climo_"+mons+".nc|cat >>"+list_file)
    elif seasn == "SON":
        mons_for_seasn=["09","10","11"]
        if os.path.exists(list_file):
            os.system("rm "+list_file)
        for mons in mons_for_seasn:
            os.system("ls "+out_path+"/*climo_"+mons+".nc|cat >>"+list_file)

    with open(list_file) as f_obj:
        lines=f_obj.readlines()
    lists=listToString(lines)
    
    climo_file=caseid+"_climo_"+seasn+".nc"
    cmd="ncra "+lists+" "+out_path+"/"+climo_file
    os.system(cmd)
    os.system("mv "+list_file+" "+out_path+"/")
print("----- all finished -----")
	

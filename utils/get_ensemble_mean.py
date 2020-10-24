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
caseid=["CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl",\
        "CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl-ens1",\
        "CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl-ens2",\
        "CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl-ens3"]
#monthly_data_path="/glade/u/home/xianwen/work/archive/" #+caseid+"/atm/hist/"
monthly_data_path="/global/cscratch1/sd/xianwen/E3SM_simulations/" #+caseid+"/archive/climo_2005-2011/"
years=np.arange(2005,2011)
months_to_do=np.array(["01","02","03","04","05","06","07","08","09","10","11","12","ANN","JJA","DJF","MAM","SON"])
caseout="CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl-ens_mean_2005-2011"
out_path=monthly_data_path+caseout+"/"
if not os.path.exists(out_path):
     os.makedirs(out_path)
# create list of all input files
# Monthly ensemble mean
#for yr in years:
for mon in months_to_do:
    print("-- doing "+mon+" --")
    list_file="list_"+mon+".txt"
    if os.path.exists(list_file):
        os.system("rm "+list_file)
    for case in caseid:
       os.system("ls "+monthly_data_path+case+"/archive/climo_2005-2011/*climo*"+mon+".nc|cat >>"+list_file)

    with open(list_file) as f_obj:
        lines=f_obj.readlines()
    lists=listToString(lines)
    
    climo_file=caseout+"_climo_"+mon+".nc"
    cmd="ncra "+lists+" "+out_path+"/"+climo_file
    os.system(cmd)
    os.system("mv "+list_file+" "+out_path+"/")
print("----- all finished -----")
	

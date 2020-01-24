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
caseid="E3SMv2_standard_PresSST_UMRadALLoff"
print(caseid)
monthly_data_path="/global/cscratch1/sd/xianwen/acme_scratch/cori-knl/"+caseid+"/remap_180x360/"
#years=np.arange(0001,0004)
#years=np.array(["2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010"])
years=np.array(["0001","0002","0003","0004","0005","0006","0007","0008","0009","0010","0011"])
# Output
months_to_do=["01","02","03","04","05","06","07","08","09","10","11","12"]
#seasons_to_do=["ANN","DJF","MAM","JJA","SON"]
variables_to_ext=["TS","PSL","TREFHT","QREFHT","T","Q","Z3"]

out_path=monthly_data_path+"../extract/"
if not os.path.exists(out_path):
     os.makedirs(out_path)

#exract a variable from input files one by one
for var in variables_to_ext:
   list_file="list_"+var+".txt"
   if os.path.exists(list_file):
       os.system("rm "+list_file)
   for yr in years:
      for mon in months_to_do:
          print("-- doing "+var+" "+yr+" "+mon+"--")
          fin=monthly_data_path+caseid+".cam.h0."+yr+"-"+mon+".nc"
          fout=out_path+var+"_"+yr+"-"+mon+".nc"
          cmd="ncks -v "+var+" "+fin+" "+fout
          os.system(cmd)
          os.system("ls "+fout+"|cat >>"+list_file)          
#combine the extracted files for a specific variable
   with open(list_file) as f_obj:
       lines=f_obj.readlines()
   lists=listToString(lines)
   fext=out_path+var+"_2000-2010.nc"
   cmd="ncrcat "+lists+" "+fext 
   os.system(cmd)
   os.system("rm "+lists)
   os.system("mv "+list_file+" "+out_path)

print("----- all finished -----")
	

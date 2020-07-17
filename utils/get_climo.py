# Python program to convert monthly model output to multi-year mean climat.
# The climate consistes of ANN, DJF, MAM, JJA, SON, and 01-12 months.
# NCO is used to do the average. 

import os
import numpy as np
import glob

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
#caseid="E3SMv2_standard_PresSST_UMRadALLoff"
caseid="CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl"
print(caseid)
#monthly_data_path="/raid00/xianwen/E3SM_output/"+caseid+"/remap_180x360"
#monthly_data_path="/global/cscratch1/sd/xianwen/acme_scratch/cori-knl/"+caseid+"/remap_180x360/"
monthly_data_path="/global/cscratch1/sd/xianwen/E3SM_simulations/"+caseid+"/archive/remap_180x360/"
#years=np.arange(2000,2013)
years=np.array(["2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012"])
# Output
months_to_do=["01","02","03","04","05","06","07","08","09","10","11","12"]
seasons_to_do=["ANN","DJF","MAM","JJA","SON"]
do_monthly=True
do_seasonal=True

out_path=monthly_data_path+"../climo/"

# check input and output path
if not os.path.exists(monthly_data_path):
   print('PATH NOT FOUND: ', monthly_data_path)
   exit()
if not os.path.exists(out_path):
   os.makedirs(out_path)

# check existence of input files
for yr in years:
    for mon in months_to_do:
       file_tmp=glob.glob(monthly_data_path+"*cam.h0*"+str(yr)+"-"+mon+".nc")
       if len(file_tmp) != 1:
           print('Either path or file name prefix is incorrect:', \
                   monthly_data_path+"*cam.h0*"+str(yr)+"-"+mon+".nc")
           exit()
       file_now=file_tmp[0]
       if not os.path.exists(file_now):
           print("File not found!!! ",file_now)
           exit()
print('--- All input files are found ^_^ ---')

# check existance of old output files
for mon in months_to_do:
    climo_file=out_path+caseid+"_climo_"+mon+".nc"
    if os.path.exists(climo_file):
        print('Old file exists!!! '+climo_file)
        exit()
for seasn in seasons_to_do:
    climo_file=out_path+caseid+"_climo_"+seasn+".nc"
    if os.path.exists(climo_file):
        print('Ole file exists!!! '+climo_file)
        exit()
print('--- Output directory is clean ^_^ ---')

# comput monthly climo
if do_monthly:
   for mon in months_to_do:
     # create list of all input files
       print("-- doing "+mon+"--")
       list_file="list_"+mon+".txt"
       if os.path.exists(list_file):
           os.system("rm "+list_file)
       for yr in years:
           file_now=glob.glob(monthly_data_path+"*cam.h0*"+str(yr)+"-"+mon+".nc")[0]
           os.system("ls "+file_now+"|cat >>"+list_file)
       with open(list_file) as f_obj:
           lines=f_obj.readlines()
       lists=listToString(lines)

       # compute climo with the listed files above
       climo_file=caseid+"_climo_"+mon+".nc"
       cmd="ncra "+lists+" "+out_path+climo_file
       os.system(cmd)
       os.system("mv "+list_file+" "+out_path)

# compute annual and seasonal climo
if do_seasonal:
   for seasn in seasons_to_do:
     # creat list of all input files
       print("-- doing "+seasn+"--")
       list_file="list_"+seasn+".txt"
       if seasn == "ANN":
           mons_for_seasn=["01","02","03","04","05","06","07","08","09","10","11","12"]
           if os.path.exists(list_file):
               os.system("rm "+list_file)
           for mons in mons_for_seasn:
               file_now=glob.glob(monthly_data_path+"*cam.h0*"+str(yr)+"-"+mon+".nc")[0]
               os.system("ls "+file_now+"|cat >>"+list_file)
       elif seasn == "DJF":
           mons_for_seasn=["12","01","02"]
           if os.path.exists(list_file):
               os.system("rm "+list_file)
           for mons in mons_for_seasn:
               file_now=glob.glob(monthly_data_path+"*cam.h0*"+str(yr)+"-"+mon+".nc")[0]
               os.system("ls "+file_now+"|cat >>"+list_file)
       elif seasn == "MAM":
           mons_for_seasn=["03","04","05"]
           if os.path.exists(list_file):
               os.system("rm "+list_file)
           for mons in mons_for_seasn:
               file_now=glob.glob(monthly_data_path+"*cam.h0*"+str(yr)+"-"+mon+".nc")[0]
               os.system("ls "+file_now+"|cat >>"+list_file)
       elif seasn == "JJA":
           mons_for_seasn=["06","07","08"]
           if os.path.exists(list_file):
               os.system("rm "+list_file)
           for mons in mons_for_seasn:
               file_now=glob.glob(monthly_data_path+"*cam.h0*"+str(yr)+"-"+mon+".nc")[0]
               os.system("ls "+file_now+"|cat >>"+list_file)
       elif seasn == "SON":
           mons_for_seasn=["09","10","11"]
           if os.path.exists(list_file):
               os.system("rm "+list_file)
           for mons in mons_for_seasn:
               file_now=glob.glob(monthly_data_path+"*cam.h0*"+str(yr)+"-"+mon+".nc")[0]
               os.system("ls "+file_now+"|cat >>"+list_file)
   
       with open(list_file) as f_obj:
           lines=f_obj.readlines()
       lists=listToString(lines)
        
       # compute climo with the listed files above
       climo_file=caseid+"_climo_"+seasn+".nc"
       cmd="ncra "+lists+" "+out_path+climo_file
       os.system(cmd)
       os.system("mv "+list_file+" "+out_path)
print("----- all finished -----")
	

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
#caseid="E3SM_DECKv1b_H1.ne30"
#monthly_data_path="./E3SM_DECKv1b_H1.ne30/remap_180x360"
caseid="CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl"
print(caseid)
monthly_data_path="/global/cscratch1/sd/xianwen/E3SM_simulations/"+caseid+"/archive/remap_180x360/"
years=np.arange(2005,2013)
print(years)
months_all=["01","02","03","04","05","06","07","08","09","10","11","12"]
seasons_to_do=["ANN","DJF","MAM","JJA","SON"]

# creat output path if not existent:
out_path=monthly_data_path+"../climo/yby/"
if not os.path.exists(out_path):
     os.makedirs(out_path)

# check existence of input files
for yr in years:
    for mon in months_all:
        file_now=monthly_data_path+caseid+".cam.h0."+str(yr)+"-"+mon+".nc"
        if not os.path.exists(file_now):
            print("File not found!!! ",file_now)
            exit()
print('--- All input files are found ^_^ ---')

# check existence of old output file:
for yr in years:
    for seasn in seasons_to_do:
        climo_file=out_path+caseid+"_"+seasn+"_"+str(yr)+".nc"
        if os.path.exists(climo_file):
            print('Old file exists!!! ',climo_file)
            exit()
print('--- Output directory is clean ^_^ ---')

# compute seasonal and annual mean for each year
for seasn in seasons_to_do:
    print("-- doing "+seasn+" for each year --")
    if seasn == "ANN":
        mons_for_seasn=["01","02","03","04","05","06","07","08","09","10","11","12"]
        for yr in years:
          # create list of input files
            list_file="list_"+seasn+"_"+str(yr)+".txt"
            if os.path.exists(list_file):
                os.system("rm "+list_file)
            for mons in mons_for_seasn:
                os.system("ls "+monthly_data_path+"*"+str(yr)+"-"+mons+".nc|cat >>"+list_file)

            with open(list_file) as f_obj:
                lines=f_obj.readlines()
            lists=listToString(lines)

          # compute means from the listed files above    
            climo_file=caseid+"_"+seasn+"_"+str(yr)+".nc"
            cmd="ncra "+lists+" "+out_path+climo_file
            os.system(cmd)
            os.system("mv "+list_file+" "+out_path)

    elif seasn == "DJF":
        mons_for_seasn=["12","01","02"]
        for yr in years:
          # create list of input files
            list_file="list_"+seasn+"_"+str(yr)+".txt"
            if os.path.exists(list_file):
                os.system("rm "+list_file)
            for mons in mons_for_seasn:
                os.system("ls "+monthly_data_path+"*"+str(yr)+"-"+mons+".nc|cat >>"+list_file)

            with open(list_file) as f_obj:
                lines=f_obj.readlines()
            lists=listToString(lines)
           
          # compute means from the listed files above    
            climo_file=caseid+"_"+seasn+"_"+str(yr)+".nc"
            cmd="ncra "+lists+" "+out_path+climo_file
            os.system(cmd)
            os.system("mv "+list_file+" "+out_path)

    elif seasn == "MAM":
        mons_for_seasn=["03","04","05"]
        for yr in years:
          # create list of input files
            list_file="list_"+seasn+"_"+str(yr)+".txt"
            if os.path.exists(list_file):
                os.system("rm "+list_file)
            for mons in mons_for_seasn:
                os.system("ls "+monthly_data_path+"*"+str(yr)+"-"+mons+".nc|cat >>"+list_file)

            with open(list_file) as f_obj:
                lines=f_obj.readlines()
            lists=listToString(lines)
           
          # compute means from the listed files above    
            climo_file=caseid+"_"+seasn+"_"+str(yr)+".nc"
            cmd="ncra "+lists+" "+out_path+climo_file
            os.system(cmd)
            os.system("mv "+list_file+" "+out_path)

    elif seasn == "JJA":
        mons_for_seasn=["06","07","08"]
        for yr in years:
          # create list of input files
            list_file="list_"+seasn+"_"+str(yr)+".txt"
            if os.path.exists(list_file):
                os.system("rm "+list_file)
            for mons in mons_for_seasn:
                os.system("ls "+monthly_data_path+"*"+str(yr)+"-"+mons+".nc|cat >>"+list_file)

            with open(list_file) as f_obj:
                lines=f_obj.readlines()
            lists=listToString(lines)
           
          # compute means from the listed files above    
            climo_file=caseid+"_"+seasn+"_"+str(yr)+".nc"
            cmd="ncra "+lists+" "+out_path+climo_file
            os.system(cmd)
            os.system("mv "+list_file+" "+out_path)

    elif seasn == "SON":
        mons_for_seasn=["09","10","11"]
        for yr in years:
          # create list of input files
            list_file="list_"+seasn+"_"+str(yr)+".txt"
            if os.path.exists(list_file):
                os.system("rm "+list_file)
            for mons in mons_for_seasn:
                os.system("ls "+monthly_data_path+"*"+str(yr)+"-"+mons+".nc|cat >>"+list_file)

            with open(list_file) as f_obj:
                lines=f_obj.readlines()
            lists=listToString(lines)
           
          # compute means from the listed files above    
            climo_file=caseid+"_"+seasn+"_"+str(yr)+".nc"
            cmd="ncra "+lists+" "+out_path+climo_file
            os.system(cmd)
            os.system("mv "+list_file+" "+out_path)


#======================================================
#
#
#======================================================
print ("============== starting ==============")
import os
import numpy as np
from get_parameters import *
from d1_lon_lat_contour_model_vs_model import *
from d2_polar_contour_model_vs_model import *
from d3_northw_energy_transport_model_vs_model import *
from d4_time_series_areamean_model_vs_model import *
from d5_polar_stream_model_vs_model import *
from d6_polar_contour_lts_model_vs_model import *
from d7_polar_stream_contour_model_vs_model import *
from d8_polar_prof_contour_model_vs_model import *
from d9_polar_prof_lines_model_vs_model import *
from d10_polar_stream_contour_vert_model_vs_model import *


#+++++++++++++++++++++++
# common set-ups
#+++++++++++++++++++++++ 

#-- diag_package directories
os.environ["DIAG_HOME"]=os.getcwd()
#os.environ["SRCDIR"]=os.environ["DIAG_HOME"]+"/src" # not connected to main.py at pre    sent.
os.environ["WORKDIR"]=os.environ["DIAG_HOME"]+"/workdir"
os.environ["HTMLDIR"]=os.environ["DIAG_HOME"]+"/html_template"
os.environ["DATADIR"]=os.environ["DIAG_HOME"]+"/inputdata"

#-- case dependent setups
os.environ["exp_name"]="MC6_Scat"
os.environ["ctl_name"]="MC6_noScat"
#os.environ["fpath_ctl"]=os.environ["DATADIR"]+"/"+os.environ["ctl_name"]+"_test.nc"
#os.environ["fpath_exp"]=os.environ["DATADIR"]+"/"+os.environ["exp_name"]+"_test.nc"
#os.environ["fpath_ctl"]="/raid00/xianwen/Yi-Hsuan/E3SM_DECKv1b_H1.ne30/remap_180x360/climo"
os.environ["fpath_ctl"]="/raid00/xianwen/Yi-Hsuan/E3SM_coupled_restart_20TR_Yr2000-noScat.Year2000_2014/remap_180x360/climo"
os.environ["fpath_exp"]="/raid00/xianwen/Yi-Hsuan/E3SM_coupled_restart_20TR_Yr2000-Scat.Year2000_2014/remap_180x360/climo"
#os.environ["ctl_run_id"]="E3SM_DECKv1b_H1.ne30"
os.environ["ctl_run_id"]="E3SM_coupled_restart_20TR_Yr2000-noScat.Year2000_2014"
os.environ["exp_run_id"]="E3SM_coupled_restart_20TR_Yr2000-Scat.Year2000_2014"
os.environ["diagcase_name"]=os.environ["exp_name"]+"-"+os.environ["ctl_name"]
os.environ["OUTDIR"]=os.environ["WORKDIR"]+"/"+os.environ["diagcase_name"]

#varnms_2d=["FLUT","FLUTC","FSNTOA","FSNTOAC","SOLIN","LWCF","SWCF","FLNS","FLNSC","FLDS","FSDS","FSDSC","FSNS","FSNSC","SHFLX","LHFLX","TS"]
#varnms_3d=["QRL","QRS","CLOUD","CLDLIQ","CLDICE","T","P","RH","Q","U","V","W"]
#seasons=["ANN","DJF","JJA"]
varnms_2d=["TGCLDIWP"]
#varnms_2d=["TMQ"]
#varnms_2d=["TVH"]
#varnms_2d=["CLDMED"]
varnms_3d=["CLDLIQ"]
seasons=["DJF"]
scale_ctl=1000.
scale_exp=1000.

##-- one time test --
##varnms_2d=["FLUT"]
##seasons=["ANN"]

#-- select which diag to do (1:do, 0:not do)
do_lon_lat_contour=1
do_polar_contour_N=0
do_polar_contour_S=0
do_northw_energy_transport=0
do_time_series_areamean=0
do_polar_stream_N=0
do_polar_contour_lts_N=0
do_polar_stream_contour_N=0
do_polar_prof_contour_N=0
do_polar_prof_lines_N=0
do_polar_stream_contour_vert_N=0

#-- set format of figure files
os.environ["fig_show"]="True"
os.environ["fig_save"]="False"
os.environ["fig_suffix"]="png" # supported format: png, eps, pdf, etc.

#-- create work directory (for figures and html)
if not os.path.exists(os.environ["OUTDIR"]):
    os.makedirs(os.environ["OUTDIR"])
    os.makedirs(os.environ["OUTDIR"]+"/figures")
    os.makedirs(os.environ["OUTDIR"]+"/tables")
else:
    print("Use existing directory:")
    print(os.environ["OUTDIR"])

#pre-open a table file for the saving of means and differences achieved from each diag.
header=" "*18+" "*(10-len(os.environ["ctl_name"]))+os.environ["ctl_name"]+" "*5+\
       " "*(10-len(os.environ["exp_name"]))+os.environ["exp_name"]+" "*5+\
       " "*(10-len("Diff"))+"Diff"+" "*5+\
       " "*(10-len("Diff%"))+"Diff%"
if "ANN" in seasons:
    table_file_ANN=open(os.environ["OUTDIR"]+"/tables/"+os.environ["diagcase_name"]+"_time_means_sl_ANN.txt","w+")
    table_file_ANN.write(header+"\r\n")
if "DJF" in seasons:
    table_file_DJF=open(os.environ["OUTDIR"]+"/tables/"+os.environ["diagcase_name"]+"_time_means_sl_DJF.txt","w+")
    table_file_DJF.write(header+"\r\n")
if "JJA" in seasons:
    table_file_JJA=open(os.environ["OUTDIR"]+"/tables/"+os.environ["diagcase_name"]+"_time_means_sl_JJA.txt","w+")
    table_file_JJA.write(header+"\r\n")

if do_lon_lat_contour:
    print("--do lon lat contour--")
    if "ANN" in seasons: 
        table_file_ANN.write("Global"+"\r\n")
    if "DJF" in seasons: 
        table_file_DJF.write("Global"+"\r\n")
    if "JJA" in seasons: 
        table_file_JJA.write("Global"+"\r\n")
    for varnm in varnms_2d:
        print(varnm)
        for seasn in seasons:
            if seasn == "ANN":
                lon_lat_contour_model_vs_model(varnm,seasn,scale_ctl,scale_exp,table_file_ANN)
                #table_file_ANN.write(line+"\r\n")
            if seasn == "DJF":
                lon_lat_contour_model_vs_model(varnm,seasn,scale_ctl,scale_exp,table_file_DJF)
                #table_file_DJF.write(line+"\r\n")
            if seasn == "JJA":
                lon_lat_contour_model_vs_model(varnm,seasn,scale_ctl,scale_exp,table_file_JJA)
                #table_file_JJA.write(line+"\r\n")

if do_polar_contour_N:
    print("--do polar contour N--")
    if "ANN" in seasons: 
        table_file_ANN.write("North Pole"+"\r\n")
    if "DJF" in seasons: 
        table_file_DJF.write("North Pole"+"\r\n")
    if "JJA" in seasons: 
        table_file_JJA.write("North Pole"+"\r\n")
    for varnm in varnms_2d:
        print(varnm)
        for seasn in seasons:
            if seasn == "ANN":
                polar_contour_model_vs_model(varnm,seasn,scale_ctl,scale_exp,"N",table_file_ANN)
                #table_file_ANN.write(line+"\r\n")
            if seasn == "DJF":
                polar_contour_model_vs_model(varnm,seasn,scale_ctl,scale_exp,"N",table_file_DJF)
                #table_file_DJF.write(line+"\r\n")
            if seasn == "JJA":
                polar_contour_model_vs_model(varnm,seasn,scale_ctl,scale_exp,"N",table_file_JJA)
                #table_file_JJA.write(line+"\r\n")

if do_polar_contour_S:
    print("--do polar contour S--")
    if "ANN" in seasons: 
        table_file_ANN.write("South Pole"+"\r\n")
    if "DJF" in seasons: 
        table_file_DJF.write("South Pole"+"\r\n")
    if "JJA" in seasons: 
        table_file_JJA.write("South Pole"+"\r\n")
    for varnm in varnms_2d:
        print(varnm)
        for seasn in seasons:
            if seasn == "ANN":
                polar_contour_model_vs_model(varnm,seasn,scale_ctl,scale_exp,"S",table_file_ANN)
                #table_file_ANN.write(line+"\r\n")
            if seasn == "DJF":
                polar_contour_model_vs_model(varnm,seasn,scale_ctl,scale_exp,"S",table_file_DJF)
                #table_file_DJF.write(line+"\r\n")
            if seasn == "JJA":
                polar_contour_model_vs_model(varnm,seasn,scale_ctl,scale_exp,"S",table_file_JJA)
                #table_file_JJA.write(line+"\r\n")

if do_northw_energy_transport:
    print("--do northward energy transport--")
    for seasn in seasons:
        northw_energy_transport_model_vs_model(seasn)


if do_time_series_areamean:
    print("--do time series of areamean --")
    seasn="09"
    varnm="ICEFRAC"
    time_series_areamean_model_vs_model(varnm,seasn,scale_ctl,scale_exp)

if do_polar_stream_N:
    print("--do polar stream --")
    seasn="DJF"
    polar_stream_model_vs_model(seasn,scale_ctl,scale_exp,"N")

if do_polar_contour_lts_N:
    print("--do polar contour lts --")
    seasn="DJF"
    varnm="LTS"
    polar_contour_lts_model_vs_model(varnm,seasn,scale_ctl,scale_exp,"N")

if do_polar_stream_contour_N:
    print("--do polar stream contour --")
    seasn="DJF"
    polar_stream_contour_model_vs_model(seasn,scale_ctl,scale_exp,"N")

if do_polar_prof_contour_N:
    print("--do polar prof contour --")
    seasn="DJF"
    varnm="FICE"
    polar_prof_contour_model_vs_model(varnm,seasn,scale_ctl,scale_exp,"N")

if do_polar_prof_lines_N:
    print("--do polar lines --")
    seasn="DJF"
    varnm="FICE"
    polar_prof_lines_model_vs_model(varnm,seasn,scale_ctl,scale_exp,"N")

if do_polar_stream_contour_vert_N:
    print("--do stream contour polar lines --")
    seasn="DJF"
    #varnm="Q"
    polar_stream_contour_vert_model_vs_model(seasn,scale_ctl,scale_exp,"N")
#--------------------------------
# creat html file
#--------------------------------
# copy html template to work directory
os.system("cp "+os.environ["HTMLDIR"]+"/E3SM_diag.html "+os.environ["OUTDIR"]+"/")
os.system("cp "+os.environ["HTMLDIR"]+"/E3SM_logo.png "+os.environ["OUTDIR"]+"/")
os.system("echo '<H2><font color=navy>  "+os.environ["diagcase_name"]+"  <A></H3>' >> " \
        +os.environ["OUTDIR"]+"/E3SM_diag.html")

## write tables first
os.system("echo '<H3><font color=navy>------ time averaged tables ------ <A></H3>'     >> " \
        +os.environ["OUTDIR"]+"/E3SM_diag.html")
os.system("echo '<H3><font color=navy>&emsp;&emsp;&emsp;ANN  DJF  JJA  <A></H3>'     >> " \
        +os.environ["OUTDIR"]+"/E3SM_diag.html")
os.system("echo '<H3><font color=navy>&emsp;&emsp;&emsp; "+\
          " <A HREF=\"tables/"+os.environ["diagcase_name"]+"_time_means_sl_ANN.txt"+"\">table</A>"+\
          " <A HREF=\"tables/"+os.environ["diagcase_name"]+"_time_means_sl_DJF.txt"+"\">table</A>"+\
          " <A HREF=\"tables/"+os.environ["diagcase_name"]+"_time_means_sl_JJA.txt"+"\">table</A>"+\
          "</H3>' >> "+os.environ["OUTDIR"]+"/E3SM_diag.html")


## diagnosis 1
if do_lon_lat_contour:
    os.system("echo '<H3><font color=navy>------ Lon-Lat contour ------ <A></H3>'     >> " \
            +os.environ["OUTDIR"]+"/E3SM_diag.html")
    os.system("echo '<H3><font color=navy>&emsp;&emsp;&emsp;ANN  DJF  JJA  <A></H3>'     >> " \
            +os.environ["OUTDIR"]+"/E3SM_diag.html")
    #for i in range(nvrs_ml):
    for varnm in varnms_2d:
        figname_ANN="d1_lon_lat_contour_"+varnm+"_ANN."+os.environ["fig_suffix"]
        figname_DJF="d1_lon_lat_contour_"+varnm+"_DJF."+os.environ["fig_suffix"]
        figname_JJA="d1_lon_lat_contour_"+varnm+"_JJA."+os.environ["fig_suffix"]
        os.system("echo '<H3><font color=navy> "+varnm+\
                  " <A HREF=\"figures/"+figname_ANN+"\">plots</A>"+\
                  " <A HREF=\"figures/"+figname_DJF+"\">plots</A>"+\
                  " <A HREF=\"figures/"+figname_JJA+"\">plots</A>"+\
                  "</H3>' >> "+os.environ["OUTDIR"]+"/E3SM_diag.html")
     #figname="d1_lon_lat_contour_N_"+varnm+"."+os.environ["fig_suffix"]
     #os.system("echo '<H3><font color=navy> "+varnm+" <A HREF=\"figures/"+figname+\
     #     "\">plots</A></H3>' >> "+os.environ["OUTDIR"]+"/E3SM_diag.html")
## diagnosis 1
if do_polar_contour_N:
    os.system("echo '<H3><font color=navy>------ North Pole contour ------ <A></H3>'     >> " \
            +os.environ["OUTDIR"]+"/E3SM_diag.html")
    os.system("echo '<H3><font color=navy>&emsp;&emsp;&emsp;ANN  DJF  JJA  <A></H3>'     >> " \
            +os.environ["OUTDIR"]+"/E3SM_diag.html")
    #for i in range(nvrs_ml):
    #figname="d2_polar_contour_N_"+varnm+"."+os.environ["fig_suffix"]
    #os.system("echo '<H3><font color=navy> "+varnm+" <A HREF=\"figures/"+figname+\
    #     "\">plots</A></H3>' >> "+os.environ["OUTDIR"]+"/E3SM_diag.html")

    for varnm in varnms_2d:
        figname_ANN="d2_polar_contour_N_"+varnm+"_ANN."+os.environ["fig_suffix"]
        figname_DJF="d2_polar_contour_N_"+varnm+"_DJF."+os.environ["fig_suffix"]
        figname_JJA="d2_polar_contour_N_"+varnm+"_JJA."+os.environ["fig_suffix"]
        os.system("echo '<H3><font color=navy> "+varnm+\
                  " <A HREF=\"figures/"+figname_ANN+"\">plots</A>"+\
                  " <A HREF=\"figures/"+figname_DJF+"\">plots</A>"+\
                  " <A HREF=\"figures/"+figname_JJA+"\">plots</A>"+\
                  "</H3>' >> "+os.environ["OUTDIR"]+"/E3SM_diag.html")
## diagnosis 1
if do_polar_contour_S:
    os.system("echo '<H3><font color=navy>------ South Pole contour ------ <A></H3>'     >> " \
            +os.environ["OUTDIR"]+"/E3SM_diag.html")
    os.system("echo '<H3><font color=navy>&emsp;&emsp;&emsp;ANN  DJF  JJA  <A></H3>'     >> " \
            +os.environ["OUTDIR"]+"/E3SM_diag.html")
    for varnm in varnms_2d:
        figname_ANN="d2_polar_contour_S_"+varnm+"_ANN."+os.environ["fig_suffix"]
        figname_DJF="d2_polar_contour_S_"+varnm+"_DJF."+os.environ["fig_suffix"]
        figname_JJA="d2_polar_contour_S_"+varnm+"_JJA."+os.environ["fig_suffix"]
        os.system("echo '<H3><font color=navy> "+varnm+\
                  " <A HREF=\"figures/"+figname_ANN+"\">plots</A>"+\
                  " <A HREF=\"figures/"+figname_DJF+"\">plots</A>"+\
                  " <A HREF=\"figures/"+figname_JJA+"\">plots</A>"+\
                  "</H3>' >> "+os.environ["OUTDIR"]+"/E3SM_diag.html")

## diagnosis 1
if do_northw_energy_transport:
    os.system("echo '<H3><font color=navy>------ Northward Energy Transport ------ <A></H3>'     >> " \
            +os.environ["OUTDIR"]+"/E3SM_diag.html")
    os.system("echo '<H3><font color=navy>&emsp;&emsp;&emsp;ANN  DJF  JJA  <A></H3>'     >> " \
            +os.environ["OUTDIR"]+"/E3SM_diag.html")
    figname_ANN="d3_northw_energy_transport_"+varnm+"_ANN."+os.environ["fig_suffix"]
    figname_DJF="d3_northw_energy_transport_"+varnm+"_DJF."+os.environ["fig_suffix"]
    figname_JJA="d3_northw_energy_transport_"+varnm+"_JJA."+os.environ["fig_suffix"]
    os.system("echo '<H3><font color=navy> "+"NHT"+\
              " <A HREF=\"figures/"+figname_ANN+"\">plots</A>"+\
              " <A HREF=\"figures/"+figname_DJF+"\">plots</A>"+\
              " <A HREF=\"figures/"+figname_JJA+"\">plots</A>"+\
              "</H3>' >> "+os.environ["OUTDIR"]+"/E3SM_diag.html")
    
print ("============== Finished ==============")

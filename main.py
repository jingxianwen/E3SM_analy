#======================================================
#
#
#======================================================
print ("============== starting ==============")
import os
import numpy as np
from d1_lon_lat_model_vs_model import *
from d2_polar_contour_model_vs_model import *


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
os.environ["ctl_name"]="Standard"
os.environ["fpath_ctl"]=os.environ["DATADIR"]+"/"+os.environ["ctl_name"]+"_test.nc"
os.environ["fpath_exp"]=os.environ["DATADIR"]+"/"+os.environ["exp_name"]+"_test.nc"
os.environ["diagcase_name"]=os.environ["exp_name"]+"-"+os.environ["ctl_name"]
os.environ["OUTDIR"]=os.environ["WORKDIR"]+"/"+os.environ["diagcase_name"]

varnm="FLUT"
season="ANN"
scale_ctl=1.
scale_exp=1.
#-- select which diag to do
do_lon_lat=False
do_polar_N=True
do_polar_S=False

#-- set format of figure files
os.environ["fig_show"]="False"
os.environ["fig_save"]="True"
os.environ["fig_suffix"]="png" # supported format: png, eps, pdf, etc.

#-- create work directory (for figures and html)
if not os.path.exists(os.environ["OUTDIR"]):
    os.makedirs(os.environ["OUTDIR"])
    os.makedirs(os.environ["OUTDIR"]+"/figures")
else:
    print("Use existing directory:")
    print(os.environ["OUTDIR"])
#os.environ["EXP_DIR"]=os.environ["WORKDIR"]+"/"+os.environ["diagcase_name"]

if do_lon_lat:
   lon_lat_model_vs_model(varnm,season,scale_ctl,scale_exp)
if do_polar_N:
   polar_contour_model_vs_model(varnm,season,scale_ctl,scale_exp,"N")
if do_polar_S:
   polar_contour_model_vs_model(varnm,season,scale_ctl,scale_exp,"S")


#--------------------------------
# creat html file
#--------------------------------
# copy html template to work directory
os.system("cp "+os.environ["HTMLDIR"]+"/E3SM_diag.html "+os.environ["OUTDIR"]+"/")
os.system("cp "+os.environ["HTMLDIR"]+"/E3SM_logo.png "+os.environ["OUTDIR"]+"/")
os.system("echo '<H2><font color=navy>  "+os.environ["diagcase_name"]+"  <A></H3>' >> " \
        +os.environ["OUTDIR"]+"/E3SM_diag.html")
## diagnosis 1
if do_polar_N:
     os.system("echo '<H2><font color=navy>------ Polar contour ------ <A></H3>'     >> " \
             +os.environ["OUTDIR"]+"/E3SM_diag.html")
     #for i in range(nvrs_ml):
     figname="d2_polar_contour_N_"+varnm+"."+os.environ["fig_suffix"]
     os.system("echo '<H3><font color=navy> "+varnm+" <A HREF=\"figures/"+figname+\
          "\">plots</A></H3>' >> "+os.environ["OUTDIR"]+"/E3SM_diag.html")


print ("============== Finished ==============")

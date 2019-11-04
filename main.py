#======================================================
#
#
#======================================================
print ("============== starting ==============")
import os
import numpy as np
from d1_lon_lat_model_vs_model import *


#+++++++++++++++++++++++
# common set-ups
#+++++++++++++++++++++++ 

#-- diag_package directories
os.environ["DIAG_HOME"]=os.getcwd()
#os.environ["SRCDIR"]=os.environ["DIAG_HOME"]+"/src" # not connected to main.py at pre    sent.
os.environ["WORKDIR"]=os.environ["DIAG_HOME"]+"/workdir"
os.environ["HTMLDIR"]=os.environ["DIAG_HOME"]+"/html_template"
os.environ["DATADIR"]=os.environ["DIAG_HOME"]+"/inputdata"

os.environ["exp_name"]="MC6_Scat"
os.environ["ctl_name"]="Standard"
os.environ["fpath_ctl"]=os.environ["DATADIR"]+"/"+os.environ["ctl_name"]+"_test.nc"
os.environ["fpath_exp"]=os.environ["DATADIR"]+"/"+os.environ["exp_name"]+"_test.nc"

varnm="FLUT"
season="ANN"
scale_ctl=1.
scale_exp=1.
lon_lat_model_vs_model(varnm,season,scale_ctl,scale_exp)

print ("============== Finished ==============")

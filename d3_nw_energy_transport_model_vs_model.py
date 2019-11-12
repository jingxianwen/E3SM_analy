#=====================================================
#
#=====================================================
# os
import os
#import netCDF4
from netCDF4 import Dataset as netcdf_dataset
# cartopy
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point
# matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors
# numpy
import numpy as np
# parameters
from get_parameters import *

def northw_energy_transport_model_vs_model(varnm,season,scale_ctl,scale_exp,table):
    # data path
    ctl_name=os.environ["ctl_name"]
    exp_name=os.environ["exp_name"]
    fpath_ctl=os.environ["fpath_ctl"]+"/"+os.environ["ctl_run_id"]+"_climo_"+season+".nc"
    fpath_exp=os.environ["fpath_exp"]+"/"+os.environ["exp_run_id"]+"_climo_"+season+".nc"
    
    # open data file
    file_ctl=netcdf_dataset(fpath_ctl,"r")
    file_exp=netcdf_dataset(fpath_exp,"r")

    # read lat and lon
    lat=file_ctl.variables["lat"]
    lon=file_ctl.variables["lon"]

    # read required variables
    # variable list: FLUT, FSNTOA, FSNS, FLNS, LHFLX, SHFLX
    #-- TOA
    flut_ctl=flie_ctl.variables["FLUT"]
    fsnt_ctl=flie_ctl.variables["FSNTOA"]
    flut_exp=flie_exp.variables["FLUT"]
    fsnt_exp=flie_exp.variables["FSNTOA"]
    #-- Surface
    fsns_ctl=file_ctl.variables["FSNS"]
    flns_ctl=file_ctl.variables["FLNS"]
    lhflx_ctl=file_ctl.variables["lhflx"]
    shflx_ctl=file_ctl.variables["shflx"]
    fsns_exp=file_exp.variables["FSNS"]
    flns_exp=file_exp.variables["FLNS"]
    lhflx_exp=file_exp.variables["lhflx"]
    shflx_exp=file_exp.variables["shflx"]

    # compute zonal mean and change sign 
    # TOA:positive downward
    # Surface: positive upward
    flut_ctl_zm=np.mean(flut_ctl,axis=2)
    fsnt_ctl_zm=np.mean(fsnt_ctl,axis=2)
    fsns_ctl_zm=np.mean(fsns_ctl,axis=2)
    flns_ctl_zm=np.mean(flns_ctl,axis=2)
    lhflx_ctl_zm=np.mean(lhflx_ctl,axis=2)
    shflx_ctl_zm=np.mean(shflx_ctl,axis=2)
    flut_exp_zm=np.mean(flut_exp,axis=2)
    fsnt_exp_zm=np.mean(fsnt_exp,axis=2)
    fsns_exp_zm=np.mean(fsns_exp,axis=2)
    flns_exp_zm=np.mean(flns_exp,axis=2)
    lhflx_exp_zm=np.mean(lhflx_exp,axis=2)
    shflx_exp_zm=np.mean(shflx_exp,axis=2)
    exit()
#===============================

#def get_parameters(varnm,season):
#    #list_rad=["FLUT","FLUTC","FLNT","FLNTC","FSNT","FSNTC","FSDS","FSDSC","FSNS","FSNSC"]
#    if varnm == "FLUT":
#        parameters={"units":"W/m2",\
#		   "contour_levs":[120, 140, 160, 180, 200, 220, 240, 260, 280, 300],\
#		   "diff_levs":[-50, -40, -30, -20, -10, -5, 5, 10, 20, 30, 40, 50],\
#                   "colormap":"PiYG_r",\
#                   "colormap_diff":"bwr"\
#		   }
#    return parameters
#
#def get_area_mean_range(varnm,lat):
#   # 1. area weighted average 
#    #convert latitude to radians
#    latr=np.deg2rad(lat)
#    #use cosine of latitudes as weights for the mean
#    weights=np.cos(latr)
#    #first calculate zonal mean
#    zonal_mean=varnm.mean(axis=2)
#    #then calculate weighted global mean
#    area_mean=np.average(zonal_mean,axis=1,weights=weights)
#   # 2. min and max
#    minval=varnm.min()
#    maxval=varnm.max()
#    return area_mean,minval,maxval

#=====================================================
#
#=====================================================
# os
import os
import netCDF4 as nc
from netCDF4 import Dataset as netcdf_dataset
from netCDF4 import MFDataset
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

def time_series_areamean_model_vs_model(varnm,season,scale_ctl,scale_exp):
    # data path
    ctl_name=os.environ["ctl_name"]
    exp_name=os.environ["exp_name"]
    fpath_ctl=os.environ["fpath_ctl"] #+"/"+os.environ["ctl_run_id"]+"_climo_"+season+".nc"
    fpath_exp=os.environ["fpath_exp"] #+"/"+os.environ["exp_run_id"]+"_climo_"+season+".nc"
     
    # open data file
    #file_ctl=MFDataset(fpath_ctl+"/../*"+season+".nc","r")
    #file_exp=MFDataset(fpath_exp+"/../*"+season+".nc","r")
    file_ctl=MFDataset(fpath_ctl+"/*2014_20*"+season+".nc","r")
    file_exp=MFDataset(fpath_exp+"/*2014_20*"+season+".nc","r")

    # read lat and lon
    lat=file_ctl.variables["lat"]
    lon=file_ctl.variables["lon"]

    # read data and calculate mean/min/max
    dtctl=file_ctl.variables[varnm] #*scale_ctl
    dtexp=file_exp.variables[varnm] #*scale_exp
    dtdif=dtexp[:,:,:]-dtctl[:,:,:]

    # calculate mean for globe
    #latbound=np.min(np.where(lat[:]>60))+1
    stats_ctl=get_area_mean_min_max(dtctl[:,:,:],lat[:])
    stats_exp=get_area_mean_min_max(dtexp[:,:,:],lat[:])
    stats_dif=get_area_mean_min_max(dtdif[:,:,:],lat[:])
    stats_difp=stats_dif[0]/stats_ctl[0]*100.

    # calculate mean for North of 60N
    #latbound=np.min(np.where(lat[:]>60))+1
    #stats_ctl=get_area_mean_min_max(dtctl[:,latbound:,:],lat[latbound:])
    #stats_exp=get_area_mean_min_max(dtexp[:,latbound:,:],lat[latbound:])
    #stats_dif=get_area_mean_min_max(dtdif[:,latbound:,:],lat[latbound:])
    #stats_difp=stats_dif[0]/stats_ctl[0]*100.

    # calculate mean for sourth of 60S
    #latbound=np.max(np.where(lat[:]<-60))+1
    #stats_ctl=get_area_mean_min_max(dtctl[:,0:latbound,:],lat[0:latbound])
    #stats_exp=get_area_mean_min_max(dtexp[:,0:latbound,:],lat[0:latbound])
    #stats_dif=get_area_mean_min_max(dtdif[:,0:latbound,:],lat[0:latbound])
    #stats_difp=stats_dif[0]/stats_ctl[0]*100.

    years=np.arange(2000,2015)
    fig = plt.figure(figsize=[7.0,11.0],dpi=150.)
    ax=fig.add_subplot(111)
    print(stats_ctl[0])
    print(years)
    line_ctl,=ax.plot(years,stats_ctl[0],ls="-",c="k",marker="D",label=ctl_name)
    line_exp,=ax.plot(years,stats_exp[0],ls="-",c="r",marker="D",label=exp_name)
    ax.set_xlabel("Year")
    ax.set_ylabel(varnm+" (fraction)")
 
    ax.legend(loc="best",fontsize=10)

##    if os.environ["fig_save"]=="True":
##        fname="d1_lon_lat_contour_"+varnm+"_"+season+"."+os.environ["fig_suffix"]
##        plt.savefig(os.environ["OUTDIR"]+"/figures/"+fname)
    if os.environ["fig_show"]=="True":
        plt.show()
    plt.close()


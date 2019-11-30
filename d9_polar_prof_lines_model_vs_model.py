#====================================================
#
#====================================================
# os
import os
#import netCDF4
from netCDF4 import Dataset as netcdf_dataset
# cartopy
import os
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point
# matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.path as mpath
import matplotlib.colors as colors
# numpy
import numpy as np
# parameters
from get_parameters import *

def polar_prof_lines_model_vs_model(varnm,season,scale_ctl,scale_exp,pole):
    # data path
    ctl_name=os.environ["ctl_name"]
    exp_name=os.environ["exp_name"]
    #fpath_ctl=os.environ["fpath_ctl"]
    #fpath_exp=os.environ["fpath_exp"]
    fpath_ctl=os.environ["fpath_ctl"]+"/"+os.environ["ctl_run_id"]+"_climo_"+season+".nc"
    fpath_exp=os.environ["fpath_exp"]+"/"+os.environ["exp_run_id"]+"_climo_"+season+".nc"

    # open data file
    file_ctl=netcdf_dataset(fpath_ctl,"r")
    file_exp=netcdf_dataset(fpath_exp,"r")

    # read lat and lon
    lat=file_ctl.variables["lat"]
    lon=file_ctl.variables["lon"]
    lev=file_ctl.variables["lev"]
    nlat=lat[:].size
    nlon=lon[:].size
    nlev=lev[:].size

    # read data and calculate difference
    dtctl=file_ctl.variables[varnm][:,:,:,:]*scale_ctl
    dtexp=file_exp.variables[varnm][:,:,:,:]*scale_exp
    dtdif=dtexp[:,:,:,:]-dtctl[:,:,:,:]

    if pole == "N":
        latbound1=np.min(np.where(lat[:]>60))
        latbound2=nlat-1
    elif pole == "S":
        latbound1=0
        latbound2=np.max(np.where(lat[:]<-60))+1
    
    # compute area mean profiles
    dtctl_am=get_area_mean_prof(dtctl[:,:,latbound1:latbound2,:],lat[latbound1:latbound2])
    dtexp_am=get_area_mean_prof(dtexp[:,:,latbound1:latbound2,:],lat[latbound1:latbound2])
    dtdif_am=get_area_mean_prof(dtdif[:,:,latbound1:latbound2,:],lat[latbound1:latbound2])

    fig = plt.figure(figsize=[5.0,6.0],dpi=150.)
    
    #ax = fig.add_subplot(111)
    rect=(0.25,0.15,0.6,0.7)
    ax = fig.add_axes(rect)
    ax.set_aspect("auto")
    line_ctl= ax.plot(dtctl_am[0,:],lev,ls="-",lw=2,c="k",label="noScat")
    line_exp= ax.plot(dtexp_am[0,:],lev,ls="--",lw=2,c="r",label="Scat")
    line_dif= ax.plot(dtdif_am[0,:],lev,ls="-",lw=2,c="c",label="Scat-noScat")
    line_zero= ax.plot(dtdif_am[0,:]*0.0,lev,ls="-",lw=1,c="gray")

    ax.set_ylim(1000,1)
    ax.set_xlabel(varnm,fontsize=12,fontweight="bold",fontstyle="normal")
    ax.set_ylabel("Pressure (hPa)",fontsize=12,fontweight="bold",fontstyle="normal")
    xtlocs=np.float32(ax.xaxis.get_ticklocs())
    ytlocs=np.int64(ax.yaxis.get_ticklocs())
    ax.set_xticklabels(xtlocs,fontsize=12,fontweight="bold")
    ax.set_yticklabels(ytlocs,fontsize=12,fontweight="bold")

    ax.legend(loc="best",fontsize=12)

    #save figure as file
    if os.environ["fig_save"]=="True":
        fname="d9_polar_prof_lines_"+pole+"_"+varnm+"_"+season+"."+os.environ["fig_suffix"]
        plt.savefig(os.environ["OUTDIR"]+"/figures/"+fname)
    if os.environ["fig_show"]=="True":
        plt.show()
    plt.close()

    return()


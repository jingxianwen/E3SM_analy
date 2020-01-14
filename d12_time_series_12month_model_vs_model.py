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
from pylab import *
# numpy
import numpy as np
# parameters
from get_parameters import *

def time_series_12month_model_vs_model():
    # data path
    ctl_name=os.environ["ctl_name"]
    exp_name=os.environ["exp_name"]
    months=["01","02","03","04","05","06","07","08","09","10","11","12"]
    MSE_ctl=np.zeros((12))  # energy transport
    MSE_exp=np.zeros((12))
    MSE_dif=np.zeros((12))

    FSRF_ctl=np.zeros((12)) # net energy interaction at surface
    FSRF_exp=np.zeros((12))
    FSRF_dif=np.zeros((12))

    FTOA_ctl=np.zeros((12)) # net radiation at TOA
    FTOA_exp=np.zeros((12))
    FTOA_dif=np.zeros((12))

    ATM_ctl=np.zeros((12))  # changes in atmospheric energy storage 
    ATM_exp=np.zeros((12))
    ATM_dif=np.zeros((12))

    for im in range(0,12):
        fpath_ctl=os.environ["fpath_ctl"]+"/"+os.environ["ctl_run_id"]+"_climo_"+months[im]+".nc"
        fpath_exp=os.environ["fpath_exp"]+"/"+os.environ["exp_run_id"]+"_climo_"+months[im]+".nc"
    
        # open data file
        file_ctl=netcdf_dataset(fpath_ctl,"r")
        file_exp=netcdf_dataset(fpath_exp,"r")
    
        # read lat and lon
        lat=file_ctl.variables["lat"]
        lon=file_ctl.variables["lon"]
        nlat=lat[:].size
        nlon=lon[:].size
    
     #Arctic area mean 
        latlim=70.  #region limitation for latitude
        latbound1=np.min(np.where(lat[:]>latlim))
        latbound2=nlat
     #Antarctic area mean 
        #latbound1=0
        #latbound2=np.max(np.where(lat[:]<-60))+1

     # read data and calculate mean/min/max
     # for TOA energy budget-->
        dtctl=file_ctl.variables["FSNTOA"][:,:,:]-file_ctl.variables["FLUT"][:,:,:]
        dtexp=file_exp.variables["FSNTOA"][:,:,:]-file_exp.variables["FLUT"][:,:,:]
        dtdif=dtexp[:,:,:]-dtctl[:,:,:]
        stats_ctl=get_area_mean_min_max(dtctl[:,latbound1:latbound2,:],lat[latbound1:latbound2])
        stats_exp=get_area_mean_min_max(dtexp[:,latbound1:latbound2,:],lat[latbound1:latbound2])
        stats_dif=get_area_mean_min_max(dtdif[:,latbound1:latbound2,:],lat[latbound1:latbound2])
        FTOA_ctl[im]=stats_ctl[0]
        FTOA_exp[im]=stats_exp[0]
        FTOA_dif[im]=stats_dif[0]
    
     # for Surface energy budget-->
        dtctl=file_ctl.variables["FLNS"][:,:,:]-file_ctl.variables["FSNS"][:,:,:]\
             +file_ctl.variables["LHFLX"][:,:,:]+file_ctl.variables["SHFLX"][:,:,:]
        dtexp=file_exp.variables["FLNS"][:,:,:]-file_exp.variables["FSNS"][:,:,:]\
             +file_exp.variables["LHFLX"][:,:,:]+file_exp.variables["SHFLX"][:,:,:]
        dtdif=dtexp[:,:,:]-dtctl[:,:,:]
        stats_ctl=get_area_mean_min_max(dtctl[:,latbound1:latbound2,:],lat[latbound1:latbound2])
        stats_exp=get_area_mean_min_max(dtexp[:,latbound1:latbound2,:],lat[latbound1:latbound2])
        stats_dif=get_area_mean_min_max(dtdif[:,latbound1:latbound2,:],lat[latbound1:latbound2])
        stats_difp=np.array(stats_dif[0]/stats_ctl[0]*100.)
        FSRF_ctl[im]=stats_ctl[0]
        FSRF_exp[im]=stats_exp[0]
        FSRF_dif[im]=stats_dif[0]

     # for meridional energy transport-->
        rearth=6.37122e6  # radius of the Earth [m]
        dtctl=file_ctl.variables["TVH"][:,latbound1,:]
        dtexp=file_exp.variables["TVH"][:,latbound1,:]
        lat_rad=lat[latbound1]/180.*np.pi  #latitude from degree to radians
        #length of the latitude * mean TVH
        tvh_tot_ctl=2.*np.pi*(rearth*np.cos(lat_rad))*\
            np.mean(dtctl,axis=1)
        tvh_tot_exp=2.*np.pi*(rearth*np.cos(lat_rad))*\
            np.mean(dtexp,axis=1)
        #north cap area
        area_N=2.*np.pi*(rearth**2)*(1-np.sin(lat_rad))
        #area mean tvh
        MSE_ctl[im]=tvh_tot_ctl/area_N
        MSE_exp[im]=tvh_tot_exp/area_N
        MSE_dif[im]=MSE_exp[im]-MSE_ctl[im]
   
     # for atmospheric energy storage-->
        ATM_ctl[im]=FSRF_ctl[im]+FTOA_ctl[im]+MSE_ctl[im]
        ATM_exp[im]=FSRF_exp[im]+FTOA_exp[im]+MSE_exp[im]
        ATM_dif[im]=ATM_exp[im]-ATM_ctl[im]
    #print(FSRF_ctl)
    #print(FSRF_exp)
    #print(FTOA_ctl)
    #print(FTOA_exp)
    #print(MSE_ctl)
    #print(MSE_exp)
    #print(ATM_ctl)
    #print(ATM_exp)
  # make plot
    xlocs=range(1,13)
    yzero=np.zeros(12)
    xlabels=["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"]
    fig=plt.figure(figsize=(7,9))
    ax1=axes([0.1,0.55,0.8,0.4])
    ax2=axes([0.1,0.05,0.8,0.4])
    line_FSRF_ctl=ax1.plot(xlocs,FSRF_ctl,ls="-",lw=2,label="FSRF_ctl",c="y")
    line_FSRF_exp=ax1.plot(xlocs,FSRF_exp,ls=":",lw=2,label="FSRF_exp",c="y")
    line_FTOA_ctl=ax1.plot(xlocs,FTOA_ctl,ls="-",lw=2,label="FTOA_ctl",c="c")
    line_FTOA_exp=ax1.plot(xlocs,FTOA_exp,ls=":",lw=2,label="FTOA_exp",c="c")
    line_MSE_ctl=ax1.plot(xlocs,MSE_ctl,ls="-",lw=2,label="MSE_ctl",c="m")
    line_MSE_exp=ax1.plot(xlocs,MSE_exp,ls=":",lw=2,label="MSE_exp",c="m")
    line_ATM_ctl=ax1.plot(xlocs,ATM_ctl,ls="-",lw=2,label="ATM_ctl",c="darkorange")
    line_ATM_exp=ax1.plot(xlocs,ATM_exp,ls=":",lw=2,label="ATM_exp",c="darkorange")
    line_yzero=ax1.plot(xlocs,yzero,ls="-",lw=1,c="gray")
    ax1.set(xlim=[1,12])
    ax1.set_xticks(xlocs)
    ax1.set_xticklabels(xlabels)
    ax1.set_ylabel("W/m2")
    ax1.legend(loc="best")
    ax1.set_title("Atmospheric energy budgets in Arctic (>"+str(int(latlim))+"N)",fontsize=12,c="k")

    line_FSRF_dif=ax2.plot(xlocs,FSRF_dif,ls="-",lw=2,label="FSRF_dif",c="y")
    line_FTOA_dif=ax2.plot(xlocs,FTOA_dif,ls="-",lw=2,label="FTOA_dif",c="c")
    line_MSE_dif=ax2.plot(xlocs,MSE_dif,ls="-",lw=2,label="MSE_dif",c="m")
    line_ATM_dif=ax2.plot(xlocs,ATM_dif,ls="-",lw=2,label="ATM_dif",c="darkorange")
    line_yzero_dif=ax2.plot(xlocs,yzero,ls="-",lw=1,c="gray")
    ax2.set(xlim=[1,12])
    ax2.set_xticks(xlocs)
    ax2.set_xticklabels(xlabels)
    ax2.set_ylabel("W/m2")
    ax2.legend(loc="best")
    ax2.set_title("Diff in atmospheric energy budgets in Arctic (>"+str(int(latlim))+"N)",fontsize=12,c="k")
    
    fig.savefig("./workdir/others/time_series_12month_energybudget_arctic_"+str(int(latlim))+"N.png",dpi=150)
    show()

    #print(MSE_ctl)
    #print(MSE_exp)
    #print(MSE_dif)
    return()


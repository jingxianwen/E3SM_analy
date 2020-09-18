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

# scipy
from scipy import stats

# parameters
from get_parameters import get_area_mean_min_max

#def lon_lat_contour_model_vs_model(varnm,season,scale_ctl,scale_exp,table):
# data path
ctl_name="Abs" #os.environ["ctl_name"]
exp_name="Scat" #os.environ["exp_name"]
ctl_pref="CMIP_RRTMG_UMRad_scat_offline.ne30_ne30.cori-knl-ens0"
exp_pref="CMIP_RRTMG_UMRad_scat_offline.ne30_ne30.cori-knl-ens0"

fpath_ctl="/global/cscratch1/sd/xianwen/E3SM_simulations/"+ctl_pref+"/archive/remap_180x360/"
fpath_exp="/global/cscratch1/sd/xianwen/E3SM_simulations/"+exp_pref+"/archive/remap_180x360/"
 
years=np.arange(2001,2002) 
#years=np.arange(2010,2021) 
months_all=["01","02","03","04","05","06","07","08","09","10","11","12"]

# variable group 1:
varnm_ctl=np.array(["FLUT_OFF","FLDS_OFF"])
varnm_exp=np.array(["FLUT","FLDS"])

#varnms=np.array(["FSNTOA","FSNS","TS"])
var_long_name="Net Flux at TOA"
figure_name="net_flux_toa"
#units="K"
#var_long_name="Surface Net SW"
#figure_name="Surface_Net_SW_zonal_ANN"
#var_long_name="TOA Net SW"
#figure_name="TOA_Net_SW_zonal_ANN_VIS_icealb"
units=r"W/m$^2$"

nlat=np.int64(180)
#means_yby_ctl=np.zeros((years.size,varnms.size,nlat)) #year by year mean for each variable
#means_yby_exp=np.zeros((years.size,varnms.size,nlat)) #year by year mean for each variable
means_ctl=np.zeros((nlat)) #multi-year mean for each variable
means_exp=np.zeros((nlat)) #multi-year mean for each variable
diffs=np.zeros((2,nlat)) #multi-year exp-ctl diff for each variable

# open data file
fctl=fpath_ctl+"CMIP_RRTMG_UMRad_scat_offline.ne30_ne30.cori-knl-ens0.cam.h0.2000-07.nc"
fexp=fpath_exp+"CMIP_RRTMG_UMRad_scat_offline.ne30_ne30.cori-knl-ens0.cam.h0.2000-07.nc"
file_ctl=netcdf_dataset(fctl,"r")
file_exp=netcdf_dataset(fexp,"r")

# read lat and lon
lat=file_ctl.variables["lat"]
lon=file_ctl.variables["lon"]

#stats_ctl=np.zeros((14))
#stats_exp=np.zeros((14))
#stats_dif=np.zeros((14))
#stats_difp=np.zeros((14))
#print(stats_ctl)
# read data and calculate mean/min/max

# first variable
dtctl=file_ctl.variables[varnm_ctl[0]]
dtexp=file_exp.variables[varnm_exp[0]] 

means_ctl[:]=np.mean(dtctl[:,:,:],axis=2)[0,:]*-1.
means_exp[:]=np.mean(dtexp[:,:,:],axis=2)[0,:]*-1.
diffs[0,:]=means_exp-means_ctl

# second variable
dtctl=file_ctl.variables[varnm_ctl[1]]
dtexp=file_exp.variables[varnm_exp[1]] 

means_ctl[:]=np.mean(dtctl[:,:,:],axis=2)[0,:]
means_exp[:]=np.mean(dtexp[:,:,:],axis=2)[0,:]
diffs[1,:]=means_exp-means_ctl

zeros=np.zeros(diffs.shape)

# make the plot
fig=plt.figure(figsize=(7,8))
ax1=fig.add_axes([0.12,0.55,0.8,0.36])

ax1.plot(lat[:],diffs[0,:],color="k",lw=2,ls="-")
#ax1.plot(lat[:],means_exp[0,:],color="k",lw=2,ls=":",label="TSIS")
#ax1.legend(loc="upper left",fontsize=12)
#ax1.legend(fontsize=12)
ax1.set_title("Diff in TOA Net Flux (Scat-noScat, July)",fontsize=14)
ax1.set_ylabel(units,fontsize=12)
ax1.set_xlabel("Latitude",fontsize=12)
ax1.set_xlim(-90,90)
ax1.set_ylim(0.,3.1)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

ax2=fig.add_axes([0.12,0.09,0.8,0.36])

ax2.plot(lat[:],diffs[1,:],color="k",lw=2)
ax2.plot(lat[:],zeros[0,:],color="lightgray",lw=1)
ax2.set_title("Diff in SFC Net Flux (Scat-noScat, July)",fontsize=14) #+var_long_name,fontsize=12)
ax2.set_ylabel(units,fontsize=12)
ax2.set_xlabel("Latitude",fontsize=12)
ax2.set_xlim(-90,90)
ax2.set_ylim(-0.1,0.5)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

#plt.savefig(figure_name+".eps")
#plt.savefig(figure_name+"3.png",dpi=(200))
plt.show()

exit()

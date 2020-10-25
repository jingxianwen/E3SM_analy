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
import matplotlib.collections as collections

# numpy
import numpy as np

# scipy
from scipy import stats

# parameters
from get_parameters import get_area_mean_min_max

#def lon_lat_contour_model_vs_model(varnm,season,scale_ctl,scale_exp,table):
# data path
ctl_name="noScat" #os.environ["ctl_name"]
exp_name="Scat" #os.environ["exp_name"]
ctl_pref="CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl"
exp_pref="CMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl"

fpath_ctl="/global/cscratch1/sd/xianwen/E3SM_simulations/"+\
        ctl_pref+"/archive/climo_2010-2019/"
fpath_exp="/global/cscratch1/sd/xianwen/E3SM_simulations/"+\
        exp_pref+"/archive/climo_2010-2019/"

#years=np.arange(2010,2020) 
#months_all=["01","02","03","04","05","06","07","08","09","10","11","12"]

figure_name="zonal_tpw"
units=r"kg m$^-$$^2$"

# open data file
fctl=fpath_ctl+ctl_pref+"_climo_ANN.nc"
fexp=fpath_exp+exp_pref+"_climo_ANN.nc"
file_ctl=netcdf_dataset(fctl,"r")
file_exp=netcdf_dataset(fexp,"r")

# read lat and lon
lat=file_ctl.variables["lat"]
lon=file_ctl.variables["lon"]

# read data and calculate mean/min/max
vn="TMQ"
dtctl_wv=file_ctl.variables[vn][0,:,:]
dtexp_wv=file_exp.variables[vn][0,:,:] 
means_ctl_wv =  np.mean(dtctl_wv[:,:],axis=1)
means_exp_wv =  np.mean(dtexp_wv[:,:],axis=1)
diff_wv = means_exp_wv - means_ctl_wv

vn="TS"
dtctl_ts=file_ctl.variables[vn][0,:,:]
dtexp_ts=file_exp.variables[vn][0,:,:] 
means_ctl_ts = np.mean(dtctl_ts[:,:],axis=1)
means_exp_ts = np.mean(dtexp_ts[:,:],axis=1)
diff_ts = means_exp_ts - means_ctl_ts

#print(np.mean(gm_yby_ctl_wv,axis=1))
#print(np.mean(gm_yby_ctl_wv,axis=1))
#print(np.mean(gm_yby_ctl_net,axis=1))
#print(np.mean(gm_yby_exp_wv,axis=1))
#print(np.mean(gm_yby_exp_wv,axis=1))
#print(np.mean(gm_yby_exp_net,axis=1))
#exit()

# compute multi-year mean and ttest

zeros=np.zeros(diff_wv.shape)

diff_wv_cc=diff_ts[:]*means_ctl_wv[:]*0.075

# make the plot
fig=plt.figure(figsize=(8,7))

ax2=fig.add_axes([0.15,0.15,0.7,0.6])
ax2.plot(lat[:],diff_wv[:],color="k",lw=4,alpha=1.0,label="model")
ax2.plot(lat[:],diff_wv_cc[:],color="k",lw=2,ls=":",label="C-C estimate") #,label="\u0394TPW"
ax2.plot(lat[:],zeros[:],color="gray",lw=1)
ax2.legend(fontsize=14)
ax2.set_title("Differences (Scat - noScat)",fontsize=14) #+var_long_name,fontsize=12)
ax2.set_xlabel("Latitude",fontsize=14)
ax2.set_ylabel("\u0394TPW ("+units+")",fontsize=14)
ax2.set_xlim(-90,90)
#ax2.set_ylim(-0.6,0.6)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# add shading 

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.savefig(figure_name+".png")
plt.show()

exit()

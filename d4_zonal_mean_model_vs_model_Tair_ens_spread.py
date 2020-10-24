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
ctl_name="noScat" #os.environ["ctl_name"]
exp_name="Scat" #os.environ["exp_name"]
nlat=180

#---------------------------
#      CMIP 4 ensembles
#---------------------------
ctl_pref_CMIP=["CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl",\
               #"CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl-ens1",\
               "CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl-ens2"] #,\
               #"CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl-ens3"]
exp_pref_CMIP=["CMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl",\
               "CMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl-ens1"] #,\
               #"CMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl-ens2",\
               #"CMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl-ens3"]
num_ens=2
zm_tair_4ens_ctl_DJF=np.zeros((num_ens,nlat))
zm_tair_4ens_exp_DJF=np.zeros((num_ens,nlat))
zm_tair_4ens_ctl_JJA=np.zeros((num_ens,nlat))
zm_tair_4ens_exp_JJA=np.zeros((num_ens,nlat))

for iens in range(0,num_ens):
    fpath_ctl="/global/cscratch1/sd/xianwen/E3SM_simulations/"+\
            ctl_pref_CMIP[iens]+"/archive/climo_2010-2019/"
    fpath_exp="/global/cscratch1/sd/xianwen/E3SM_simulations/"+\
            exp_pref_CMIP[iens]+"/archive/climo_2010-2019/"
    for iss in ["DJF","JJA"]:
        fctl=fpath_ctl+ctl_pref_CMIP[iens]+"_climo_"+iss+".nc"
        fexp=fpath_exp+exp_pref_CMIP[iens]+"_climo_"+iss+".nc"
        file_ctl=netcdf_dataset(fctl,"r")
        file_exp=netcdf_dataset(fexp,"r")
        dtctl=file_ctl.variables["TREFHT"]
        dtexp=file_exp.variables["TREFHT"] 
        if iss is "DJF":
           zm_tair_4ens_ctl_DJF[iens,:]=np.mean(dtctl[:,:,:],axis=2)[0,:]
           zm_tair_4ens_exp_DJF[iens,:]=np.mean(dtexp[:,:,:],axis=2)[0,:]
        elif iss is "JJA":
           zm_tair_4ens_ctl_JJA[iens,:]=np.mean(dtctl[:,:,:],axis=2)[0,:]
           zm_tair_4ens_exp_JJA[iens,:]=np.mean(dtexp[:,:,:],axis=2)[0,:]

zm_tair_4ens_diff_DJF=zm_tair_4ens_exp_DJF - zm_tair_4ens_ctl_DJF
zm_tair_4ens_diff_JJA=zm_tair_4ens_exp_JJA - zm_tair_4ens_ctl_JJA

zm_ensmean_tair_diff_DJF=np.mean(zm_tair_4ens_diff_DJF,axis=0)
zm_ensmin_tair_diff_DJF=np.min(zm_tair_4ens_diff_DJF,axis=0)
zm_ensmax_tair_diff_DJF=np.max(zm_tair_4ens_diff_DJF,axis=0)

zm_ensmean_tair_diff_JJA=np.mean(zm_tair_4ens_diff_JJA,axis=0)
zm_ensmin_tair_diff_JJA=np.min(zm_tair_4ens_diff_JJA,axis=0)
zm_ensmax_tair_diff_JJA=np.max(zm_tair_4ens_diff_JJA,axis=0)

#---------------------------
#      AMIP 1 ensemble
#---------------------------
ctl_pref_AMIP=["AMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl"]
exp_pref_AMIP=["AMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl"]

for iss in ["DJF","JJA"]:
    fctl="/global/cscratch1/sd/xianwen/E3SM_simulations/"+\
            ctl_pref_AMIP[0]+"/archive/climo/"+ctl_pref_AMIP[0]+"_climo_"+iss+".nc"
    fexp="/global/cscratch1/sd/xianwen/E3SM_simulations/"+\
            exp_pref_AMIP[0]+"/archive/climo/"+exp_pref_AMIP[0]+"_climo_"+iss+".nc"
    file_ctl=netcdf_dataset(fctl,"r")
    file_exp=netcdf_dataset(fexp,"r")
    dtctl=file_ctl.variables["TREFHT"]
    dtexp=file_exp.variables["TREFHT"] 

    lat=file_ctl.variables["lat"]
    lon=file_ctl.variables["lon"]
    if iss is "DJF":
       zm_tair_AMIP_ctl_DJF=np.mean(dtctl[:,:,:],axis=2)[0,:]
       zm_tair_AMIP_exp_DJF=np.mean(dtexp[:,:,:],axis=2)[0,:]
       zm_tair_AMIP_diff_DJF=zm_tair_AMIP_exp_DJF - zm_tair_AMIP_ctl_DJF
    elif iss is "JJA":
       zm_tair_AMIP_ctl_JJA=np.mean(dtctl[:,:,:],axis=2)[0,:]
       zm_tair_AMIP_exp_JJA=np.mean(dtexp[:,:,:],axis=2)[0,:]
       zm_tair_AMIP_diff_JJA=zm_tair_AMIP_exp_JJA - zm_tair_AMIP_ctl_JJA


print(zm_tair_AMIP_diff_JJA)
#print(zm_ensmin_tair_diff)
#print(zm_ensmax_tair_diff)


zeros=np.zeros((nlat))

# make the plot
fig=plt.figure(figsize=(9,6))
ax1=fig.add_axes([0.10,0.15,0.3,0.4])

ax1.plot(lat[:],zm_ensmean_tair_diff_DJF[:],color="r",lw=2,ls="-",label="Scat - noScat, CMIP")
ax1.plot(lat[:],zm_tair_AMIP_diff_DJF[:],color="b",lw=2,ls="-",label="Scat - noScat, AMIP")
ax1.plot(lat[:],zeros[:],color="gray",lw=1,ls="-")
ax1.fill_between(lat[:],zm_ensmin_tair_diff_DJF[:],zm_ensmax_tair_diff_DJF[:],facecolor='r', alpha=0.3)
#ax1.plot(lat[:],means_exp[0,:],color="k",lw=2,ls=":",label="TSIS")
#ax1.legend(loc="upper left",fontsize=12)
ax1.legend(loc="upper left",fontsize=9)
ax1.set_title("DJF",fontsize=14)
ax1.set_ylabel("\u0394SAT (K)",fontsize=10)
ax1.set_xlabel("Latitude",fontsize=10)
ax1.set_xlim(-90,90)
ax1.set_ylim(-1.5,2.5)
xloc=[-90,-60,-30,0,30,60,90]
plt.xticks(xloc,fontsize=10)
plt.yticks(fontsize=10)

ax2=fig.add_axes([0.5,0.15,0.3,0.4])
ax2.plot(lat[:],zm_ensmean_tair_diff_JJA[:],color="r",lw=2,ls="-",label="Scat - noScat, CMIP")
ax2.plot(lat[:],zm_tair_AMIP_diff_JJA[:],color="b",lw=2,ls="-",label="Scat - noScat, AMIP")
ax2.plot(lat[:],zeros[:],color="gray",lw=1,ls="-")
ax2.fill_between(lat[:],zm_ensmin_tair_diff_JJA[:],zm_ensmax_tair_diff_JJA[:],facecolor='r', alpha=0.3)
ax2.legend(loc="upper left",fontsize=9)

ax2.set_title("JJA",fontsize=14) #+var_long_name,fontsize=12)
ax2.set_ylabel("\u0394SAT (K)",fontsize=10)
ax2.set_xlabel("Latitude",fontsize=10)
ax2.set_xlim(-90,90)
ax2.set_ylim(-1.5,2.5)
xloc=[-90,-60,-30,0,30,60,90]
plt.xticks(xloc,fontsize=10)
plt.yticks(fontsize=10)

#plt.savefig(figure_name+".eps")
plt.savefig("zm_tair_CMIP_AMIP_2ens_spread_2010-2019.png",dpi=(200))
plt.show()

exit()

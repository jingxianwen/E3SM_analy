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
ctl_name="Standard" #os.environ["ctl_name"]
exp_name="Modified_noEmis" #os.environ["exp_name"]
ctl_pref="E3SMv2_standard_PresSST_UMRadALLoff"
exp_pref="E3SMv2_modified_PresSST_UMRad_noEmis"

fpath_ctl="/global/cscratch1/sd/xianwen/acme_scratch/cori-knl/"+ctl_pref+"/"+"remap_180x360/"
fpath_exp="/global/cscratch1/sd/xianwen/acme_scratch/cori-knl/"+exp_pref+"/"+"remap_180x360/"
 
years=np.arange(2000,2011) 
months_all=["01","02","03","04","05","06","07","08","09","10","11","12"]

var_group_todo=2
# variable group 1:
#varnms=np.array(["ICEFRAC"])
varnms="T"
#pole='S'
#if pole is 'N':
#   var_long_name="Arctic Sea Ice Fraction"
#   figure_name="Arctic_Sea_Ice_Monthly"
#elif pole is 'S':
#   var_long_name="Antarctic Sea Ice Fraction"
#   figure_name="Antarctic_Sea_Ice_Monthly"
#units=" " #"Fraction"
#units=r"W/m$^2$"
units="K"
figure_name="time_series_"+varnms+"250mb_"+exp_name+"-"+ctl_name+".png"



# variable group 2:
#varnms=np.array(["FSSUS13","FSSUS12","FSSUS11","FSSUS10","FSSUS09",\
#        "FSSUS08","FSSUS07","FSSUS06","FSSUS05","FSSUS04",\
#        "FSSUS03","FSSUS02","FSSUS01","FSSUS14"])
#varnms_sub=np.array(["FSSDS13","FSSDS12","FSSDS11","FSSDS10","FSSDS09",\
#        "FSSDS08","FSSDS07","FSSDS06","FSSDS05","FSSDS04",\
#        "FSSDS03","FSSDS02","FSSDS01","FSSDS14"])
#var_long_name="Band-by-Band Surface net Upward SW"
#figure_name="Band_by_Band_surface_net_Upward_SW_ANN"
#units=r"W/m$^2$"

#f1=fpath_ctl+"solar_TSIS_cesm211_standard-ETEST-f19_g17-ens1.cam.h0.0001-01.nc"
#f2=fpath_exp+"tsis_ctl_cesm211_standard-ETEST-f19_g17-ens1.cam.h0.0001-01.nc"
nyrs=np.int64(years.size)
nmon=np.int64(12)
nall=nyrs*nmon
means_yby_ctl=np.zeros((nall)) #year by year mean for each variable
means_yby_exp=np.zeros((nall)) #year by year mean for each variable
#means_ctl=np.zeros((nall)) #multi-year mean for each variable
#means_exp=np.zeros(()) #multi-year mean for each variable
diffs=np.zeros((nall)) #multi-year exp-ctl diff for each variable
zeros=np.zeros(nall)
nn=np.int64(0)
for iy in range(0,years.size): 
    # read data and calculate mean/min/max
    for im in range(0,nmon):
        # open data file
        fctl=fpath_ctl+ctl_pref+".cam.h0."+str(years[iy])+"-"+months_all[im]+".nc"
        fexp=fpath_exp+exp_pref+".cam.h0."+str(years[iy])+"-"+months_all[im]+".nc"
        file_ctl=netcdf_dataset(fctl,"r")
        file_exp=netcdf_dataset(fexp,"r")
        
        # read lat and lon
        lat=file_ctl.variables["lat"]
        lon=file_ctl.variables["lon"]
        lev=file_ctl.variables["lev"]
        if var_group_todo is 1:
           dtctl=file_ctl.variables[varnms]
           dtexp=file_exp.variables[varnms] 
        elif var_group_todo is 2:
           lev250=np.min(np.where(lev[:]>250.))
           dtctl=file_ctl.variables[varnms][:,lev250,:,:] #*1000. #[:,:,:]-file_ctl.variables[varnms_sub][:,:,:]
           dtexp=file_exp.variables[varnms][:,lev250,:,:] #*1000. #[:,:,:]-file_exp.variables[varnms_sub][:,:,:]
        nlat=180
        #if pole == "N":
        #   latbound1=np.min(np.where(lat[:]>50))
        #   latbound2=nlat
        #elif pole == "S":
        #   latbound1=0
        #   latbound2=np.max(np.where(lat[:]<-50))+1
        #print(lat[latbound1:latbound2])
        #exit()
        #dtdif=dtexp[:,:,:]-dtctl[:,:,:]
        means_yby_ctl[nn]=get_area_mean_min_max(dtctl[:,:,:],lat[:])[0]
        means_yby_exp[nn]=get_area_mean_min_max(dtexp[:,:,:],lat[:])[0]
        nn=nn+1
        #stats_dif[i]=get_area_mean_min_max(dtdif[:,:,:],lat[:])[0]
        #stats_difp[i]=stats_dif[0]/stats_ctl[0]*100.

diffs=means_yby_exp-means_yby_ctl
xlocs=np.arange(1,nall+1) 


# make the plot
fig=plt.figure(figsize=(7,8))
ax1=fig.add_axes([0.18,0.57,0.78,0.35])
#ax2=fig.add_axes([0.13,0.10,0.78,0.35])
#x=np.array([1,2,3,4,5,6,7,8,9,10,11,12])
#labels_fig=np.array(["J","F","M","A","M","J","J","A","S","O","N","D"])
ax1.plot(xlocs[:],diffs[:],color="k",lw=2)
#ax1.set_title("Diff in "+varnms+" ("+exp_name+"-"+ctl_name+")",fontsize=12)
ax1.set_title("Diff in "+varnms+" 250mb ("+exp_name+"-"+ctl_name+")",fontsize=12)
ax1.set_ylabel(units,fontsize=12)
ax1.set_xlabel("months",fontsize=12)
#ax1.grid(True)
#ax1.set_axisbelow(True)
#ax1.xaxis.grid(color='gray', linestyle=':')
#ax1.yaxis.grid(color='gray', linestyle=':')
#ax1.set_xticks(xlocs)
#ax1.set_xticklabels(labels=labels_fig,rotation=0,fontsize=10)
ax1.set_xlim(1,nn)
#ax1.set_ylim=([0,means_ctl.max*1.1])

#bars=[None]*diffs_sig.size
#ax2.plot(x[:],diffs[:],color="tab:blue")
#ax2.plot(x[:],diffs[:],color="k",lw=2)
#ax2.plot(x[:],diffs_sig[:],color="darkorange",lw=4,alpha=1.0)
#ax2.plot(x[:],zeros[:],color="gray",lw=1)
#ax2.set_title("Diff in "+var_long_name+" (ANN)",fontsize=12)
#ax2.set_ylabel(units,fontsize=12)
#ax2.set_xlabel("month",fontsize=12)
#ax2.grid(True)
#ax2.set_axisbelow(True)
#ax2.xaxis.grid(color='gray', linestyle=':')
#ax2.yaxis.grid(color='gray', linestyle=':')
#ax2.set_xticks(x)
#ax2.set_xticklabels(labels=labels_fig,rotation=0,fontsize=10)
#ax2.set_xlim(1,12)
plt.savefig(figure_name+".png")
plt.show()

exit()

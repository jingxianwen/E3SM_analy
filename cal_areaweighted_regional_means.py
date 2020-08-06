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
from get_parameters import get_area_mean_min_max

# data path
fpath_ctl='/global/cscratch1/sd/xianwen/E3SM_simulations/CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl/archive/climo/'
f1=fpath_ctl+"CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl_climo_JJA.nc"
print(f1)

varlst=[ \
        "SOLIN","FSNTOA","FSNTOAC","FLUT","FLUTC", \
        "FSDS","FSDSC","FSNS","FSNSC","FLDS","FLDSC","FLNS","FLNSC","LHFLX","SHFLX",\
        "CLDLOW","CLDMED","CLDHGH","TGCLDLWP","TGCLDIWP"]

# open data file
file_ctl=netcdf_dataset(f1,"r")

# read lat and lon
lat=file_ctl.variables["lat"]
lon=file_ctl.variables["lon"]
#lev=file_ctl.variables["lev"]
#lev500=np.min(np.where(lev[:]>500.))
lat_N=np.min(np.where(lat[:]>66.5))
for var in varlst:
    dtctl=file_ctl.variables[var][:,:,:] #[time,lat,lon]
    stat=get_area_mean_min_max(dtctl[:,lat_N:,:],lat[lat_N:])
    #print(var,round(stat[0],2))
    print(var,stat[0])

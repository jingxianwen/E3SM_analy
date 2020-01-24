
#import netCDF4
from netCDF4 import Dataset as netcdf_dataset
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

dataset=netcdf_dataset("test.nc","r")
flut=dataset.variables["FLUT"][0,:,:]
lat=dataset.variables["lat"][:]
lon=dataset.variables["lon"][:]

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
plt.contourf(lon, lat, flut, 60,transform=ccrs.PlateCarree())
# Save the plot by calling plt.savefig() BEFORE plt.show()
#plt.savefig('coastlines.pdf')
#plt.savefig('coastlines.png')

plt.show()

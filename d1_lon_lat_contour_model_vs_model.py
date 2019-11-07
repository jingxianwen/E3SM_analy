
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

def lon_lat_contour_model_vs_model(varnm,season,scale_ctl,scale_exp):
    # data path
    ctl_name=os.environ["ctl_name"]
    exp_name=os.environ["exp_name"]
    fpath_ctl=os.environ["fpath_ctl"]
    fpath_exp=os.environ["fpath_exp"]
    
    # open data file
    file_ctl=netcdf_dataset(fpath_ctl,"r")
    file_exp=netcdf_dataset(fpath_exp,"r")

    # read lat and lon
    lat=file_ctl.variables["lat"]
    lon=file_ctl.variables["lon"]

    # read data and calculate mean/min/max
    dtctl=file_ctl.variables[varnm] #*scale_ctl
    dtexp=file_exp.variables[varnm] #*scale_exp
    dtdif=dtexp[:,:,:]-dtctl[:,:,:]
    stats_ctl=get_area_mean_range(dtctl[:,:,:],lat[:])
    stats_exp=get_area_mean_range(dtexp[:,:,:],lat[:])
    stats_dif=get_area_mean_range(dtdif[:,:,:],lat[:])
    
    # add cyclic
    dtctl=add_cyclic_point(dtctl[:,:,:])
    dtexp=add_cyclic_point(dtexp[:,:,:])
    dtdif=add_cyclic_point(dtdif[:,:,:])
    lon=np.append(lon[:],360.)

    # make plot
    parameters=get_parameters(varnm,season)
    projection = ccrs.PlateCarree(central_longitude=180)

    fig = plt.figure(figsize=[7.0,11.0],dpi=150.)
    plotTitle = {'fontsize': 13.}
    plotSideTitle = {'fontsize': 9.}
    plotText = {'fontsize': 8.}
    panel = [(0.1691, 0.6810, 0.6465, 0.2258), \
             (0.1691, 0.3961, 0.6465, 0.2258), \
             (0.1691, 0.1112, 0.6465, 0.2258), \
             ]
    labels=[exp_name,ctl_name,exp_name+"-"+ctl_name] 
    units=parameters["units"]
    for i in range(0,3):
       #1. first plot
        levels = None
        norm = None
        if i != 2:
            cnlevels=parameters["contour_levs"]
        else:
            cnlevels=parameters["diff_levs"]

        #if len(cnlevels) >0:
        #        levels = [-1.0e8] + cnlevels + [1.0e8]
        #        norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)

        ax = fig.add_axes(panel[i],projection=projection)
        ax.set_global()
        cmap=parameters["colormap"]
        #p1 = ax.contourf(lon[:],lat[:],dtexp[0,:,:])
        if i == 0:
            dtplot=dtexp[:,:,:]
            cmap=parameters["colormap"]
            stats=stats_exp
        elif i == 1:
            dtplot=dtctl[:,:,:]
            cmap=parameters["colormap"]
            stats=stats_ctl
        else:
            dtplot=dtdif[:,:,:]
            cmap=parameters["colormap_diff"]
            stats=stats_dif

        p1 = ax.contourf(lon[:],lat[:],dtplot[0,:,:],\
                    transform=projection,\
                    #norm=norm,\
                    levels=cnlevels,\
                    cmap=cmap,\
                    extend="both",\
            	    )
        ax.set_aspect("auto")
        ax.coastlines(lw=0.3)
        # title
        ax.set_title(labels[i],loc="left",fontdict=plotSideTitle)
        #ax.set_title("exp",fontdict=plotTitle)
        ax.set_title(units,loc="right",fontdict=plotSideTitle)
        ax.set_xticks([0, 60, 120, 180, 240, 300, 359.99], crs=ccrs.PlateCarree())
        ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
        lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='.0f')
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.tick_params(labelsize=8.0, direction='out', width=1)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        # color bar
        cbax = fig.add_axes((panel[i][0] + 0.6635, panel[i][1] + 0.0215, 0.0326, 0.1792))
        #cbax = fig.add_axes((panel[i][0] + 0.6635, panel[i][1] + 0.0215, 0.0326, 0.2850))
        cbar = fig.colorbar(p1, cax=cbax, ticks=cnlevels)
        #w, h = get_ax_size(fig, cbax)
        cbar.ax.tick_params(labelsize=9.0, length=0)

        # Mean, Min, Max
        fig.text(panel[i][0] + 0.6635, panel[i][1] + 0.2107,
                 "Mean\nMin\nMax", ha='left', fontdict=plotText)
        fig.text(panel[i][0] + 0.7835, panel[i][1] + 0.2107, "%.2f\n%.2f\n%.2f" %
                 stats[0:3], ha='right', fontdict=plotText)

    #fig.suptitle(varnm, x=0.5, y=0.96, fontsize=14)
    fig.suptitle(varnm, x=0.5, y=0.96, fontdict=plotTitle)
    plt.show()
    plt.close()
    
    
def get_parameters(varnm,season):
    #list_rad=["FLUT","FLUTC","FLNT","FLNTC","FSNT","FSNTC","FSDS","FSDSC","FSNS","FSNSC"]
    if varnm == "FLUT":
        parameters={"units":"W/m2",\
		   "contour_levs":[120, 140, 160, 180, 200, 220, 240, 260, 280, 300],\
		   "diff_levs":[-50, -40, -30, -20, -10, -5, 5, 10, 20, 30, 40, 50],\
                   "colormap":"PiYG_r",\
                   "colormap_diff":"bwr"\
		   }
    return parameters

def get_area_mean_range(varnm,lat):
   # 1. area weighted average 
    #convert latitude to radians
    latr=np.deg2rad(lat)
    #use cosine of latitudes as weights for the mean
    weights=np.cos(latr)
    #first calculate zonal mean
    zonal_mean=varnm.mean(axis=2)
    #then calculate weighted global mean
    area_mean=np.average(zonal_mean,axis=1,weights=weights)
   # 2. min and max
    minval=varnm.min()
    maxval=varnm.max()
    return area_mean,minval,maxval


#Program to get a wrf output file with u10,v10,hgt and get the quantity U dot nabla b.

#Load required libraries
from netCDF4 import Dataset
from wrf import to_np, getvar, smooth2d, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
#Set filename
ncfile = Dataset("../dataFiles/wrfout_d03_2017-02-24_00:00:00")
#Open new file for result
newfile = Dataset("netcdf_div/DIV-wrfout_d03_2017-02-24_00:00:00.nc","w",format="NETCDF4")    

n_tim = len(ncfile.dimensions["Time"])
n_lat = len(ncfile.dimensions["south_north"])
n_lon = len(ncfile.dimensions["west_east"])

#u10 = getvar(ncfile,'U10',timeidx=0)
#v10 = getvar(ncfile,'V10',timeidx=0)
b = getvar(ncfile,"ter")
hgt = getvar(ncfile,'HGT')

gradb = np.gradient(b)
dx = gradb[0]
dy = gradb[1]

#Get variables for future printing
xlat = getvar(ncfile,'lat',meta=False)
xlong = getvar(ncfile,'lon',meta=False)
xTimes = getvar(ncfile,'Times')

time = newfile.createDimension("time", None)
lat = newfile.createDimension("south_north", n_lat)
lon = newfile.createDimension("west_east", n_lon) 
# Add variables for dimensions
times = newfile.createVariable("time","f4",("time",))
latitudes = newfile.createVariable("lat","f4",("south_north","west_east"))
longitudes = newfile.createVariable("lon","f4",("south_north","west_east"))
div = newfile.createVariable("UdotGradB","f4",("time","south_north","west_east"))


#Put the actual data
latitudes[:,:] = xlat
longitudes[:,:] = xlong


#Begin time loop
timevec=range(0,n_tim)
#timevec=range(0,1)

for tim in timevec:

    u = getvar(ncfile,'U10',meta=False,timeidx=tim)
    v = getvar(ncfile,'V10',meta=False,timeidx=tim)

    #Now make vectors of what we want to operate on, using reshape
    dimA = u.shape[0]
    dimB = u.shape[1]
    dimC = dimA * dimB
    urs = np.reshape(u,dimC)
    vrs = np.reshape(v,dimC)
    dxrs = np.reshape(dx,dimC)
    dyrs = np.reshape(dy,dimC)

    vel = np.stack((urs,vrs),axis=-1)
    grb = np.stack((dxrs,dyrs),axis=-1)

    vdotgb = []
    for line in range(len(vel)):
        dotp = np.dot(vel[line],grb[line])
        vdotgb.append(dotp)
        
    vdotb = np.reshape(vdotgb,[dimA,dimB])


    # PLOTTING SECTION

    lats, lons = latlon_coords(b)
    cart_proj = get_cartopy(b)

    # Create a figure
    fig = plt.figure(figsize=(12,9))
    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj) 

    # Download and add the states and coastlines
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                name='admin_1_states_provinces_shp')
    ax.add_feature(states, linewidth=.5)
    ax.coastlines('50m', linewidth=0.8)

    # Make the contour outlines and filled contours for the smoothed sea level pressure.
    plt.contour(to_np(lons), to_np(lats), vdotb, 10, colors="black",
                transform=crs.PlateCarree())
    plt.contourf(to_np(lons), to_np(lats), vdotb, 10, transform=crs.PlateCarree(),
                cmap=get_cmap("jet"))

    # Add a color bar
    plt.colorbar(ax=ax, shrink=.62)

    # Set the map limits.  Not really necessary, but used for demonstration.
    ax.set_xlim(cartopy_xlim(b))
    ax.set_ylim(cartopy_ylim(b))

    # Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")

    plt.title(r'$U \cdot \nabla b$, tiempo {}Z'.format(tim))
    
    figstr = '../figures/udotnb_time_{}'.format(tim+100)
    plt.savefig(figstr)
    plt.clf()

    #plt.show()

    #Now lets write the resulting data into a NetCDF file for future use. 


    times[tim] = tim*3600
    div[tim,:,:] = vdotb
    
newfile.close()

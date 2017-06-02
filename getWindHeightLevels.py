#===========================================================================================
#Program Name: getWindHeightLevels.py

#Description: Program that takes WRF output as input, and then transforms wind variable from eta    
#             coordinates to height coordinates. 
#Language: Python 2.7

#Programmer Oscar Jurado (ojurado@ciencias.unam.mx)
#Creation date: June-2017

#------------------Use------------------------------------------------------------------
#anaconda
#source activate pyWRF
#python getWindHeightLevels.py

#------------------Requisites------------------------------------------------------------
#WRF output data in ../../../salidas/

#-----------------Version---------------------------------------------------------------
#v1.0 June/17 Program is created

#----------------Known issues-----------------------------------------------------------

#----------------Dependencies and Libraries----------------------------------------------
#wrf-python package. 
#Python 2.7
#numpy (1.9 or later)
#wrapt (1.10 or later)
#xarray (0.7.0 or later)
#PyNIO (1.4.3 or later)
#netCDF4-python (1.2.0 or later)
#PyNGL (1.4.3 or later)
#matplotlib (1.4.3 or later)
#cartopy (0.13 or later)
#basemap (1.0.8 or later)
#csv

#It is highly recommended to create a conda environment with this packages. An environment.yml with all the requisites is available from the author.

#---------------Import-------------------------------------------------------------------
from __future__ import print_function
from netCDF4 import Dataset
from wrf import getvar, interplevel, to_np, vinterp, ll_to_xy
from math import pi, atan2
import numpy as np
import csv

#-----------------Local directories----------------------------------------------------- 
file_dir = '../../../salidas/'
latlon_dir = '/home/esli/Documentos/WRF/experimentos/Toluca/pluma'

#=====================================================================
#       FUNCTION = uv2sd
#=====================================================================
#This function takes u and v as 2d arrays and returns two 2d arrays, one for the wind speed
#and another for the wind direction. 
def uv2sd(u,v):
    "Transform u,v to wind speed and direction."
    wind_mag = np.sqrt(u**2 + v**2)
    wind_dir = np.zeros(u.shape)

    for i in range(u.shape[0]):
        for j in range(u.shape[1]):
            wind_dir[i][j] = 270 - ((180/pi) * atan2(v[i][j],u[i][j]))

    return wind_mag, wind_dir

#======================================================================
#       FUNCTION = getLatLon
#======================================================================
#This function takes input from a text file to read lat/lon points for the stations.
def getLatLon(filename):
    "Read lat and lon from filename."
    geofile = open(filename , "r")
    lines = geofile.readlines()
    geofile.close()
    #Remove header
    lines = lines[1:-1]
    #Make empty arrays for lat and lon
    lats = []
    lons = []
    #Get only 3d and 4th column, making them float numbers. 
    for line in lines:
        lats.append(float(line.split('\t')[2]))
        lons.append(float(line.split('\t')[3]))
        
    return lats, lons
    
#======================================================================
#       FUNCTION = writeOutput
#======================================================================
#This function writes the output file
def writeOutput(lat,lon,levels,hour1,wmag1,wdir1,hour2,wmag2,wdir2,hour3,wmag3,wdir3,locat):
    "Writes output to file, lat lon are scalars, wmag wdir are vectors"
    filename = './latlonpairs/{}.csv'.format(locat+100)
    wrf_points = open(filename,"w")
    mywriter = csv.writer(wrf_points, delimiter=' ')
    latstr = 'Lat={}'.format(lat)
    lonstr = 'Long={}'.format(lon)
    ll_line = [latstr,lonstr]
    header = ["Altitude","Hour","WD[deg]","WS[m/s]","Hour","WD[deg]","WS[m/s]","Hour","WD[deg]","WS[m/s]"]
    blank = ''
    mywriter.writerow(ll_line)
    mywriter.writerow(blank)
    mywriter.writerow(header)
    mywriter.writerow(blank)
    for indx in range(len(wmag1)):
        row = []
        row.append(levels[indx])
        row.append(hour1)
        row.append('{0:0.2f}'.format(wdir1[indx]))
        row.append('{0:0.2f}'.format(wmag1[indx]))
        row.append(hour2)
        row.append('{0:0.2f}'.format(wdir2[indx]))
        row.append('{0:0.2f}'.format(wmag2[indx]))
        row.append(hour3)
        row.append('{0:0.2f}'.format(wdir3[indx]))
        row.append('{0:0.2f}'.format(wmag3[indx]))
        mywriter.writerow(row)
        mywriter.writerow(blank)
        
    wrf_points.close()
    
    
#-----------------BEGIN PROGRAM------------------------------------------------------------

#First open input WRF-output file
#Make string with file name and then assing to nc file

filestr = '{}times.nc'.format(file_dir)
ncfile = Dataset(filestr)

#Read Lat/Lon data too
#Make string for files
geostr = '{}/Toluca_Vis10_20170119.txt'.format(latlon_dir)
#Call function
lats, lons = getLatLon(geostr)

#Also open the output file


#Get the model pressure
prs = getvar(ncfile,"pressure")

#Get the model height (in masl)
HGT = getvar(ncfile,"z")
hgt = to_np(HGT)
#Get total number of levels
levels = hgt.shape[0]

#Now get terrain height
THT = getvar(ncfile,"HGT")
tht = to_np(THT)

#Convert model height to meters above surface, we need a loop cause hgt has several levels
zht = np.zeros(hgt.shape)
for i in range(levels):
    zht[i] = hgt[i][:][:] - tht

#Get PBLH data for each timestep
PBL = getvar(ncfile,"PBLH",timeidx=0)
pbl = to_np(PBL)
    
#Get the max PBLH for the field.
pbl_point = np.max(pbl)

#Create the levels for interpolation, using pblh of station as upper limit
interp_levels_m = range(10,pbl_point,50)
#interp_field function expects the levels to be in Km, so we need to convert
interp_levels = [float(x) / 1000 for x in interp_levels_m]
    
    
#Time loop
timevec=[0,1,2]
hour = [18,19,20]
uint = np.zeros((len(timevec),len(interp_levels),30,39))
vint = np.zeros((len(timevec),len(interp_levels),30,39))

for tim in timevec:

    #Unstagger winds to mass points
    uwn = getvar(ncfile,"ua",timeidx=tim)
    vwn = getvar(ncfile,"va",timeidx=tim)



    #Now do the actual vertical interpolation.
    #U component
    interp_u = vinterp(ncfile,
                field=uwn,
                vert_coord="ght_agl",
                interp_levels=interp_levels,
                extrapolate=False,
                field_type="none",
                timeidx=tim,
                log_p=False)

    #V component
    interp_v = vinterp(ncfile,
                field=vwn,
                vert_coord="ght_agl",
                interp_levels=interp_levels,
                extrapolate=False,
                field_type="none",
                timeidx=tim,
                log_p=False)

    #Get the wind profile for the desired station

    uint[tim] = to_np(interp_u)
    vint[tim] = to_np(interp_v)

    
    

#Iterate over every lat and lon value
    
for locat in range(len(lats)):
        
    #Declare the lat/lon for the desired point
    lat = lats[locat]
    lon = lons[locat]
        
    #Get the corresponding domain index for the point
    X_Y = ll_to_xy(ncfile, lat, lon)
    x_y = to_np(X_Y)
    xindx = x_y[0]
    yindx = x_y[1]
    if (xindx > 29):
        xindx = 29
        
    if (yindx > 29):
        yinx = 29
    
    ws1_pnt = []
    wd1_pnt = []
    ws2_pnt = []
    wd2_pnt = []
    ws3_pnt = []
    wd3_pnt = []
        

    #Loop in heights 
    for i in range(uint.shape[1]):
        u1 = uint[0][i][:][:]
        v1 = vint[0][i][:][:]
        #Transform from U V to speed and direction (meteorological)
        ws1, wd1 = uv2sd(u1,v1)
        ws1_pnt.append(ws1[xindx][yindx])
        wd1_pnt.append(wd1[xindx][yindx])
        #Time2
        u2 = uint[1][i][:][:]
        v2 = vint[1][i][:][:]
        #Transform from U V to speed and direction (meteorological)
        ws2, wd2 = uv2sd(u2,v2)
        ws2_pnt.append(ws2[xindx][yindx])
        wd2_pnt.append(wd2[xindx][yindx])
        #Time3
        u3 = uint[2][i][:][:]
        v3 = vint[2][i][:][:]
        #Transform from U V to speed and direction (meteorological)
        ws3, wd3 = uv2sd(u3,v3)
        ws3_pnt.append(ws3[xindx][yindx])
        wd3_pnt.append(wd3[xindx][yindx])
            
        #Write the output file
        
        writeOutput(lat,lon,interp_levels_m,hour[0],ws1_pnt,wd1_pnt,hour[1],ws2_pnt,wd2_pnt,hour[2],ws3_pnt,wd3_pnt,locat)










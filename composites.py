# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 14:14:32 2017

@author: Ambidados
"""
# import necessary libraries and rename them so they are easier to use
import numpy as np
#import netCDF4 as nc
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
from netCDF4 import Dataset

# create a Dataset object called Dset using the Dataset constructor from the
# netCDF4 library
# To use this program simply replace <netCDF_File> with the name of your file.
Dset = Dataset("20110106090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc", "r", format ="NETCDF4")

# Get information from our Dataset object, Dset and store it in varialbes
# my netCDF file contains information on global sea surface temperature. 
# Dset.variables calls my Dataset object and goes to the variables section of
# the information that it holds. We find the names of the variables that we need
# and assign them to variables in python. The array that stores all of the SST
# information is called “analysed_sst” and our first line gets that data and
# stores it in the variable sst.
#%%
sst_m = comp_ng
sst_units = Dset.variables["analysed_sst"].units

# Here we do the same thing except now we are getting information about the lon and lat variables.
lon_scale = Dset.variables["lon"][:] 
lon_units = Dset.variables["lon"].units
lon_max = Dset.variables["lon"].valid_max 
lon_min = Dset.variables["lon"].valid_min

lat_scale = Dset.variables["lat"][:] 
lat_units = Dset.variables["lat"].units
lat_max = Dset.variables["lat"].valid_max 
lat_min = Dset.variables["lat"].valid_min
plt.figure()
#____________________
# These lines create the basemap for the data.
Mbase = bm.Basemap(projection = 'cyl', llcrnrlat = -24.0, llcrnrlon = -44.0, urcrnrlat = -21.0, urcrnrlon = -40.0, resolution='f')

Mbase.drawparallels(np.array([-24, -23, -22, -21]), labels=[0,1,0,0])
Mbase.drawmeridians(np.array([-44, -43, -42, -41, -40]), labels = [0,0,1,0])

#Mbase.drawsmask(land_color='0.8', ocean_color='w', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution='l', grid=5)
Mbase.drawcoastlines()
Mbase.fillcontinents(color = 'gray') # blue can be changed to any color that you would like
#____________________

mymap = plt.contour(lon_scale, lat_scale, sst_m, 15,levels=[16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.], colors = 'k')
mymap2 = plt.contourf(lon_scale, lat_scale, sst_m, 15,levels=[16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.], cmap =  plt.cm.jet)
plt.colorbar(mymap2, orientation = 'horizontal') 
plt.plot(-41.78,-22.37,'.r', markersize=10,markeredgewidth=2)
plt.plot(-43.19,-22.99,'.b', markersize=10,markeredgewidth=2)
plt.axis([-44,-40,-24,-21])
plt.xlabel(u'Longitude (°)')
plt.ylabel(u'Latitude (°)')
#plt.title(u'Água Fria')
plt.text(-43.4, -22.95, u'RJ')
plt.text(-42.4, -22.4, u'Macaé')
plt.text(-43.9, -21.5, u'Composite TSM - (abaixo de 22.5°)')
plt.show()
plt.savefig('figuras/composite_ng.png')

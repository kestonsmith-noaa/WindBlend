


import numpy as np
from scipy.spatial import KDTree

def latlon_to_cartesian(lat, lon, R=6371):
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    x = R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return np.stack([x, y, z], axis=-1)

# 1. Define your data (Pointset and target Point)
pointset_latlon = np.array([[34.05, -118.24], [40.71, -74.00]]) # LA, NY
target_latlon = np.array([37.77, -122.41]) # SF

# 2. Convert to 3D
points_3d = latlon_to_cartesian(pointset_latlon[:,0], pointset_latlon[:,1])
target_3d = latlon_to_cartesian(target_latlon[0], target_latlon[1])

# 3. Build tree and query nearest
tree = KDTree(points_3d)
euclidean_dist, index = tree.query(target_3d)

# 4. Convert Euclidean chord distance to Great Circle distance
# Formula: d_gc = 2 * R * arcsin(chord_dist / (2 * R))
R = 6371
great_circle_dist = 2 * R * np.arcsin(euclidean_dist / (2 * R))

print(f"Minimum distance: {great_circle_dist:.2f} km")



import numpy as np

def haversine_min_distance(origin_lat, origin_lon, lats_array, lons_array):
    # Earth radius in km
    R = 6371.0 
    
    # Convert all coordinates from degrees to radians
    lat1, lon1 = np.radians(origin_lat), np.radians(origin_lon)
    lats2, lons2 = np.radians(lats_array), np.radians(lons_array)
    
    # Differences
    dlat = lats2 - lat1
    dlon = lons2 - lon1
    
    # Haversine formula
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lats2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    distances = R * c
    return np.min(distances)

def haversine_min_distanceX(lat1, lon1, lats2, lons2):
    # Earth radius in km
    R = 6371.0 
    
    # Convert all coordinates from degrees to radians
    #DONE BEFORE INPUT
#    lat1, lon1 = np.radians(origin_lat), np.radians(origin_lon)
#    lats2, lons2 = np.radians(lats_array), np.radians(lons_array)
    
    # Differences
    dlat = lats2 - lat1
    dlon = lons2 - lon1
    
    # Haversine formula
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lats2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    distances = R * c
    return np.min(distances)

def QuickDistance(lat1, lon1, lats2, lons2):
    deg2kmY=111.
    deg2kmX=np.cos( np.pi * lat1 / 180.)
    d= np.min(  np.sqrt( (  (lat1-lats2)*deg2kmY)**2 + ((lon1-lons2)*deg2kmX)**2 )  )
    return d

# Example: Distance from one point to a set of 3 others
origin = (40.7128, -74.0060) # New York City
lats = np.array([34.0522, 51.5074, 48.8566]) # LA, London, Paris
lons = np.array([-118.2437, -0.1278, 2.3522])

min_km = haversine_min_distance(origin[0], origin[1], lats, lons)
print(f"Minimum Distance: {min_km:.2f} km")

import math

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0 # Radius of Earth in km
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dlat, dlon = math.radians(lat2 - lat1), math.radians(lon2 - lon1)
    a = math.sin(dlat / 2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlon / 2)**2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

import datetime
import numpy as np
import netCDF4 as nc
import sys
import re
import InterpNWPSUtility as nwps

import xarray as xr
import numpy as np
import xesmf as xe

flin=sys.argv[1]
flout=sys.argv[2]

mshfl="meshes/RWPS.V0a.msh"
xi, yi, ei = nwps.loadWW3Mesh(mshfl)
nn=len(xi)
xi=xi+360

data = nc.Dataset(flin,"r")
t1=np.asarray(data["time"][:])

nn=len(xi)
nt=len(t1)

x1=np.asarray(data["longitude"][:])
y1=np.asarray(data["latitude"][:])
print("x1,y1 shape")
nx=x1.shape[0]
ny=x1.shape[1]
print(x1.shape)
print(y1.shape)
print(x1)

##################################################################################
# Compute distance to boundary for each node in mesh:
# only needs to be done for nodes within interpolater
##################################################################################

xb=np.hstack((x1[1,:],x1[:,ny-1].T,x1[nx-1,:],x1[:,1].T))
yb=np.hstack((y1[1,:],y1[:,ny-1].T,y1[nx-1,:],y1[:,1].T))

#convert to radians for more efficient distance calculation
ybr, xbr = np.radians(yb), np.radians(xb)
yir, xir = np.radians(yi), np.radians(xi)

np.savetxt('xbyb.txt', np.vstack((xb,yb)))
np.savetxt('xiyi.txt', np.vstack((xi,yi)))

print("xb")
print(xb)
print("xb shape")
print(xb.shape)
print("xi shape")
print(xi.shape)

dist2bnd=np.zeros(nn)
for k in range(nn):
#    dist2bnd[k]=haversine_min_distance(yi[k],xi[k],np.array(yb),np.array(xb))
#    dist2bnd[k]=haversine_min_distanceX(yir[k],xir[k],np.array(ybr),np.array(xbr))
    dist2bnd[k]=QuickDistance(yi[k],xi[k],yb,xb)
    if k%10000==0:
        print(str(k)+":"+ str(nn)+":"+str(k/nn) )
        print(dist2bnd[k])
##################################################################################


##################################################################################
# Interpolate from curvilinear mesh to unstructured mesh nodes
##################################################################################

# 1. Load your RRFS data
# RRFS uses Lambert Conformal or curvilinear grids; xESMF handles these via 2D lat/lon arrays
ds_source = xr.open_dataset(flin)

# 2. Define target points (e.g., specific weather stations)
target_lats = yi
target_lons = xi

# Create a target dataset formatted for 'locstream'
ds_target = xr.Dataset(
    {
        "lat": (["locations"], target_lats),
        "lon": (["locations"], target_lons),
    }
)

# 3. Create the Regridder
# 'bilinear' is typically used for continuous variables like wind speed
# locstream_out=True tells xESMF to interpolate to the specific points in ds_target
regridder = xe.Regridder(
    ds_source, 
    ds_target, 
    method="bilinear", 
    locstream_out=True,
    unmapped_to_nan=True
)

# 4. Perform the interpolation
# If interpolating wind direction, convert to U/V components first to avoid angular issues
ds_interpolated = regridder(ds_source)

u=ds_interpolated["UGRD_10maboveground"].values
v=ds_interpolated["VGRD_10maboveground"].values
##################################################################################

fill_value0=-99999
ne=ei.shape[0]

with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:

    ncout.createDimension('level' , 1)  
    ncout.createDimension('node' , nn)
    ncout.createDimension('element' , ne)
    ncout.createDimension('time', nt)
    ncout.createDimension('noel', 3)

    lon_var=ncout.createVariable('longitude', 'f8', ('node',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'X'
    lon_var[:]=xi[:]

    lat_var=ncout.createVariable('latitude', 'f8', ('node',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'Y'
    lat_var[:]=yi[:]

    time_var=ncout.createVariable('time', 'f8', ('time',))
    time_var.units         = 'seconds since 1970-01-01 00:00:00.0 0:00'
    time_var.long_name     = 'verification time generated by wgrib2 function verftime()'
    time_var.standard_name = 'time'
    time_var.axis          = 'T'
    time_var.reference_time = 1777334400
    time_var.reference_date = '2026.04.28 00:00:00 UTC'
    time_var[:]=t1[:]

    tri_var=ncout.createVariable('tri', 'i4', ('noel','element'))
    tri_var.long_name     = 'element list'
    tri_var.standard_name = 'element list'
    tri_var[:]=np.transpose(ei)

    d_var=ncout.createVariable('dist2bnd', 'f4', ('node',))
    d_var.long_name     = 'distance to boundary'
    d_var.units         = 'km'
    d_var.standard_name = 'distance to boundary'
    d_var[:]=dist2bnd[:]

    u_var=ncout.createVariable('uwnd', 'f4', ('node','time'),fill_value    = fill_value0)
    u_var.long_name     = 'eastward_wind'
    u_var.units         = 'm/s'
    u_var.standard_name = 'eastward_wind'
    u_var.level = '10 m above ground'
    u_var[:,:]=np.transpose(u[:,:])

    v_var=ncout.createVariable('vwnd', 'f4', ('node','time'),fill_value    = fill_value0)
    v_var.long_name     = 'northward_wind'
    v_var.units         = 'm/s'
    v_var.standard_name = 'northward_wind'
    v_var.level = '10 m above ground'
    v_var[:,:]=np.transpose(v[:,:])

    ncout.close

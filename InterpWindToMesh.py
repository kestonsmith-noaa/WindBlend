


import numpy as np
import os

def QuickDistance(lat1, lon1, lats2, lons2):
    deg2kmY=111.
    deg2kmX=np.cos( np.pi * lat1 / 180.)*deg2kmY
#    deg2kmX=np.cos( lat1 )*deg2kmY
    d= np.min(  np.sqrt( (  (lat1-lats2)*deg2kmY)**2 + ((lon1-lons2)*deg2kmX)**2 )  )
    return d

import datetime
import netCDF4 as nc
import sys
import re
import InterpNWPSUtility as nwps

import xarray as xr
import numpy as np
#import xesmf as xe
import esmpy

import scipy.sparse as sp

def VarianceLinearDistanceToBndy(InteriorNodeList, DistanceToBoundary, InteriorVariance, VarianceOnBoundary, LengthScale):
    Variance=np.zeros(len(DistanceToBoundary))+np.inf
    SpatialFunction=DistanceToBoundary/LengthScale
    j=np.where(SpatialFunction>1.)
    SpatialFunction[j]=1.
    Variance[InteriorNodeList] = VarianceOnBoundary + ( InteriorVariance - VarianceOnBoundary ) * SpatialFunction[InteriorNodeList]
    return Variance

def VarianceInverseDistanceToBndy(InteriorNodeList, DistanceToBoundary, InteriorVariance, LengthScale):
    Variance=np.zeros(len(DistanceToBoundary))+np.inf
    SpatialFunction=LengthScale / DistanceToBoundary
    j=np.where(SpatialFunction>1.)
    SpatialFunction[j]=1.
    Variance[InteriorNodeList] = InteriorVariance  * SpatialFunction[InteriorNodeList]
    return Variance


def CurvilinearGridCreateInterpWeights(xi,yi,x1,y1, weights_file):
# Compute interpolation weights to interpolate from curvilinear grid (x1,y1) to points (xi,yi)
    nx,ny=x1.shape
    nn=len(xi)
    n1=nx*ny
    src_lon = x1
    src_lat = y1
    dst_lon=np.zeros((1,nn))
    dst_lat=np.zeros((1,nn))
    dst_lon[0,:] = xi[:]
    dst_lat[0,:] = yi[:]

    esmpy.Manager()

    src_grid = esmpy.Grid(
      max_index=np.array([src_lon.shape[0], src_lon.shape[1]]),
      staggerloc=esmpy.StaggerLoc.CENTER,
      coord_sys=esmpy.CoordSys.SPH_DEG
    )

    dst_grid = esmpy.Grid(
      max_index=np.array([nn, 1]),
      staggerloc=esmpy.StaggerLoc.CENTER,
      coord_sys=esmpy.CoordSys.SPH_DEG
    )

    # 4. Populate the coordinate data
    src_lon_ptr = src_grid.get_coords(0)
    src_lat_ptr = src_grid.get_coords(1)
    src_lon_ptr[...] = src_lon
    src_lat_ptr[...] = src_lat

    dst_lon_ptr = dst_grid.get_coords(0)
    dst_lat_ptr = dst_grid.get_coords(1)
    dst_lon_ptr[...] = dst_lon.T
    dst_lat_ptr[...] = dst_lat.T

    # 5. Create esmpy Fields
    src_field = esmpy.Field(src_grid, name="src_field")
    dst_field = esmpy.Field(dst_grid, name="dst_field")
    src_field.data[...]=np.sqrt(x1/180)/(90+y1) # arbitrary function of x,y

    np.savetxt('F.txt', src_field.data[...])
    np.savetxt('X.txt', x1)
    np.savetxt('Y.txt', y1)
    print(f"Creating weights: {weights_file}")
    regrid = esmpy.Regrid(
      src_field,
      dst_field,
      filename=weights_file,
      regrid_method=esmpy.RegridMethod.BILINEAR,
      unmapped_action=esmpy.UnmappedAction.IGNORE # Optional: Ignores missing/masked points
    )
    np.savetxt('Fi.txt', dst_field.data[...])
    np.savetxt('xi.txt', xi)
    np.savetxt('yi.txt', yi)
    
    return

# Main program

nargin = len(sys.argv) - 1

if nargin < 5 :
    print("missing argument for specifying spatial forecast variance")
    print("call InterpWindToMesh.crvln.py flin mesh flout variance lambda")

flin=sys.argv[1]
mshfl=sys.argv[2]
flout=sys.argv[3]
ForecastErrorVariance=float(sys.argv[4]) #(m m /s /s)
LambdaBndyErrorVariance=float(sys.argv[5]) #(km)

xi, yi, ei = nwps.loadWW3Mesh(mshfl)

nn=len(xi)

if np.mean(xi)<0:
	xi=xi+360.

data = nc.Dataset(flin,"r")
#read spaital dimensions and determine if input mesh is curvilinear or regular
x1=np.asarray(data["longitude"][:])
y1=np.asarray(data["latitude"][:])
if (len(x1.shape)==2 and len(y1.shape)==2):
    IsCrvLn=True
elif (len(x1.shape)==1 and len(y1.shape)==1):
    IsCrvLn=False
else:
    print("input file spatial dimension is not right. ending program")
    sys.exit()
# shift to common longitude
if np.mean(np.mean(x1))<0:
    x1=x1+360.

#######################################
# === Create weights or resuse it  ===#
#######################################

meshslash=mshfl.rfind('/')+1
dom=flin.split(".")
dom=dom[len(dom)-2]
weights_file = "InterpolationWeights."+mshfl[meshslash:len(mshfl)-3]+dom+".nc"
print("interpolation weights file = "+ weights_file)

if IsCrvLn:
    nx=x1.shape[0]
    ny=x1.shape[1]
else:
    nx=len(x1)
    ny=len(y1)

if not os.path.isfile(weights_file):
    print("No existing interpolation weights. Computing weights and saving to file: "+ weights_file)
    if not IsCrvLn:
        x1 = np.tile(x1,(ny,1))
        y1 = np.tile(y1,(nx,1)).T
    CurvilinearGridCreateInterpWeights(xi, yi, x1, y1, weights_file)

n1=nx*ny

##################################################################################
# START: Interpolate wind data using predefined weights
##################################################################################

#######################################
# === resuse weights  ===#
#######################################
with xr.open_dataset(weights_file) as ds_s:
   # Standard sparse storage uses 'row', 'col', and 'data' variables
   row = ds_s['row'].values
   col = ds_s['col'].values
   weights = ds_s['S'].values
matrix = sp.coo_matrix((weights, (row-1, col-1)), shape=(nn,n1)).tocsr()
print("sparse interpolation matrix")
print(matrix)
######################################
t1=np.asarray(data["time"][:])
nt=len(t1)
nan=float("nan")

# read and reshape input data for interpolation
Up=np.zeros((nt,n1))+nan
Vp=np.zeros((nt,n1))+nan
fill_value0=data["UGRD_10maboveground"]._FillValue
print("fill value="+str(fill_value0))
for k in range(nt):
    U=data["UGRD_10maboveground"][k,:,:]
    Up[k,:]=np.transpose(U).reshape(n1)
    V=data["VGRD_10maboveground"][k,:,:]
    Vp[k,:]=np.transpose(V).reshape(n1)

j=np.where( Up==fill_value0 )
Up[j]=nan
Vp[j]=nan
j=np.where( Vp==fill_value0 )
Up[j]=nan
Vp[j]=nan

#Cary out interpolation with sparse matrix multiplication
u = matrix @ Up.T
v = matrix @ Vp.T

#Use fill value where the interpolator has no coverage
row_sum = matrix.sum(axis=1)
j0=np.where(row_sum==0)
u[j0,:]=nan
v[j0,:]=nan

##################################################################################
# FINISHED: Interpolate wind data using predefined weights
##################################################################################

##################################################################################
# START: Compute distance to boundary for each node in mesh:
##################################################################################
#dist2bnd_file = "DistToBndy."+mshfl[meshslash:len(mshfl)-3]+dom+".txt"
dist2bnd_file = "DistToBndy."+mshfl[meshslash:len(mshfl)-3]+dom+".nc"
if os.path.isfile(dist2bnd_file):
    print(f"Reusing existing weights: {weights_file}")
    #dist2bnd = np.loadtxt(dist2bnd_file)
    with xr.open_dataset(dist2bnd_file) as ds_s:
        dist2bnd = ds_s['dist2bnd'].values
else:
    nx=x1.shape[0]
    ny=x1.shape[1]
    xb=np.hstack((x1[1,:],x1[:,ny-1].T,x1[nx-1,:],x1[:,1].T))
    yb=np.hstack((y1[1,:],y1[:,ny-1].T,y1[nx-1,:],y1[:,1].T))
    dist2bnd=np.zeros(nn)
    for k in range(nn):
        dist2bnd[k]=QuickDistance(yi[k],xi[k],yb,xb)
        if k%10000==0:
            print("calculating distance to boundary, "+str(k)+":"+ str(nn)+":"+str(k/nn) )
    #np.savetxt(dist2bnd_file, dist2bnd, '%f')
    with nc.Dataset(dist2bnd_file, 'w', format='NETCDF4') as ncout:
        ncout.createDimension('node' , nn)
        d_var=ncout.createVariable('dist2bnd', 'f4', ('node',))
        d_var.long_name     = 'distance to boundary'
        d_var.units         = 'km'
        d_var.standard_name = 'distance to boundary'
        d_var[:]=dist2bnd[:]
##################################################################################
# FINISHED: Compute distance to boundary for each node in mesh:
##################################################################################


##################################################################################
# START: Assign forecast error variance based on distance to boundary and inputs
##################################################################################
if LambdaBndyErrorVariance > 0. :
    um0=u[:,0]
    InteriorNodeList=np.where(um0**2>=0) #find points that are valid floats
    Variance0=VarianceLinearDistanceToBndy(InteriorNodeList, dist2bnd, ForecastErrorVariance, 10.*ForecastErrorVariance,LambdaBndyErrorVariance )
else :
    #Assign spatially constant error variance if LambdaBndyErrorVariance<=0
    Variance0=0.*xi + ForecastErrorVariance 
Variance=np.zeros((nt,nn))
for k in range(nt):
    Variance[k,:]=Variance0
##################################################################################
# FINISHED: Assign forecast error variance based on distance to boundary and inputs
##################################################################################



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

    u_var=ncout.createVariable('uwnd', 'f4', ('time','node'),fill_value    = fill_value0)
    u_var.long_name     = 'eastward_wind'
    u_var.units         = 'm/s'
    u_var.standard_name = 'eastward_wind'
    u_var.level = '10 m above ground'
    u_var[:,:]=u[:,:].T

    v_var=ncout.createVariable('vwnd', 'f4', ('time','node'),fill_value    = fill_value0)
    v_var.long_name     = 'northward_wind'
    v_var.units         = 'm/s'
    v_var.standard_name = 'northward_wind'
    v_var.level = '10 m above ground'
    v_var[:,:]=v[:,:].T

    ErrorVariance_var=ncout.createVariable('ErrorVariance', 'f4', ('time','node'),fill_value    = fill_value0)
    ErrorVariance_var.long_name     = 'forecast error variance'
    ErrorVariance_var.units         = 'm m /s /s'
    ErrorVariance_var.standard_name = 'variance'
    ErrorVariance_var.level = '10 m above ground'
    ErrorVariance_var[:,:]=Variance[:,:]

    ncout.close

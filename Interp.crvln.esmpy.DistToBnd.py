


import numpy as np

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

#import xarray as xr
import numpy as np
#import xesmf as xe
import esmpy as ESMF

import scipy.sparse as sp

flin=sys.argv[1]
mshfl=sys.argv[2]
flout=sys.argv[3]

#mshfl="meshes/RWPS.V0a.msh"
xi, yi, ei = nwps.loadWW3Mesh(mshfl)

nn=len(xi)
print(nn)

if np.mean(xi)<0:
	xi=xi+360.

data = nc.Dataset(flin,"r")
t1=np.asarray(data["time"][:])

nn=len(xi)
nt=len(t1)

x1=np.asarray(data["longitude"][:])
y1=np.asarray(data["latitude"][:])
if np.mean(np.mean(x1))<0:
        x1=x1+360.
src_lon=data["longitude"]
src_lat=data["latitude"]

print("x1,y1 shape")
nx=x1.shape[0]
ny=x1.shape[1]
print(x1.shape)
print(y1.shape)
print(x1)

n1=nx*ny
src_lon = x1
src_lat = y1

dst_lon=np.zeros((1,nn))
dst_lat=np.zeros((1,nn))
dst_lon[0,:] = xi[:]
dst_lat[0,:] = yi[:]


# 3. Create ESMF Grid objects
# For curvilinear grids, specify StaggerLoc.CENTER and provide 2D arrays
src_grid = ESMF.Grid(
    max_index=np.array([src_lon.shape[0], src_lon.shape[1]]),
    staggerloc=ESMF.StaggerLoc.CENTER,
    coord_sys=ESMF.CoordSys.SPH_DEG
)

#####
dst_grid = ESMF.Grid(
    max_index=np.array([nn, 1]),
    staggerloc=ESMF.StaggerLoc.CENTER,
    coord_sys=ESMF.CoordSys.SPH_DEG
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

# 5. Create ESMF Fields
src_field = ESMF.Field(src_grid, name="src_field")
dst_field = ESMF.Field(dst_grid, name="dst_field")

src_field.data[...] = data["UGRD_10maboveground"][0,:,:]
np.savetxt('u0.txt', src_field.data[...])

regrid = ESMF.Regrid(
    src_field, 
    dst_field, 
    regrid_method=ESMF.RegridMethod.BILINEAR,
    unmapped_action=ESMF.UnmappedAction.IGNORE, # Optional: Ignores missing/masked points
    factors=True # Enables the extraction of weights
)

factors, factors_index = regrid.get_factors()
print("Interpolation Factors (Weights):", factors)
print("Source Indices (Columns), Destination Indices (Rows):", factors_index)
# 6. Verify and read the created NetCDF Weights file
print("Successfully interpolated data. Weights exported to esmpy_weights.nc")

row=factors_index[:,1]
col=factors_index[:,0]

nw=len(factors)
np.savetxt('factors.txt', factors)
np.savetxt('factors_index.txt', factors_index)

WeightFile='SparseMatrixInterpWeights.nc'

matrix = sp.coo_matrix((factors, (row-1, col-1)), shape=(nn,n1)).tocsr()

nnz=len(factors) # number of non-zero elements of sparse interpolation matrix

with nc.Dataset(WeightFile, 'w', format='NETCDF4') as ncout:
    ncout.createDimension('K' , nnz)
    n_var = ncout.createVariable('N','i4', ())
    n_var[:]=nn
    m_var = ncout.createVariable('M','i4', ())
    m_var[:]=n1
    row_index=ncout.createVariable('row', 'i4', ('K',))
    row_index[:]=row[:]
    col_index=ncout.createVariable('col', 'i4', ('K',))
    col_index[:]=col[:]
    val_var=ncout.createVariable('value', 'f8', ('K',))
    val_var[:]=factors[:]

nan=float("nan")

# define and reshape input data for interpolation
Up=np.zeros((nt,n1))+nan
Vp=np.zeros((nt,n1))+nan
fill_value0=data["UGRD_10maboveground"]._FillValue
print("fill value="+str(fill_value0))

for k in range(nt):
    U=data["UGRD_10maboveground"][k,:,:]
    Up[k,:]=np.transpose(U).reshape(n1)
    V=data["VGRD_10maboveground"][k,:,:]
    Vp[k,:]=np.transpose(V).reshape(n1)
row_sum = matrix.sum(axis=1)
u = matrix @ Up.T
v = matrix @ Vp.T

#Use fill value where the interpolator has no coverage
j0=np.where(row_sum==0)
u[j0,:]=nan
v[j0,:]=nan
#############################################################################
# Compute distance to boundary for each node in mesh:
# only needs to be done for nodes within interpolater
##################################################################################

xb=np.hstack((x1[1,:],x1[:,ny-1].T,x1[nx-1,:],x1[:,1].T))
yb=np.hstack((y1[1,:],y1[:,ny-1].T,y1[nx-1,:],y1[:,1].T))

np.savetxt('xbyb.txt', np.vstack((xb,yb)))
np.savetxt('xiyi.txt', np.vstack((xi,yi)))

dist2bnd=np.zeros(nn)
for k in range(nn):
    dist2bnd[k]=QuickDistance(yi[k],xi[k],yb,xb)
    if k%10000==0:
        print("calculating distance to boundary, "+str(k)+":"+ str(nn)+":"+str(k/nn) )
##################################################################################

#fill_value0=-99999
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
    u_var[:,:]=u[:,:]

    v_var=ncout.createVariable('vwnd', 'f4', ('node','time'),fill_value    = fill_value0)
    v_var.long_name     = 'northward_wind'
    v_var.units         = 'm/s'
    v_var.standard_name = 'northward_wind'
    v_var.level = '10 m above ground'
    v_var[:,:]=v[:,:]

    ncout.close

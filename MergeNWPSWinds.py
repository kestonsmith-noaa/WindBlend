

import datetime
import numpy as np
import netCDF4 as nc
import sys
import re
import InterpNWPSUtility as nwps

nwpsdom=nwps.GetDomainStrings(31)
print(nwpsdom)
mshfl="RWPS.PIXXL.PI1km.PP.WW3bEle.msh"
xi, yi, ei = nwps.loadWW3Mesh(mshfl)
nn=len(xi)
nt=171 #this needs to be set and time added to field

u=np.zeros((nn,nt))
v=np.zeros((nn,nt))

Ndomains=33
"""
for k in range(Ndomains):
    nwpsdom=nwps.GetDomainStrings(k)
    flu=nwpsdom+".ui.txt"
    flv=nwpsdom+".vi.txt"
    flj=nwpsdom+".ji.txt"
    flt=nwpsdom+".time.txt"
    u0 = np.loadtxt(flu)
    v0 = np.loadtxt(flv)
    j0 = np.loadtxt(flj, dtype=int)
    t0 = np.loadtxt(flt)
    u[j0,:]=u0
    v[j0,:]=v0
"""
flbg="../NBM/NBM20260108/NBM.wind.01.08.WW3.nc"
data = nc.Dataset(flbg,"r")
x1=np.asarray(data["longitude"][:])
x1=x1-360 # shift grid to RWPS coordinates 
y1=np.asarray(data["latitude"][:])
t1=np.asarray(data["time"][:])
u1=np.asarray(data["UGRD_10maboveground"][:,:,:])
v1=np.asarray(data["VGRD_10maboveground"][:,:,:])

from scipy.interpolate import RegularGridInterpolator

InterpNBMu = RegularGridInterpolator((x1, y1), np.transpose(u1), method='linear', bounds_error=False, fill_value=None)
InterpNBMv = RegularGridInterpolator((x1, y1), np.transpose(v1), method='linear', bounds_error=False, fill_value=None)
uit=InterpNBMu((xi,yi))
vit=InterpNBMv((xi,yi))
#np.savetxt("NBM.ui.txt",uit,"%f")
#np.savetxt("NBM.vi.txt",vit,"%f")

k=31
nwpsdom=nwps.GetDomainStrings(k)
flu=nwpsdom+".ui.txt"
flv=nwpsdom+".vi.txt"
flj=nwpsdom+".ji.txt"
flt=nwpsdom+".time.txt"
u0 = np.loadtxt(flu)
v0 = np.loadtxt(flv)
j0 = np.loadtxt(flj, dtype=int)
t0 = np.loadtxt(flt)

nt=len(t1)
uit[j0,:]=u0[:,0:nt]
vit[j0,:]=v0[:,0:nt]
#np.savetxt("NBMsew.ui.txt",uit,"%f")
#np.savetxt("NBMsew.vi.txt",vit,"%f")

j=np.where(np.isnan(uit+vit))
fill_value0=-999999
uit[j]=fill_value0
vit[j]=fill_value0

nn=len(xi)
nt=len(t1)
flout="NBMmixNWPS.wind.WW3.nc"
MAPSTA=np.ones(nn, dtype=int)
print(ei.shape)
print(ei)
ne=ei.shape[0]

time=np.zeros(nt)
for k in range(nt):
    time[k]= nwps.UnixTimeToDaysSince1990(t1[k])

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
    time_var.units         = 'days since 1990-01-01 00:00:00'
    time_var.long_name     = 'julian day (UT)'
    time_var.standard_name = 'time'
    time_var.axis          = 'T'
    time_var[:]=time[:]

    tri_var=ncout.createVariable('tri', 'i4', ('noel','element'))
    tri_var.long_name     = 'element list'
    tri_var.standard_name = 'element list'
    tri_var[:]=np.transpose(ei)
    
    map_var=ncout.createVariable('MAPSTA', 'i2', ('node',))
    map_var.units         = '1'
    map_var.long_name     = 'status map'
    map_var.standard_name = 'status map'
    map_var.axis          = 'node'
    map_var[:]=MAPSTA[:]

    u_var=ncout.createVariable('uwnd', 'f4', ('node','time'),fill_value    = fill_value0)
    u_var.long_name     = 'eastward_wind'
    u_var.units         = 'm/s'
    u_var.standard_name = 'eastward_wind'
    u_var.level = '10 m above ground'
    u_var[:,:]=uit[:,:]

    v_var=ncout.createVariable('vwnd', 'f4', ('node','time'),fill_value    = fill_value0)
    v_var.long_name     = 'northward_wind'
    v_var.units         = 'm/s'
    v_var.standard_name = 'northward_wind'
    v_var.level = '10 m above ground'
    v_var[:,:]=vit[:,:]

    ncout.close

"""

fl='/mnt/sda/keston/WW3/WW3/regtests/ww3_tp2.17/input/wind.nc'
ncdisp(fl)


Source:
           /mnt/sda/keston/WW3/WW3/regtests/ww3_tp2.17/input/wind.nc
Format:
           netcdf4
Global Attributes:
           WAVEWATCH_III_version_number    = '6.03'
           WAVEWATCH_III_switches          = 'F90 NOGRB TRKNC NOPA LRB4 DIST MPI SCRIP MLIM PR3 UQ NC4 FLX0 LN1 ST4 STAB0 NL1 BT4 DB1 TR0 BS0 IS0 IC0 REF0 XX0 WCUR WNT2 WNX1 RWND CRT1 CRX1 O0 O1 O2 O2a O2b O2c O3 O4 O5 O6 O7'
           SIN4 namelist parameter BETAMAX = 1.33
           product_name                    = 'ww3.201512_wnd.nc'
           area                            = 'Inlet'
           latitude_resolution             = 'n/a'
           longitude_resolution            = 'n/a'
           southernmost_latitude           = '40.38446'
           northernmost_latitude           = '40.99023'
           westernmost_longitude           = '-72.92410'
           easternmost_longitude           = '-72.03251'
           minimum_altitude                = '-12000 m'
           maximum_altitude                = '9000 m'
           altitude_resolution             = 'n/a'
           start_date                      = '2015-12-14 00:00:00'
           stop_date                       = '2015-12-15 00:00:00'
Dimensions:
           level   = 1
           node    = 3070
           element = 5780
           time    = 25    (UNLIMITED)
           noel    = 3
Variables:
    longitude
           Size:       3070x1
           Dimensions: node
           Datatype:   single
           Attributes:
                       units         = 'degree_east'
                       long_name     = 'longitude'
                       standard_name = 'longitude'
                       valid_min     = -180
                       valid_max     = 360
                       axis          = 'X'
    latitude 
           Size:       3070x1
           Dimensions: node
           Datatype:   single
           Attributes:
                       units         = 'degree_north'
                       long_name     = 'latitude'
                       standard_name = 'latitude'
                       valid_min     = -90
                       valid_max     = 180
                       axis          = 'Y'
    time     
           Size:       25x1
           Dimensions: time
           Datatype:   double
           Attributes:
                       long_name     = 'julian day (UT)'
                       standard_name = 'time'
                       units         = 'days since 1990-01-01 00:00:00'
                       conventions   = 'relative julian days with decimal part (as parts of the day )'
                       axis          = 'T'
    tri      
           Size:       3x5780
           Dimensions: noel,element
           Datatype:   int32
    MAPSTA   
           Size:       3070x1
           Dimensions: node
           Datatype:   int16
           Attributes:
                       long_name     = 'status map'
                       standard_name = 'status map'
                       units         = '1'
                       valid_min     = -32
                       valid_max     = 32
    uwnd     
           Size:       3070x25
           Dimensions: node,time
           Datatype:   single
           Attributes:
                       long_name     = 'eastward_wind'
                       standard_name = 'eastward_wind'
                       globwave_name = 'eastward_wind'
                       units         = 'm s-1'
                       _FillValue    = 9.969209968386869e+36
                       scale_factor  = 1
                       add_offset    = 0
                       valid_min     = -990
                       valid_max     = 990
                       comment       = 'wind=sqrt(U10**2+V10**2)'
    vwnd     
           Size:       3070x25
           Dimensions: node,time
           Datatype:   single
           Attributes:
                       long_name     = 'northward_wind'
                       standard_name = 'northward_wind'
                       globwave_name = 'northward_wind'
                       units         = 'm s-1'
                       _FillValue    = 9.969209968386869e+36
                       scale_factor  = 1
                       add_offset    = 0
                       valid_min     = -990
                       valid_max     = 990
                       comment       = 'wind=sqrt(U10**2+V10**2)'


#np.savetxt("NBM.ui.txt",np.transpose(uit),"%f")
#np.savetxt("NBM.vi.txt",np.transpose(vit),"%f")
"""

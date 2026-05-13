

import datetime
import numpy as np
import netCDF4 as nc
import sys
import re
#import InterpNWPSUtility as nwps

datestr=sys.argv[1]
cycl=sys.argv[2]

winddir="forecasts/wind."+datestr+"."+cycl
outdir="rwps_winds."+datestr+"."+cycl


rwps_pr=outdir+"/rrfs."+datestr+"."+cycl+".wind10m.pr.nc"
rwps_hi=outdir+"/rrfs."+datestr+"."+cycl+".wind10m.hi.nc"
rwps_na=outdir+"/rrfs."+datestr+"."+cycl+".wind10m.na.nc"
rwps_ak=outdir+"/rrfs."+datestr+"."+cycl+".wind10m.ak.nc"
rwps_conus=outdir+"/rrfs."+datestr+"."+cycl+".wind10m.conus.nc"

#file that contains actual rwps times
rwps_oc_ti=outdir+"/nbm."+datestr+"."+cycl+".wind10m.oc.ti.nc"

rwps_wind_out=outdir+"/rwps.windblend."+datestr+"."+cycl+".wind10m.nc"

BackgroundVariance=100. #variance for nbm

LocalFS  = [ rwps_pr, rwps_hi, rwps_ak, rwps_conus, rwps_na] #file names
VarFS    = [ 4.     , 4.    , 9.      , 16.       , 25.    ] # m m /s /s
LambdaFS = [ 150.   , 200.  , 500.    , 1000.     , 1500.  ] #km


#LocalFS  = [ rwps_pr, rwps_hi, rwps_ak ] #file names
#VarFS    = [ 1.      , 1.    , 4.      ] # m m /s /s
#LambdaFS = [ 150.     , 250.  , 1000.    ] #km

print("updating background forecast "+rwps_oc_ti+ "with forecasts:")
print(LocalFS)


data = nc.Dataset(rwps_oc_ti,"r")
x=np.asarray(data["longitude"][:])
#x1=x1-360 # shift grid to RWPS coordinates 
y=np.asarray(data["latitude"][:])
t=np.asarray(data["time"][:])
u=np.asarray(data["uwnd"][:,:])
v=np.asarray(data["vwnd"][:,:])
e=np.asarray(data["tri"][:,:])

nt=len(t)
nn=len(x)
ne=e.shape[1]

var0=BackgroundVariance+np.zeros((nn,nt)) # prior variance for BG field 

nFS=len(VarFS)
for jFS in range(nFS):
    flm=LocalFS[jFS]
    print("updating forecast winds based on: "+flm)
    datam = nc.Dataset(flm,"r")
    um=np.asarray(datam["uwnd"][:,:])
    vm=np.asarray(datam["vwnd"][:,:])
    dm=np.asarray(datam["dist2bnd"][:])
    tm=np.asarray(datam["time"][:])

    #Set fill values to nan
    nan=float('nan')
    fill_valueu = datam["uwnd"]._FillValue
    j=np.where(um==fill_valueu)
    um[j]=nan
    fill_valuev = datam["vwnd"]._FillValue
    j=np.where(vm==fill_valuev)
    vm[j]=nan

    #Find which spatial points have valid forecasts
    #these are assumed to not change in time 
    um0=um[:,0]
    ng0=np.where(um0**2>=0) #find points that are valid floats

    varm=np.zeros(len(dm))+np.inf
    lambdam=LambdaFS[jFS]
    varm0=VarFS[jFS]
    ForecastVariance=VarFS[jFS]

    #Set spatial variance structure based on distance
    #to local domain boundary (in domain specified ny ng0)
#    varm[ng0]=varm0*lambdam/dm[ng0]
#    j=np.where(varm<varm0)
#    varm[j]=varm0
    SpatialFunction=dm[ng0]/lambdam
    j=np.where(SpatialFunction>1.)
    SpatialFunction[j]=1.
    EdgeVar=10.*BackgroundVariance
    varm[ng0]=EdgeVar + (ForecastVariance-EdgeVar)*SpatialFunction

    np.savetxt("um0."+str(jFS)+".txt",um0,"%f")
    np.savetxt("ng0."+str(jFS)+".txt",ng0,"%i")
    np.savetxt("varm."+str(jFS)+".txt",varm,"%f")

    for k in range(nt):
        j=np.where(t[k]==tm) # find common merge point in data
        j=j[0].tolist()
        if len(j)==1:
            print("j="+str(j)+" : k="+str(k))
            u[ng0,k]=u[ng0,k]+(var0[ng0,k]/(var0[ng0,k]+varm[ng0]))*(um[ng0,j]-u[ng0,k])
            v[ng0,k]=v[ng0,k]+(var0[ng0,k]/(var0[ng0,k]+varm[ng0]))*(vm[ng0,j]-v[ng0,k])
    #        var0[ng0,k]=var0[ng0,k] * ( 1. - ( var0[ng0,k] / (var0[ng0,k]+varm[ng0] ) ) )
            var0[ng0,k]=var0[ng0,k] * ( varm[ng0] / ( var0[ng0,k]+varm[ng0] ) )


  Grid name : Inlet                         

  Comment character is '$'


#  Description of inputs
# --------------------------------------------------
#       Input type        : winds         
#       Format type       : pre-processed file  

#           File name         : ../../forcing/rwps.windblend.20260502.00.wind10m.nc                             
#           Dimension along x : time
#           Dimension along y : 
#           Field component 1 : uwnd
#           Field component 2 : vwnd

# *** WAVEWATCH III WARNING IN W3PRNC : 
#     calendar ATTRIBUTE NOT DEFINED
#     DEFAULTING TO "standard" CALENDAR
#     INPUT FILE MUST RESPECT STANDARD/GREGORIAN CALENDAR

# *** WAVEWATCH III ERROR IN W3PRNC : 
#     _FillValue ATTRIBUTE NOT DEFINED FOR : uwnd
fill_value0=-999999

with nc.Dataset(rwps_wind_out, 'w', format='NETCDF4') as ncout:
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
    lon_var[:]=x[:]

    lat_var=ncout.createVariable('latitude', 'f8', ('node',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'Y'
    lat_var[:]=y[:]

#    d_var=ncout.createVariable('priorvariance', 'f4', ('node',))
#    d_var.long_name     = 'variance'
#    d_var.units         = 'm m / s / s'
#    d_var.standard_name = 'variance'
#    d_var[:]=varm[:]

    time_var=ncout.createVariable('time', 'f8', ('time',))
    time_var.units         = 'seconds since 1970-01-01 00:00:00.0 0:00'
    time_var.long_name     = 'verification time generated by wgrib2 function verftime()'
    time_var.standard_name = 'time'
    time_var.axis          = 'T'
    time_var.reference_time = 1777334400
    time_var.reference_date = '2026.04.28 00:00:00 UTC'
    time_var[:]=t[:]

    tri_var=ncout.createVariable('tri', 'i4', ('noel','element'))
    tri_var.long_name     = 'element list'
    tri_var.standard_name = 'element list'
    tri_var[:]=e

    u_var=ncout.createVariable('uwnd', 'f4', ('node','time'),fill_value    = fill_value0)
#    u_var=ncout.createVariable('uwnd', 'f4', ('node','time'))
    u_var.long_name     = 'eastward_wind'
    u_var.units         = 'm/s'
    u_var.standard_name = 'eastward_wind'
    u_var.level = '10 m above ground'
    u_var[:,:]=u[:,:]

    v_var=ncout.createVariable('vwnd', 'f4', ('node','time'),fill_value    = fill_value0)
#    v_var=ncout.createVariable('vwnd', 'f4', ('node','time'))
    v_var.long_name     = 'northward_wind'
    v_var.units         = 'm/s'
    v_var.standard_name = 'northward_wind'
    v_var.level = '10 m above ground'
    v_var[:,:]=v[:,:]

    vari_var=ncout.createVariable('postvariance', 'f4', ('node','time'))
    vari_var.long_name     = 'variance of wind estimate'
    vari_var.units         = 'm m / s / s'
    vari_var.standard_name = 'variance of wind estimate'
    vari_var.level = '10 m above ground'
    vari_var[:,:]=var0[:,:]

    ncout.close

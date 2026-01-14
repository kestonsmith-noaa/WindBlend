

from scipy.interpolate import RegularGridInterpolator
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

import datetime
import numpy as np
import netCDF4 as nc
import sys
import re
import InterpNWPSUtility as nwps


knts2mps=1./1.94384001
fill_value=-999999


flin="sew_202512010102_WIND.nc"
#flin="ajk_202501010113_WIND.nc"

mshfl="RWPS.PIXXL.PI1km.PP.WW3bEle.msh"
xi,yi = nwps.loadWW3MeshCoords(mshfl)

nwpsdom="sew"
u,v,x,y,t=nwps.ReadNWPSWind(flin)

xu=np.max(x)
yu=np.max(y)
xl=np.min(x)
yl=np.min(y)

print("x min="+str(xl) + ", x max="+str(xu))
print("xi min="+str(np.min(xi)))
print("xi max="+str(np.max(xi)))
j=np.where( np.logical_and(   np.logical_and(xi>xl,xi<xu), np.logical_and(yi>yl,yi<yu) ) )
j=np.array(j).flatten()
xit=xi[j]
yit=yi[j]
print(j.shape)
print(u.shape)

print("start u interp")
uip=nwps.interpolate_curvilinear_to_pointsMD(x,y, np.transpose(u), xit, yit)
print("done u interp")
print("start v interp")
vip=nwps.interpolate_curvilinear_to_pointsMD(x,y, np.transpose(v), xit, yit)
print("done v interp")

j0=np.where(~np.isnan(uip[:,0]))
j=j[j0]
uip=uip[j0,:] # remove NaN rows
vip=vip[j0,:] 
xit=xi[j]
yit=yi[j]

np.savetxt(nwpsdom+".ui.txt",uip,"%f")  
np.savetxt(nwpsdom+".vi.txt",vip,"%f")  
np.savetxt(nwpsdom+".xi.txt",xit,"%f")  #Not nesesary 
np.savetxt(nwpsdom+".yi.txt",yit,"%f")  #Not nesesary
np.savetxt(nwpsdom+".time.txt",t,"%f")  #Not nesesary
np.savetxt(nwpsdom+".ji.txt",j,"%i")  


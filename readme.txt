

This is a set of routines for retreiving wind forecasts from various source and interpolating them to the nodes of an unstructured mesh (for wavewatch III). It is assumed that there is a "background" forecast which spans the entire domain.  Each forecast is assigned a spatially variable error variance based on proximity to the forecast domain's boundary. A common time frame for the forecast is established, whose range is the background forecast and also includes times from the other forecasts that fall within the range.

In current form these routines use the National Blend of Models (NBM) oceanic (OC) domain as the background forecast. Subsequent forecasts are taken from the Rapid Refresh Forecast System (RRFS), Alaska (AK), Hawaii (HI), Puerto Rico(PR), Continental US (CONUS) and North America (NA) domains.  the  Tools here take are The structure of the different forecasts is A common time frame for the 

To run for May 31 2026 forecast cycle 00 run,

$sh GetWinds.sh 20260531 00
$sh ProcessFcasts.sh 20260531 00 meshes/RWPS.V0a.small.msh

The last argument of ProcessFcasts.sh is a WW3 unstructured mesh
this will create a directory, rwps_winds.RWPS.V0a.small.20260602.00 containing files:

nbm.RWPS.V0a.small.20260602.00.wind10m.oc.ti.nc <- NBM wind forecast interpolated to 
rwps.est.RWPS.V0a.small.20260602.00.wind10m.nc  <- NBM updated with the five RRFS forecasts.

There are also several intermidiate files.

For computational efficiency interpolation weights are saved to files (and reused if available).  Distance to boundary is handled in the same manner,

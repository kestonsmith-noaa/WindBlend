#!/bin/bash



#module purge
#module load python/3.11
#module load rdhpcs-conda/25.11.0
conda activate xesmf_envb

pip list -v

date="20260507"
cycl="00"
winddir="forecasts/wind.$date.$cycl"
outdir="rwps_winds.$date.$cycl"

nbm_oc="$winddir/nbm.$date.$cycl.wind10m.oc.nc"
nbm_oc_uv="$winddir/nbm.$date.$cycl.wind10m.oc.uv.nc"

rrfs_pr="$winddir/rrfs.$date.$cycl.wind10m.pr.nc"
rrfs_hi="$winddir/rrfs.$date.$cycl.wind10m.hi.nc"
rrfs_na="$winddir/rrfs.$date.$cycl.wind10m.na.nc"
rrfs_ak="$winddir/rrfs.$date.$cycl.wind10m.ak.nc"
rrfs_conus="$winddir/rrfs.$date.$cycl.wind10m.conus.nc"

rwps_oc="$outdir/nbm.$date.$cycl.wind10m.oc.nc"
rwps_oc_ti="$outdir/nbm.$date.$cycl.wind10m.oc.ti.nc"
rwps_pr="$outdir/rrfs.$date.$cycl.wind10m.pr.nc"
rwps_hi="$outdir/rrfs.$date.$cycl.wind10m.hi.nc"
rwps_na="$outdir/rrfs.$date.$cycl.wind10m.na.nc"
rwps_ak="$outdir/rrfs.$date.$cycl.wind10m.ak.nc"
rwps_conus="$outdir/rrfs.$date.$cycl.wind10m.conus.nc"

mkdir $outdir

# convert NBM spd,dir to u,v
rm $rwps_oc_uv
echo "python3 SpdDir2UVnbm.py $nbm_oc $nbm_oc_uv"
python3 SpdDir2UVnbm.py $nbm_oc $nbm_oc_uv

#interpolate(spatial) nbm wind forecasts to RWPS nodes
rm $rwps_oc
echo "python3 Interp.reg.DistToBnd.py $nbm_oc_uv $rwps_oc"
python3 Interp.reg.DistToBnd.nbm.py $nbm_oc_uv $rwps_oc

#interpolate(spatial) rrfs wind forecasts to RWPS nodes
rm $rwps_pr
echo "python3 Interp.reg.DistToBnd.py $rrfs_pr $rwps_pr"
python3 Interp.reg.DistToBnd.py $rrfs_pr $rwps_pr

#interpolate(temporal) nbm forecast wind to rrfs pr forecast times
rm $rwps_oc_ti
echo "python3 InterpTimeNBM.py $rwps_oc $rwps_oc_ti"
python3 InterpTimeNBM.py $rwps_oc $rwps_pr $rwps_oc_ti

#interpolate(spatial) rrfs wind forecasts to RWPS nodes
rm $rwps_hi
echo "python3 Interp.reg.DistToBnd.py $rrfs_hi $rwps_hi"
python3 Interp.reg.DistToBnd.py $rrfs_hi $rwps_hi

#interpolate(spatial) rrfs wind forecasts to RWPS nodes
rm $rwps_ak
echo "python3 Interp.crvln.esmf.DistToBnd.py $rrfs_ak $rwps_ak"
python3 Interp.crvln.esmf.DistToBnd.py $rrfs_ak $rwps_ak

#interpolate(spatial) rrfs wind forecasts to RWPS nodes
rm $rwps_na
echo "python3 Interp.crvln.esmf.DistToBnd.py $rrfs_na $rwps_na"
python3 Interp.crvln.esmf.DistToBnd.py $rrfs_na $rwps_na

#interpolate(spatial) rrfs wind forecasts to RWPS nodes
rm $rwps_conus
echo "python3 Interp.crvln.esmf.DistToBnd.py $rrfs_conus $rwps_conus"
python3 Interp.crvln.esmf.DistToBnd.py $rrfs_conus $rwps_conus


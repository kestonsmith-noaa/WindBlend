#!/bin/bash
#PBS -N ESMPy
#PBS -j oe
#PBS -S /bin/bash
#PBS -q dev
#PBS -A NWPS-DEV
#PBS -l walltime=00:05:00
#PBS -l select=1:ncpus=1:mem=8G
#PBS -l place=excl
#PBS -l debug=true

module reset
module load PrgEnv-intel/8.5.0
module load intel/19.1.3.304
module load craype/2.7.17
module load cray-mpich/8.1.19
module load hdf5-C/1.14.0
module load netcdf-C/4.9.2
module load esmf-C/8.6.0
module load ve/hafs/2.1

pip list -v

##date="20260527"
##cycl="00"
date=$1
cycl=$2
winddir="forecasts/wind.$date.$cycl"
#mesh="meshes/RWPS.V0a.msh"
mesh=$3

# extract mesh name from path
meshname="${mesh##*/}"
# remove .msh suffix from mesh name
meshname="${meshname: 0: -4}"

# incorporate meshname date and cycle into output directory name to avoid
# applying winds to wrong mesh
outdir="rwps_winds.$meshname.$date.$cycl"

echo "outputing files to: $outdir"

nbm_oc="$winddir/nbm.$date.$cycl.wind10m.oc.nc"
nbm_oc_uv="$winddir/nbm.$date.$cycl.wind10m.oc.uv.nc"

rrfs_pr="$winddir/rrfs.$date.$cycl.wind10m.pr.nc"
rrfs_hi="$winddir/rrfs.$date.$cycl.wind10m.hi.nc"
rrfs_na="$winddir/rrfs.$date.$cycl.wind10m.na.nc"
rrfs_ak="$winddir/rrfs.$date.$cycl.wind10m.ak.nc"
rrfs_conus="$winddir/rrfs.$date.$cycl.wind10m.conus.nc"

rwps_oc="$outdir/nbm.$meshname.$date.$cycl.wind10m.oc.nc"
rwps_oc_ti="$outdir/nbm.$meshname.$date.$cycl.wind10m.oc.ti.nc"
rwps_pr="$outdir/rrfs.$meshname.$date.$cycl.wind10m.pr.nc"
rwps_hi="$outdir/rrfs.$meshname.$date.$cycl.wind10m.hi.nc"
rwps_na="$outdir/rrfs.$meshname.$date.$cycl.wind10m.na.nc"
rwps_ak="$outdir/rrfs.$meshname.$date.$cycl.wind10m.ak.nc"
rwps_conus="$outdir/rrfs.$meshname.$date.$cycl.wind10m.conus.nc"

mkdir $outdir

# convert NBM spd,dir to u,v
rm $rwps_oc_uv
echo "python3 SpdDir2UVnbm.py $nbm_oc $nbm_oc_uv"
python3 SpdDir2UVnbm.py $nbm_oc $nbm_oc_uv

#interpolate(spatial) nbm wind forecasts to RWPS nodes
rm $rwps_oc
echo "python3 Interp.reg.DistToBnd.py $nbm_oc_uv $rwps_oc"
python3 Interp.reg.DistToBnd.nbm.py $nbm_oc_uv $mesh $rwps_oc

#interpolate(spatial) rrfs wind forecasts to RWPS nodes
rm $rwps_pr
echo "python3 Interp.reg.DistToBnd.py $rrfs_pr $rwps_pr"
python3 Interp.reg.DistToBnd.py $rrfs_pr $mesh $rwps_pr

#interpolate(temporal) nbm forecast wind to rrfs pr forecast times
rm $rwps_oc_ti
echo "python3 InterpTimeNBM.py $rwps_oc $rwps_oc_ti"
python3 InterpTimeNBM.py $rwps_oc $rwps_pr $rwps_oc_ti

#interpolate(spatial) rrfs wind forecasts to RWPS nodes
rm $rwps_hi
echo "python3 Interp.reg.DistToBnd.py $rrfs_hi $rwps_hi"
python3 Interp.reg.DistToBnd.py $rrfs_hi $mesh $rwps_hi

#interpolate(spatial) rrfs wind forecasts to RWPS nodes
rm $rwps_ak
echo "python3 Interp.crvln.esmf.DistToBnd.py $rrfs_ak $rwps_ak"
python3 Interp.crvln.esmpy.DistToBnd.py $rrfs_ak $mesh $rwps_ak

#interpolate(spatial) rrfs wind forecasts to RWPS nodes
rm $rwps_na
echo "python3 Interp.crvln.esmf.DistToBnd.py $rrfs_na $rwps_na"
python3 Interp.crvln.esmpy.DistToBnd.py $rrfs_na $mesh $rwps_na

#interpolate(spatial) rrfs wind forecasts to RWPS nodes
rm $rwps_conus
echo "python3 Interp.crvln.esmf.DistToBnd.py $rrfs_conus $rwps_conus"
python3 Interp.crvln.esmpy.DistToBnd.py $rrfs_conus $mesh $rwps_conus

#Blend rrfs winds with nbm
python3 BlendNBMwRRFS.LinVar.py $date $cycl $meshname $outdir

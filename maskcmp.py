import xarray
import sys
import numpy as np
import numpy.ma as ma

debug = False

# Grab the regridded data
ds1 = xarray.open_dataset(sys.argv[1])
elevation = ds1.elevation.to_masked_array()

print elevation.shape

# Set unmasked areas to zero, masked areas to 1
elevation[~elevation.mask] = 0
elevation[elevation.mask] = 1

if (debug):
    # Output to new file
    new = ds1.elevation
    new[:] = elevation[:]
    dsnew = new.to_dataset()
    dsnew.to_netcdf('tmp1.nc')

# Grab the topog data
ds2 = xarray.open_dataset(sys.argv[2])
depth = ds2.depth.to_masked_array()
print depth.shape

# Mask out land, set all ocean areas to zero, fill in masked land with 1
depth[depth == 0] = ma.masked
depth[depth > 0] = 0
depth[depth.mask] = 1

if (debug):
    # Output to new file
    new = ds2.depth
    new[:] = depth[:]
    dsnew = new.to_dataset()
    dsnew.to_netcdf('tmp2.nc')

diff = elevation - depth

new = ds1.elevation
new[:] = diff[:]
dsnew = new.to_dataset()
dsnew.to_netcdf('diff.nc')


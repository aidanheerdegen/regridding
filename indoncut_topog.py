import xray

ds = xray.open_dataset('sis01_mosaic/topog.nc')
dsnew = ds.isel(ny=slice(1083,1296),nx=slice(99,650))
dsnew.to_netcdf('topog_1080_90E-145E_15S_5N.nc')

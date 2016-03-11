import xray

ds = xray.open_dataset('sis01_mosaic/ocean_hgrid.nc')
dsnew = ds.isel(ny=slice(2168,2592),nx=slice(200,1300),nyp=slice(2167,2592),nxp=slice(199,1300))
dsnew.to_netcdf('ocean_hgrid_1080_90E-145E_15S_5N.nc')

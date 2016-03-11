program treegrid

  ! Regrid netcdf data file using kdtree

  use ncio, only: ncvar, nc_read, nc_open, nc_create, nc_write_dim, nc_write, nc_get_att, &
       nc_size, nc_print_attr, nc_v_init
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use iso_varying_string
  use string_functions, only: join
  use file_functions, only: exists, freeunit, stderr, stdout, open
  use precision
  use regrid_functions, only: regrid
  use pathfind_functions, only: earth_R

  implicit none

  character(len=2000) :: datafile, gridfile, fname, outfile
  character(len=200) :: latname, lonname, varname

  real, allocatable, target :: data(:,:)

  real, allocatable :: newdata(:,:)

  real, allocatable :: newgridx(:,:), newgridy(:,:)
  real, allocatable :: src_grid(:,:,:), dst_grid(:,:,:)
  real, allocatable :: lon(:), lat(:)

  integer :: error, id_data, id_grid, nlat, nlon, nlatnew, nlonnew, nx, ny

  type(ncvar) :: v, xgridv, ygridv

  logical :: supergrid

  real :: radius, oceanfrac, depth
  
  type (varying_string) :: myoptions(8)

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  myoptions(1) = 'help'
  myoptions(2) = 'var'
  myoptions(3) = 'lon'
  myoptions(4) = 'lat'
  myoptions(5) = 'radius'
  myoptions(6) = 'oceanfrac'
  myoptions(7) = 'depth'
  myoptions(8) = 'outfile'

  ! This call parses the command line arguments for command line options
  call get_options(myoptions, error)

  ! Check we weren't passed duff options -- spit the dummy if we were
  if (error > 0) then
     write(stderr,*) 'ERROR! Unknown options: ',join(bad_options()," ")
     call usage
     STOP
  end if

  ! Check if we just want to print the usage
  if (option_exists('help')) then
     call usage
     STOP
  end if

  if (num_args() < 2) then
     write(stderr,*) 'ERROR! Must supply input data file and new grid file as command-line arguments'
     call usage
     STOP
  end if

  if (option_exists('var')) then
     ! We have specified variable name
     if (.NOT. has_value('var')) then
        write(stderr,*) 'Option var must specify a value!'
        call usage
        stop
     end if
     varname = ""
     varname = get_value('var')
  else
     varname = "elevation"
  end if

  if (option_exists('lon')) then
     ! We have specified variable name for the longitude in dest grid
     if (.NOT. has_value('lon')) then
        write(stderr,*) 'Option lon must specify a value!'
        call usage
        stop
     end if
     lonname = ""
     lonname = get_value('lon')
  else
     lonname = "x"
  end if

  if (option_exists('lat')) then
     ! We have specified variable name for the longitude in dest grid
     if (.NOT. has_value('lat')) then
        write(stderr,*) 'Option lat must specify a value!'
        call usage
        stop
     end if
     latname = ""
     latname = get_value('lat')
  else
     latname = "y"
  end if

  if (option_exists('radius')) then
     ! Specified a radius (in km) over which to regrid for each pixel
     if (.NOT. has_value('radius')) then
        write(stderr,*) 'Option radius must specify a value!'
        call usage
        stop
     end if
     radius = get_value('radius')
  else
     radius = 100.
  end if
  ! Normalise by radius of earth. Input radius is defined in kilometres, earth_R is
  ! defined in metres
  radius = 1000. * (radius / earth_R)

  if (option_exists('oceanfrac')) then
     ! Specified a depth
     if (.NOT. has_value('oceanfrac')) then
        write(stderr,*) 'Option oceanfrac must specify a value!'
        call usage
        stop
     end if
     oceanfrac = get_value('oceanfrac')
  else
     oceanfrac = 0.5
  end if

  if (option_exists('depth')) then
     ! Specified a depth
     if (.NOT. has_value('depth')) then
        write(stderr,*) 'Option depth must specify a value!'
        call usage
        stop
     end if
     depth = get_value('depth')
  else
     depth = 0.
  end if

  if (option_exists('outfile')) then
     ! We have specified outfile name
     if (.NOT. has_value('outfile')) then
        write(stderr,*) 'Option outfile must specify a value!'
        call usage
        stop
     end if
     outfile = ""
     outfile = get_value('outfile')
  else
     outfile = "regridded.nc"
  end if
  ! Read in data to be re-gridded
  datafile = next_arg()

  print *,'Reading in data file: ',trim(datafile)

  call nc_open(trim(datafile), id_data, writable=.false.)

  ! Initialize the netcdf variable info and load attributes
  call nc_v_init(v,trim(varname))
  call nc_get_att(id_data,v,readmeta=.TRUE.)
  ! call nc_print_attr(v)

  nlon = v%dlen(1)
  nlat = v%dlen(2)

  allocate(lon(nlon),lat(nlat),src_grid(2,nlon,nlat),data(nlon,nlat))

  print *,nlon,'x',nlat

  ! Read longitude and latitude variables
  call nc_read(trim(fname),trim(v%dims(1)),lon,ncid=id_data)
  call nc_read(trim(fname),trim(v%dims(2)),lat,ncid=id_data)

  ! Need a lat/lon value for each point.
  ! TODO: Should test for dimensions of lat and lon and support regridding from
  ! non-rectilinear data
  src_grid(1,:,:) = spread(lon,2,nlat)
  src_grid(2,:,:) = reshape(spread(lat,1,nlon),(/nlon,nlat/))
  
  call nc_read(trim(fname),trim(varname),data,ncid=id_data)

  ! Read in data to be re-gridded
  gridfile = next_arg()
  print *,'Reading in new grid from ',trim(gridfile)
  call nc_open(trim(gridfile), id_grid, writable=.false.)

  ! Initialize the netcdf longitude variable info and load attributes
  call nc_v_init(xgridv,trim(lonname))
  call nc_get_att(id_grid,xgridv,readmeta=.TRUE.)
  call nc_print_attr(xgridv)

  ! Initialize the netcdf latitude variable info and load attributes
  call nc_v_init(ygridv,trim(latname))
  call nc_get_att(id_grid,ygridv,readmeta=.TRUE.)
  call nc_print_attr(ygridv)

  nlonnew = xgridv%dlen(1)
  nlatnew = xgridv%dlen(2)

  print *,nlonnew,' x ',nlatnew

  allocate(newgridx(nlonnew,nlatnew))
  allocate(newgridy(nlonnew,nlatnew))

  call nc_read(trim(gridfile),trim(lonname),newgridx,ncid=id_grid)
  call nc_read(trim(gridfile),trim(latname),newgridy,ncid=id_grid)

  print *,"minval newgridx: ",minval(newgridx)
  print *,"maxval newgridx: ",maxval(newgridx)
 
  ! Need to make this a command line option
  supergrid = .true.
  if (supergrid) then
     ! Supergrid, pull out every second cell
     nx = (nlonnew - 1)/2 + 1
     ny = (nlatnew - 1)/2 + 1
     allocate(dst_grid(2,nx,ny))
     dst_grid(1,:,:) = newgridx(::2,::2)
     dst_grid(2,:,:) = newgridy(::2,::2)
  else
     nx = nlonnew
     ny = nlatnew
     dst_grid(1,:,:) = newgridx
     dst_grid(2,:,:) = newgridy
  end if
  allocate(newdata(nx,ny))

  call regrid(data, src_grid, dst_grid, newdata, data<depth, oceanfrac)
  ! call regrid(data, src_grid, dst_grid, radius, newdata, data<0., oceanfrac)
  ! call regrid(data, src_grid, dst_grid, newdata)

  ! Save result
  print *,"Save result"
  call nc_create(outfile,overwrite=.TRUE.,netcdf4=.TRUE.)
  call nc_write_dim(outfile,"lon",x=dst_grid(1,:,1))
  call nc_write_dim(outfile,"lat",x=dst_grid(2,1,:))
  call nc_write(outfile,"geolon_uv",dst_grid(1,:,:),dim1="lon",dim2="lat",long_name="uv longitude",units="degrees_E")
  call nc_write(outfile,"geolat_uv",dst_grid(2,:,:),dim1="lon",dim2="lat",long_name="uv latitude",units="degrees_N")
  call nc_write(outfile,varname,newdata(:,:),dim1="lon",dim2="lat",missing_value=-1e30)

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Regrid data on new grid using kd-tree'
    write(stderr,*)
    write(stderr,*) 'Usage: treegrid [--help] data newgrid'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*)

  end subroutine usage

end program treegrid

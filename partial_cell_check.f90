program partial_cell_check

  ! Check for think partial cells in topography/bathymetry

  use ncio, only: ncvar, nc_read, nc_open, nc_close, nc_create, nc_write_dim, nc_write, nc_get_att, &
       nc_size, nc_print_attr, nc_v_init, nc_check, nc_write_attr, nc_exists_attr, nc_read_attr, nc_exists_var
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use iso_varying_string
  use string_functions, only: join
  use file_functions, only: exists, freeunit, stderr, stdout, open
  use precision

  use netcdf, only: nf90_put_var

  implicit none

  character(len=2000) :: topogfile, gridfile, fname, outfile
  character(len=200) :: vname, varname
  character(len=20000) :: buffer, cmd

  real(rd_kind), allocatable, target :: data(:,:)

  real(rd_kind), allocatable :: vgridnew(:)
  real(rd_kind), allocatable :: src_grid(:,:,:), dst_grid(:), minthickness(:)
  real, allocatable :: lon(:), lat(:)

  integer :: error, id_data, id_grid, nvgridnew, nlat, nlon, nlatnew, nlonnew, nz
  integer :: i, j, k, increased = 0, decreased = 0

  type(ncvar) :: v, vgridv

  logical :: supergrid

  real(rd_kind) :: radius, thin, depth, diff, mindepth

  real(rd_kind), parameter :: smallvalue = 0.00001

  logical :: debug, have_dim_vars = .false.
  
  type (varying_string) :: myoptions(7)

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  myoptions(1) = 'help'
  myoptions(2) = 'varname'
  myoptions(3) = 'vertname'
  myoptions(4) = 'thinfrac'
  myoptions(5) = 'mindepth'
  myoptions(6) = 'outfile'
  myoptions(7) = 'debug'

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

  debug = option_exists('debug')

  if (option_exists('varname')) then
     ! We have specified variable name
     if (.NOT. has_value('varname')) then
        write(stderr,*) 'Option var must specify a value!'
        call usage
        stop
     end if
     varname = ""
     varname = get_value('varname')
  else
     varname = "depth"
  end if

  if (option_exists('vertname')) then
     ! We have specified variable name for vertical grid variable
     if (.NOT. has_value('vertname')) then
        write(stderr,*) 'Option var must specify a value!'
        call usage
        stop
     end if
     vname = ""
     vname = get_value('vertname')
  else
     vname = "zeta"
  end if

  if (option_exists('thinfrac')) then
     ! Specify minimum partial cell depth as fraction of full cell
     if (.NOT. has_value('thinfrac')) then
        write(stderr,*) 'Option thinfrac must specify a value!'
        call usage
        stop
     end if
     thin = get_value('thinfrac')
  else
     thin = 0.2
  end if

  if (option_exists('mindepth')) then
     ! Specify absolute minimum partial cell depth
     if (.NOT. has_value('mindepth')) then
        write(stderr,*) 'Option mindepth must specify a value!'
        call usage
        stop
     end if
     mindepth = get_value('mindepth')
  else
     mindepth = 0.1
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
     outfile = "topog.nc"
  end if
  ! Read in data to be re-gridded
  topogfile = next_arg()

  print *,'Reading in topog file: ',trim(topogfile)

  call nc_open(trim(topogfile), id_data, writable=.false.)

  ! Initialize the netcdf variable info and load attributes
  call nc_v_init(v,trim(varname))
  call nc_get_att(id_data,v,readmeta=.TRUE.)
  ! call nc_print_attr(v)

  nlon = v%dlen(1)
  nlat = v%dlen(2)

  allocate(lon(nlon),lat(nlat),src_grid(2,nlon,nlat),data(nlon,nlat))

  print *,nlon,'x',nlat

  ! Read longitude and latitude variables
  if (nc_exists_var(trim(topogfile), trim(v%dims(1))) .and. nc_exists_var(trim(topogfile), trim(v%dims(1)))) then
     have_dim_vars = .true.
     call nc_read(trim(fname),trim(v%dims(1)),lon,ncid=id_data)
     call nc_read(trim(fname),trim(v%dims(2)),lat,ncid=id_data)
  end if
  call nc_read(trim(fname),trim(varname),data,ncid=id_data)

  ! Read in vertical grid
  gridfile = next_arg()
  print *,'Reading in new vertical grid from ',trim(gridfile)
  call nc_open(trim(gridfile), id_grid, writable=.false.)

  ! Initialize the netcdf vertical grid variable info and load attributes
  call nc_v_init(vgridv,trim(vname))
  call nc_get_att(id_grid,vgridv,readmeta=.TRUE.)
  call nc_print_attr(vgridv)

  nvgridnew = vgridv%dlen(1)

  print *,' Number of new vgrid points ',nvgridnew

  allocate(vgridnew(nvgridnew))

  call nc_read(trim(gridfile),trim(vname),vgridnew,ncid=id_grid)

  print *,"minval nvgridnew: ",minval(vgridnew)
  print *,"maxval nvgridnew: ",maxval(vgridnew)
 
  ! Supergrid, pull out every second cell
  nz = (nvgridnew - 1)/2 + 1
  allocate(dst_grid(nz))
  dst_grid = vgridnew(::2)

  allocate(minthickness(nz-1))
  minthickness = max((dst_grid(2:nz) - dst_grid(1:nz-1))*thin,mindepth)

  print *,'nz=',nz
  do j = 1, nlat
     do i = 1, nlon
        ! Missing values are a value of -1.0000000E+20
        if (data(i,j) <= 0) then
           data(i,j) = 0.0
           cycle
        end if
        zloop: do k = 1, nz
           ! if (data(i,j) /= 0.0) print *,i,j,data(i,j),k,dst_grid(k)
           if (dst_grid(k) > data(i,j)) then
              exit zloop
           end if
        end do zloop
        ! Upon exit of the loop, k is 1 more than the value of the level
        k = k - 1
        if (debug .and. data(i,j) /= 0.0) print *,i,j,data(i,j),k,dst_grid(k)
        ! Cells deeper than new grid
        if (k == nz) then
           if (debug) print *,i,j,data(i,j),k,dst_grid(k)
           data(i,j) = dst_grid(k)
           cycle
        end if
        diff = data(i,j) - dst_grid(k)
        ! Only need to examine cells whose thickness is less than
        ! the minimum thickness. Need a fuzzy comparison otherwise
        ! get silly decimal round off problems
        if (diff > 0. .and. (minthickness(k)-diff) > smallvalue) then
           if (diff/minthickness(k) > 0.5) then
              ! Increase cell depth to minimum thickness
              data(i,j) = dst_grid(k) + minthickness(k)
              increased = increased + 1
           else
              ! Decrease cell depth to zero
              data(i,j) = dst_grid(k)
              decreased = decreased + 1
           end if
        end if
     end do
  end do
  ! where(data<0) data = 0
  ! where(data>5500.) data = -1

  write(*,'(A,I0," (",F0.2,"%)")') 'Total number of partial cells: ',increased+decreased,100.*real(increased+decreased)/real(nlon*nlat)
  write(*,'(A,I0)') 'Increased: ',increased
  write(*,'(A,I0)') 'Decreased: ',decreased

  ! Save result
  print *,"Save result"
  call nc_create(trim(outfile),overwrite=.true.,netcdf4=.TRUE.)
  if (have_dim_vars) then
     call nc_write_dim(trim(outfile),"grid_x_T",x=lon(:),long_name="Nominal Longitude of T-cell center",units="degree_east")
     call nc_write_dim(trim(outfile),"grid_y_T",x=lat(:),long_name="Nominal Latitude of T-cell center",units="degree_north")
     call nc_write(trim(outfile),varname,data(:,:),dim1="grid_x_T",dim2="grid_y_T",long_name="Topographic depth of T-cell",units="meters")
  else
     call nc_write_dim(trim(outfile),v%dims(1),x=(/(i,i=1,nlon)/))
     call nc_write_dim(trim(outfile),v%dims(2),x=(/(i,i=1,nlat)/))
     call nc_write(trim(outfile),varname,data(:,:),dim1=v%dims(1),dim2=v%dims(2),long_name="Topographic depth of T-cell",units="meters")
     ! call nc_write(trim(outfile),varname,data(:,:),long_name="Topographic depth of T-cell",units="meters")
  end if

  call get_command(cmd)

  cmd = date_time_stamp()//' '//trim(cmd)

  if (nc_exists_attr(trim(outfile),"history")) then
     call nc_read_attr(trim(outfile),"history",buffer)
     cmd = trim(cmd)//new_line('a')//trim(buffer)
  end if
  call nc_write_attr(trim(outfile),"history",trim(cmd))

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Check for partial cells in MOM5 topography file'
    write(stderr,*)
    write(stderr,*) 'Usage: partial_cell_check [--help] topography_file vertical_grid_file'
    write(stderr,*)
    write(stderr,*) '  --help     - print this message'
    write(stderr,*) '  --varname  - variable name in topography file (default = depth)'
    write(stderr,*) '  --vertname - variable name in vertical grid file (default = zeta)'
    write(stderr,*) '  --thinfrac - minimum depth as fraction of full cell (default = 0.2)'
    write(stderr,*) '  --mindepth - absolute minimum depth (default = 0.1m)'
    write(stderr,*)
    write(stderr,*) 'The minimum partial depth will be the minimum of thinfrac*full_cell_depth'
    write(stderr,*) 'and mindepth. If the partial depth is less than this value, it will be'
    write(stderr,*) 'rounded up to this value, or down to the full cell depth of the previous'
    write(stderr,*) 'level, whichever is closer'
    write(stderr,*)

  end subroutine usage

  function date_time_stamp() result(stamp)

    character(len=19) :: stamp

    integer :: values(8)

    call date_and_time(values=values)

1   format(I0.2,"/",I0.2,"/",I0.4,X,I0.2,":",I0.2,":",I0.2)
    write(stamp,1) values(3:1:-1),values(5:7)

  end function date_time_stamp

end program partial_cell_check

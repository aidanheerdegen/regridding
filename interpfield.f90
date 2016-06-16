program interpfields

  ! Interpolate fields on to new vertical grid

  use ncio, only: ncvar, nc_read, nc_open, nc_close, nc_create, nc_write_dim, nc_write, nc_get_att, &
       nc_size, nc_print_attr, nc_v_init, nc_check, nc_write_attr, nc_exists_attr, nc_read_attr
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use iso_varying_string
  use string_functions, only: join
  use file_functions, only: exists, freeunit, stderr, stdout, open
  use precision
  use netcdf, only: nf90_inquire, nf90_inquire_variable, nf90_inquire_dimension

  implicit none

  character(len=2000) :: vgridfile, datafile, outfile, outdir, weightsfile
  character(len=200)  :: varname, dimname, depthname, zname
  
  real(rd_kind), allocatable, target :: data(:,:,:,:), newdata(:,:,:,:)
  real(rd_kind), allocatable, target :: data3d(:,:,:)
  real(rd_kind), allocatable :: vgrid(:), vgridsuper(:), lvl(:)

  real, allocatable :: weights(:)
  integer, allocatable :: neighbours(:,:)

  logical, allocatable :: land_mask(:,:)
  logical :: first_field, set_outfile, set_outdir, interp, unlimited

  integer :: id_grid, i, j, k, l, m, nlon, nlat, nlvl, error, nvgrid, nz, id_data, unit
  integer :: nn1, nn2, pos, ii, id_out, ndim, nvar, ntime, dlen

  type(ncvar) :: vgridv, fieldv
  real(rd_kind), allocatable :: lon(:), lat(:), time(:), dim(:)

  type (varying_string) :: myoptions(5)

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  myoptions(1) = 'help'
  myoptions(2) = 'depthname'
  myoptions(3) = 'zaxisname'
  myoptions(4) = 'outfile'
  myoptions(5) = 'outdir'

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

  if (num_args() < 3) then
     write(stderr,*) 'ERROR! Must supply new vgrid file, weights file and data file(s) as command-line arguments'
     call usage
     STOP
  end if

  if (option_exists('depthname')) then
     ! We have specified variable name for vertical grid variable
     if (.NOT. has_value('depthname')) then
        write(stderr,*) 'Option depthname must specify a value!'
        call usage
        stop
     end if
     depthname = ""
     depthname = get_value('depthname')
  else
     depthname = "zeta"
  end if

  if (option_exists('zaxisname')) then
     ! We have specified variable name for vertical grid variable
     if (.NOT. has_value('zaxisname')) then
        write(stderr,*) 'Option zaxisname must specify a value!'
        call usage
        stop
     end if
     zname = ""
     zname = get_value('zaxisname')
  else
     zname = "zaxis_1"
  end if

  if (option_exists('outfile')) then
     ! We have specified outfile name
     if (.NOT. has_value('outfile')) then
        write(stderr,*) 'Option outfile must specify a value!'
        call usage
        stop
     end if
     outfile = ""
     set_outfile = .true.
     outfile = get_value('outfile')
  else
     set_outfile = .false.
  end if

  if (option_exists('outdir')) then
     ! We have specified output directory name
     if (.NOT. has_value('outdir')) then
        write(stderr,*) 'Option outdir must specify a value!'
        call usage
        stop
     end if
     outdir = ""
     set_outdir = .true.
     outdir = get_value('outdir')
  else
     set_outdir = .false.
  end if

  ! Read in vertical grid
  vgridfile = next_arg()
  print *,'Reading in new vertical grid from ',trim(vgridfile)
  call nc_open(trim(vgridfile), id_grid, writable=.false.)

  ! Initialize the netcdf vertical grid variable info and load attributes
  call nc_v_init(vgridv,trim(depthname))
  call nc_get_att(id_grid,vgridv,readmeta=.TRUE.)
  call nc_print_attr(vgridv)

  ! Supercell
  nvgrid = vgridv%dlen(1)
  nz = (nvgrid - 1)/2

  print *,' Number of new vgrid points ',nz
 
  allocate(vgridsuper(nvgrid))
  allocate(vgrid(nz))

  call nc_read(trim(vgridfile),trim(depthname),vgridsuper,ncid=id_grid)

  ! Supergrid, pull out every second cell
  vgrid = vgridsuper(::2)

  ! Read in weights file
  weightsfile = next_arg()
  print *,'Reading in new weights from ',trim(weightsfile)
  unit = open(trim(weightsfile))
  allocate(neighbours(2,nz),weights(nz))
  do i = 1, nz
     read(unit,*) neighbours(:,i),weights(i)
     ! print *,neighbours(:,i),weights(i)
  end do
  close(unit)

  do while (have_args())

     datafile = next_arg()
     
     print *,'Reading in data file: ',trim(datafile)
     
     call nc_open(trim(datafile), id_data, writable=.false.)

     ! Save result

     if (.not. set_outfile) then
        if (.not. set_outdir) then
           pos = min(index(datafile,'.nc',back=.true.)-1,len_trim(datafile))
           outfile = datafile(1:pos)//'_interp.nc'
        else
           pos = index(datafile,'/',back=.true.)+1
           print *,pos,datafile(pos:len_trim(datafile))
           outfile = trim(outdir)//'/'//datafile(pos:len_trim(datafile))
        end if
     end if
     print *,"Save result in ",trim(outfile)
     call nc_create(trim(outfile),overwrite=.true.,netcdf4=.TRUE.)

     ! Grab the number of variables
     call nc_check( nf90_inquire(id_data, ndim, nvar) )

     first_field = .true.

     do i = 1, ndim
        ! Copy all non z depth dimensions to output
        call nc_check( nf90_inquire_dimension(id_data, i, name=dimname, len=dlen) )
        if (dimname == zname) cycle
        if (dimname == 'Time') then
           unlimited = .true.
        else
           unlimited = .false.
        end if
        allocate(dim(dlen))
        call nc_read(trim(datafile),trim(dimname),dim,ncid=id_data)
        print *,i,trim(dimname),shape(dim)
        call nc_write_dim(trim(outfile),dimname,x=dim,unlimited=unlimited)
        deallocate(dim)
     end do

     ! Write updated z depth dimension to output
     call nc_write_dim(trim(outfile),zname,x=real((/(ii,ii=1,nz)/),rd_kind))

     do i = 1, nvar

        ! Grab the varname for this variable id
        call nc_check( nf90_inquire_variable(id_data, i, varname) )
        
        print *,i,trim(varname)

        ! Initialize the netcdf variable info and load attributes
        call nc_v_init(fieldv,trim(varname))
        call nc_get_att(id_data,fieldv,readmeta=.TRUE.)

        ndim = size(fieldv%dims)
        interp = .true.
        if (ndim == 1 .and. fieldv%dims(1) == fieldv%name) then
           print *,'Skipping dimension'
           cycle
        else if (ndim < 4) then
           print *,'Copying this variable, insufficient dimensions: ',size(fieldv%dims)
           interp = .false.
        else if (fieldv%dims(3) /= zname) then
           print *,'Copying this variable, no matching dimension found: ',trim(zname)
           interp = .false.
        end if
        
        call nc_print_attr(fieldv)

        nlon = fieldv%dlen(1)
        nlat = fieldv%dlen(2)

        if (interp) then

           nlvl = fieldv%dlen(3)
           ntime = fieldv%dlen(4)
           
           allocate(newdata(nlon,nlat,nz,ntime),data(nlon,nlat,nlvl,ntime))
           
           print *,nlon,'x',nlat,'x',nz,'x',ntime
           
           call nc_read(trim(datafile),trim(varname),data,ncid=id_data)
           
           do m = 1, ntime
              do l = 1, nz
                 nn1 = neighbours(1,l)
                 nn2 = neighbours(2,l)
                 print *,'time: ',m,'level: ',l
                 do k = 1, nlat
                    do j = 1, nlon
                       ! New value is a weighted sum of the two closest values from original vertical grid
                       newdata(j,k,l,m) = data(j,k,nn1,m) + (data(j,k,nn2,m) - data(j,k,nn1,m))*weights(l)
                       ! print '(3I6,2X,2(F0.4,2X))',l,k,j, data(j,k,l), newdata(j,k,l) 
                    end do
                 end do
              end do
           end do

           print *,fieldv%dims(1:4)
           print *,shape(newdata)
           print *,trim(varname)
           call nc_write(trim(outfile),varname,newdata,dims=fieldv%dims(1:4),long_name=varname)
           ! call nc_write(trim(outfile),varname,newdata,dim1=trim(fieldv%dims(1)),dim2=trim(fieldv%dims(2)),dim3=trim(fieldv%dims(3)))
           ! call nc_write(trim(outfile),varname,newdata)

           deallocate(newdata,data)

        else

           nlvl = 1
           ntime = fieldv%dlen(3)
           
           allocate(data3d(nlon,nlat,ntime))
           
           print *,nlon,'x',nlat,'x',ntime
           
           call nc_read(trim(datafile),trim(varname),data3d,ncid=id_data)
           
           call nc_write(trim(outfile),varname,data3d,dims=fieldv%dims(1:3),long_name=varname)

           deallocate(data3d)

        end if
     
        first_field = .false.

     end do

     call nc_close(id_data)
        
!!$     call get_command(cmd)
!!$
!!$     cmd = date_time_stamp()//' '//trim(cmd)
!!$     
!!$     if (nc_exists_attr(trim(outfile),"history")) then
!!$        call nc_read_attr(trim(outfile),"history",buffer)
!!$        cmd = trim(cmd)//new_line('a')//trim(buffer)
!!$     end if
!!$     call nc_write_attr(trim(outfile),"history",trim(cmd))

  end do

contains
 
  subroutine usage    
  
    write(stderr,*)
    write(stderr,*) 'Regrid data on to new vertical grid (vgrid) using precalculated weights'
    write(stderr,*)
    write(stderr,*) 'Usage: interpfield [--help] vgrid weights datafile(s)'
    write(stderr,*)
    write(stderr,*) '  --help      - print this message'
    write(stderr,*) '  --help      - print this message'
    write(stderr,*) '  --depthname - name of depth variable in vgrid (default: zeta)'
    write(stderr,*)

  end subroutine usage

end program interpfields

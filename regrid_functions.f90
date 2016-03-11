module regrid_functions

  ! Regrid data using kdtree
  use kdtree2_module
  use precision
  use pathfind_functions
  use statistics, only: mean, median

  implicit none

  interface regrid
     module procedure regrid_real_2d_radius
  end interface regrid

contains

  subroutine regrid_real_2d_nearest_neighbour(datain, sourcegrid, destgrid, dataout)
  
    ! datain must be contiguous otherwise pointer remapping below throws error
    real, target, contiguous, intent(in) :: datain(:,:)
    real, intent(in)                     :: sourcegrid(:,:,:)
    real, intent(in)                     :: destgrid(:,:,:)
    real, intent(out)                    :: dataout(:,:)

    real, pointer :: data1d(:)

    type(kdtree2_result) :: results(1)
    type(kdtree2),pointer    :: tree
    
    integer :: nlonin, nlatin, npointsin
    integer :: nlonout, nlatout, npointsout

    real(kdkind) :: pos_src(3,size(datain)), pos_dest(3,size(datain))

    integer :: i, j, ij

    nlonin = size(datain,1)
    nlatin = size(datain,2)
    npointsin = nlonin*nlatin

    nlonout = size(destgrid,2)
    nlatout = size(destgrid,3)
    npointsout = nlonout*nlatout

    ! data1d = reshape(datain,shape(data1d))
    data1d(1:size(datain)) => datain

    call make_tree(sourcegrid, tree, pos_src, sort=.true.)

    ! Cycle through the points of the new grid, find n nearest neighbours, and
    ! apply function to these neighbours, and save result at the new grid point
    ij = 0
    do j = 1, nlatout
       do i = 1,  nlonout
          ij = ij + 1
          pos_dest(:,ij) = lonlat2vector(destgrid(:,i,j))
       end do
    enddo

    print *,"Traverse tree for all destination points"
    ij = 0
    do i = 1, nlatout
       do j = 1,  nlonout
          ij = ij + 1
          call kdtree2_n_nearest(tp=tree,qv=pos_dest(:,ij),nn=1,results=results)
          dataout(j,i) = data1d(results(1)%idx) 
          ! print *,i,j,ij,results(1)%idx
       end do
    end do

  end subroutine regrid_real_2d_nearest_neighbour
    
  subroutine regrid_real_2d_radius(datain, sourcegrid, destgrid, dataout, mask, oceanfraction)
  
    ! datain must be contiguous otherwise pointer remapping below throws error
    real, target, contiguous, intent(in) :: datain(:,:)
    real, intent(in)                     :: sourcegrid(:,:,:)
    real, intent(in)                     :: destgrid(:,:,:)
    real, intent(out)                    :: dataout(:,:)
    logical, intent(in), optional        :: mask(:,:)
    real, intent(in), optional           :: oceanfraction

    logical :: datamask(size(datain))
    real, pointer :: data1d(:)

    type(kdtree2_result), allocatable :: results(:)
    type(kdtree2), pointer :: tree
    
    real, allocatable :: resultsarray(:)
    real(kind=rd_kind) :: total

    integer :: nlonin, nlatin, npointsin
    integer :: nlonout, nlatout, npointsout
    integer :: nresults = 10000, nfound

    real(kdkind) :: pos_src(3,size(datain)), pos_dest(3,size(datain))

    integer :: i, j, k, ij, n, idx

    real :: oceanfrac, maxdelta, radius

    real, parameter :: sqrt2 = sqrt(2.)

    nlonin = size(datain,1)
    nlatin = size(datain,2)
    npointsin = nlonin*nlatin

    nlonout = size(destgrid,2)
    nlatout = size(destgrid,3)
    npointsout = nlonout*nlatout

    ! data1d = reshape(datain,shape(data1d))
    data1d(1:size(datain)) => datain

    if (present(mask)) then
       datamask = reshape(mask,shape(datamask))
       print *,'datamask: ',real(count(datamask))/size(datamask)

       if (present(oceanfraction)) then
          oceanfrac = oceanfraction
       else
          oceanfrac = 0.5
       end if
    else
       datamask = .true.
       oceanfrac = 0.0
    end if

    call make_tree(sourcegrid, tree, pos_src, sort=.true.)

    ! Create 3D vectors for destination grid
    ij = 0
    do j = 1, nlatout
       do i = 1,  nlonout
          ij = ij + 1
          pos_dest(:,ij) = lonlat2vector(destgrid(:,i,j))
          ! if (i == 1 .or. j == 1) print *,i,j,ij,destgrid(:,i,j),pos_dest(:,ij)
       end do
    enddo

    dataout = -1e30

    ! Cycle through the points of the new grid, find n nearest neighbours, and
    ! apply function to these neighbours, and save result at the new grid point
    print *,"Traverse tree for all destination points"
    allocate(results(nresults),resultsarray(nresults))
    ij = 0
    do j = 1, nlatout
       do i = 1,  nlonout

          ! Find the extent of the cell in degrees
          maxdelta = 0.
          if (j < nlatout) then
             maxdelta = abs(destgrid(2,i,j+1)-destgrid(2,i,j))
          end if
          if (j > 1) then
             maxdelta = max(maxdelta,abs(destgrid(2,i,j)-destgrid(2,i,j-1)))
          end if
          ! if (i < nlonout) then
          !    maxdelta = max(maxdelta,abs(destgrid(1,i+1,j)-destgrid(1,i,j)))
          ! end if
          ! if (i > 1) then
          !    maxdelta = max(maxdelta,abs(destgrid(1,i,j)-destgrid(1,i-1,j)))
          ! end if

          ! The distance between two points on a sphere is the angle between
          ! the two points (in radians) * radius. The vectors are normalised
          ! to the radius of earth, so multiply by 1. The square-root of 2
          ! assumes a square grid cell, and so includes all points out to the
          ! corner, and overlaps slightly into the neighbouring cell
          radius = deg2rad*maxdelta

          if (maxdelta > 1.) then
             print *,i,j
             print *,'1 ',abs(destgrid(2,i,j+1)-destgrid(2,i,j))
             print *,'2 ',abs(destgrid(2,i,j)-destgrid(2,i,j-1))
             ! print *,'3 ',abs(destgrid(1,i+1,j)-destgrid(1,i,j))
             ! print *,'4 ',abs(destgrid(1,i,j)-destgrid(1,i-1,j))
             print *,'maxdelta = ',maxdelta
             print *,'radius(km) = ',radius * earth_R / 1000.
          end if

          ij = ij + 1
          ! call kdtree2_n_nearest(tp=tree,qv=pos_dest(:,ij),nn=1,results=results)
          call kdtree2_r_nearest(tree, pos_dest(:,ij), real(radius,kdkind)**2, nfound, nresults, results) 

          if (nfound == 0) then
             ! print *,i,j,ij,nfound,pos_dest(:,ij)
             cycle
          end if
          ! if (nfound == 0) cycle
          n = 0
          do k = 1, min(nfound,nresults)
             idx = results(k)%idx
             if (idx < 1 .or. idx > npointsin) then
                print *,'ERROR => ',results(k)%idx
                cycle
             end if
             if (.not. datamask(idx)) cycle
             n = n + 1
             resultsarray(n) = data1d(idx)
          end do
          if (n > 0 .and. real(n)/real(nfound) > oceanfrac) then
             ! print *,i,j,size(pack(resultsarray(1:nfound),mask=datamask))
             dataout(i,j) = median(resultsarray(1:n))
          end if
       end do
    end do
    deallocate(results, resultsarray)

    print *,'ij',ij
    print *,'npointsin',npointsin

  end subroutine regrid_real_2d_radius
    
end module regrid_functions


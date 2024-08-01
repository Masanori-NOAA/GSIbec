module m_mod_fv3_lola
!$$$ module documentation block
!           .      .    .                                       .
! module:   mod_fv3_lola
!   prgmmr: parrish
!
! abstract:  This module contains routines to interpolate from a single
!             fv3 D grid tile to a rotated lat-lon analysis grid which completely
!             covers the fv3 tile.  Points beyond the fv3 tile are
!             filled with nearest fv3 edge values, but have no actual
!             impact on the analysis.
!
! program history log:
!   2017-02-24  parrish--initial documentation (patterned after
!   mod_fv3_to_a.f90)
!   2017-10-10  wu w - setup interpolation and trnsform coeff in generate_anl_grid
!                      add routines earthuv2fv3, fv3uv2earth, fv3_h_to_ll
!                        fv3_ll_to_h
!   2019-11-01  wu   - add checks in generate_anl_grid to present the mean
!                      longitude correctly to fix problem near lon=0
!   2022-03-01  X.Lu & X.Wang - add functions for HAFS dual ens capability. POC:
!   xuguang.wang@ou.edu
!   
! subroutines included:
!   sub generate_anl_grid
!   sub definecoef_regular_grids
!   sub earthuv2fv3
!   sub fv3uv2earth
!   sub fv3uv2earthens
!   sub fv3_h_to_ll
!   sub fv3_h_to_ll_ens
!   sub fv3_ll_to_h
!   sub rotate2deg 
!   sub unrotate2deg 
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

!       DIAGRAM:  D-Grid layout:
!
!   1                  nx
!   .                   .   (U,H)
!
! 1                     nx +1
! .                       .    (V)

!   U   U   U   U   U   U        + ny +1 (for U)
! V H V H V H V H V H V H V      + ny    (for V,H)
!   U   U   U   U   U   U                            xh(i) = i            dx=1
! V H V H V H V H V H V H V                          xu(i) = i
!   U   U   U   U   U   U                            xv(i) = i-0.5
! V H V H V H V H V H V H V
!   U   U   U   U   U   U                            yh(j) = j            dy=1
! V H V H V H V H V H V H V                          yu(j) = j-0.5
!   U   U   U   U   U   U                            yv(j) = j
! V H V H V H V H V H V H V
!   U   U   U   U   U   U
! V H V H V H V H V H V H V      + 1     (for V,H)
!   U   U   U   U   U   U        + 1     (for U)

! U(nx ,ny +1),V(nx +1,ny ),H(nx ,ny )

  use m_kinds, only: r_kind,i_kind
  implicit none
!
  private
  public :: m_generate_anl_grid
!  public :: generate_anl_grid,fv3_h_to_ll,fv3_ll_to_h,fv3uv2earth,earthuv2fv3
!  public :: fv3dx,fv3dx1,fv3dy,fv3dy1,fv3ix,fv3ixp,fv3jy,fv3jyp,a3dx,a3dx1,a3dy,a3dy1,a3ix,a3ixp,a3jy,a3jyp
!  public :: nxa,nya,cangu,sangu,cangv,sangv,nx,ny,bilinear
!  public :: definecoef_regular_grids,fv3_h_to_ll_ens,fv3uv2earthens
!  public :: fv3dxens,fv3dx1ens,fv3dyens,fv3dy1ens,fv3ixens,fv3ixpens,fv3jyens,fv3jypens,a3dxens,a3dx1ens,a3dyens,a3dy1ens,a3ixens,a3ixpens,a3jyens,a3jypens
!  public :: nxe,nye,canguens,sanguens,cangvens,sangvens

  logical bilinear
  integer(i_kind) nxa,nya,nx,ny
  real(r_kind) ,allocatable,dimension(:,:):: fv3dx,fv3dx1,fv3dy,fv3dy1
  integer(i_kind),allocatable,dimension(:,:)::  fv3ix,fv3ixp,fv3jy,fv3jyp
  real(r_kind) ,allocatable,dimension(:,:):: a3dx,a3dx1,a3dy,a3dy1
  real(r_kind) ,allocatable,dimension(:,:):: cangu,sangu,cangv,sangv
  integer(i_kind),allocatable,dimension(:,:)::  a3ix,a3ixp,a3jy,a3jyp
  integer(i_kind) nxe,nye
  real(r_kind) ,allocatable,dimension(:,:):: fv3dxens,fv3dx1ens,fv3dyens,fv3dy1ens
  integer(i_kind),allocatable,dimension(:,:)::  fv3ixens,fv3ixpens,fv3jyens,fv3jypens
  real(r_kind) ,allocatable,dimension(:,:):: a3dxens,a3dx1ens,a3dyens,a3dy1ens
  real(r_kind) ,allocatable,dimension(:,:):: canguens,sanguens,cangvens,sangvens
  integer(i_kind),allocatable,dimension(:,:)::  a3ixens,a3ixpens,a3jyens,a3jypens
  

contains

subroutine m_generate_anl_grid(nx,ny,grid_lon,grid_lont,grid_lat,grid_latt,gsi_lats,gsi_lons)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    generate_anl_grid
!   prgmmr: parrish
!
! abstract:  define rotated lat-lon analysis grid which is centered on fv3 tile 
!             and oriented to completely cover the tile.
!
! program history log:
!   2017-05-02  parrish
!   2017-10-10  wu   - 1. setup analysis A-grid, 
!                      2. compute/setup FV3 to A grid interpolation parameters
!                      3. compute/setup A to FV3 grid interpolation parameters         
!                      4. setup weightings for wind conversion from FV3 to earth
!   2019-11-01  wu   - add checks to present the mean longitude correctly to fix
!                       problem near lon=0
!
!   2021-08-11   lei - a fix for an upper bound of the dimnsion of  a3jyp 
!   input argument list:
!    nx, ny               - number of cells = nx*ny 
!    grid_lon ,grid_lat   - longitudes and latitudes of fv3 grid cell corners
!    grid_lont,grid_latt  - longitudes and latitudes of fv3 grid cell centers
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use m_kinds, only: r_kind,i_kind
  use constants, only: quarter,one,two,half,zero,deg2rad,rearth,rad2deg
  use gridmod,  only:grid_ratio_fv3_regional, region_lat,region_lon,nlat,nlon
  use gridmod,  only: region_dy,region_dx,region_dyi,region_dxi,coeffy,coeffx
  use gridmod,  only:init_general_transform,region_dy,region_dx 
  use m_mpimod, only: mype
  use egrid2agrid_mod, only: egrid2agrid_parm
  implicit none

  real(r_kind),allocatable,dimension(:)::xbh_a,xa_a,xa_b
  real(r_kind),allocatable,dimension(:)::ybh_a,ya_a,ya_b,yy
  real(r_kind),allocatable,dimension(:,:)::xbh_b,ybh_b
  real(r_kind) dlat,dlon,dyy,dxx,dyyi,dxxi
  real(r_kind) dyyh,dxxh


  integer(i_kind), intent(in   ) :: nx,ny                 ! fv3 tile x- and y-dimensions
  real(r_kind)   , intent(inout) :: grid_lon(nx+1,ny+1)   ! fv3 cell corner longitudes
  real(r_kind)   , intent(inout) :: grid_lont(nx,ny)      ! fv3 cell center longitudes
  real(r_kind)   , intent(inout) :: grid_lat(nx+1,ny+1)   ! fv3 cell corner latitudes
  real(r_kind)   , intent(inout) :: grid_latt(nx,ny)      ! fv3 cell center latitudes
  real(r_kind),intent(inout) :: gsi_lats(:,:),gsi_lons(:,:)

  integer(i_kind) i,j,ir,jr,n
  real(r_kind),allocatable,dimension(:,:) :: xc,yc,zc,gclat,gclon,gcrlat,gcrlon,rlon_in,rlat_in
  real(r_kind),allocatable,dimension(:,:) :: glon_an,glat_an
  real(r_kind) xcent,ycent,zcent,rnorm,centlat,centlon
  real(r_kind) adlon,adlat,alon,clat,clon
  integer(i_kind) nlonh,nlath,nxh,nyh
  integer(i_kind) ib1,ib2,jb1,jb2,jj

  integer(i_kind) nord_e2a
  real(r_kind)gxa,gya

  real(r_kind) x(nx+1,ny+1),y(nx+1,ny+1),z(nx+1,ny+1), xr,yr,zr,xu,yu,zu,rlat,rlon
  real(r_kind) xv,yv,zv,vval
  real(r_kind) cx,cy
  real(r_kind) uval,ewval,nsval
  real(r_kind) diff,sq180
  real(r_kind) d(4),ds
  integer(i_kind) kk,k


  nord_e2a=4
  bilinear=.false.


!   create xc,yc,zc for the cell centers.
  allocate(xc(nx,ny))
  allocate(yc(nx,ny))
  allocate(zc(nx,ny))
  allocate(gclat(nx,ny))
  allocate(gclon(nx,ny))
  allocate(gcrlat(nx,ny))
  allocate(gcrlon(nx,ny))
  do j=1,ny
     do i=1,nx
        xc(i,j)=cos(grid_latt(i,j)*deg2rad)*cos(grid_lont(i,j)*deg2rad)
        yc(i,j)=cos(grid_latt(i,j)*deg2rad)*sin(grid_lont(i,j)*deg2rad)
        zc(i,j)=sin(grid_latt(i,j)*deg2rad)
     enddo
  enddo

!  compute center as average x,y,z coordinates of corners of domain --

  xcent=quarter*(xc(1,1)+xc(1,ny)+xc(nx,1)+xc(nx,ny))
  ycent=quarter*(yc(1,1)+yc(1,ny)+yc(nx,1)+yc(nx,ny))
  zcent=quarter*(zc(1,1)+zc(1,ny)+zc(nx,1)+zc(nx,ny))

  rnorm=one/sqrt(xcent**2+ycent**2+zcent**2)
  xcent=rnorm*xcent
  ycent=rnorm*ycent
  zcent=rnorm*zcent
  centlat=asin(zcent)*rad2deg
  centlon=atan2(ycent,xcent)*rad2deg

!!  compute new lats, lons
  call m_rotate2deg(grid_lont,grid_latt,gcrlon,gcrlat, &
                  centlon,centlat,nx,ny)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  compute analysis A-grid  lats, lons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------obtain analysis grid dimensions nxa,nya
  nxa=1+nint((nx-one)/grid_ratio_fv3_regional)
  nya=1+nint((ny-one)/grid_ratio_fv3_regional)
  !nlat=nya
  !nlon=nxa
  nlat=ny
  nlon=nx
  if(mype==0) print *,'nlat,nlon=nya,nxa= ',nlat,nlon

!--------------------------obtain analysis grid spacing
  dlat=(maxval(gcrlat)-minval(gcrlat))/(ny-1)
  dlon=(maxval(gcrlon)-minval(gcrlon))/(nx-1)
  adlat=dlat*grid_ratio_fv3_regional
  adlon=dlon*grid_ratio_fv3_regional

!-------setup analysis A-grid; find center of the domain
  nlonh=nlon/2
  nlath=nlat/2

  if(nlonh*2==nlon)then
     clon=adlon/two
     cx=half
  else
     clon=adlon
     cx=one
  endif

  if(nlath*2==nlat)then
     clat=adlat/two
     cy=half
  else
     clat=adlat
     cy=one
  endif

!
!-----setup analysis A-grid from center of the domain
!
  if (allocated(rlat_in)) deallocate(rlat_in)
  if (allocated(rlon_in)) deallocate(rlon_in)
  allocate(rlat_in(nlat,nlon))
  allocate(rlon_in(nlat,nlon))
  do j=1,nlon
     alon=(j-nlonh)*adlon-clon
     do i=1,nlat
        rlon_in(i,j)=alon
     enddo
  enddo


  do j=1,nlon
     do i=1,nlat
        rlat_in(i,j)=(i-nlath)*adlat-clat
     enddo
  enddo

  if (allocated(region_dx )) deallocate(region_dx )
  if (allocated(region_dy )) deallocate(region_dy )
  if (allocated(region_dxi )) deallocate(region_dxi )
  if (allocated(region_dyi )) deallocate(region_dyi )
  if (allocated(coeffx )) deallocate(coeffx )
  if (allocated(coeffy )) deallocate(coeffy )
  allocate(region_dx(nlat,nlon),region_dy(nlat,nlon))
  allocate(region_dxi(nlat,nlon),region_dyi(nlat,nlon))
  allocate(coeffx(nlat,nlon),coeffy(nlat,nlon))
  dyy=rearth*adlat*deg2rad
  dyyi=one/dyy
  dyyh=half/dyy
  do j=1,nlon
     do i=1,nlat
        region_dy(i,j)=dyy
        region_dyi(i,j)=dyyi
        coeffy(i,j)=dyyh
     enddo
  enddo

  do i=1,nlat
     dxx=rearth*cos(rlat_in(i,1)*deg2rad)*adlon*deg2rad
     dxxi=one/dxx
     dxxh=half/dxx
     do j=1,nlon
        region_dx(i,j)=dxx
        region_dxi(i,j)=dxxi
        coeffx(i,j)=dxxh
     enddo
  enddo

!
!----------  setup  region_lat,region_lon in earth coord
!
  if (allocated(region_lat)) deallocate(region_lat)
  if (allocated(region_lon)) deallocate(region_lon)
  allocate(region_lat(nlat,nlon),region_lon(nlat,nlon))
  !allocate(glat_an(nlon,nlat),glon_an(nlon,nlat))

  call m_unrotate2deg(region_lon,region_lat,rlon_in,rlat_in, &
                    centlon,centlat,nlat,nlon)

  region_lat=region_lat*deg2rad
  region_lon=region_lon*deg2rad

  do j=1,nlat
     do i=1,nlon
        gsi_lats(i,j)=region_lat(j,i)
        gsi_lons(i,j)=region_lon(j,i)
     enddo
  enddo
  
  !call init_general_transform(glat_an,glon_an)
 
  !deallocate(glat_an,glon_an)

  deallocate( xc,yc,zc,gclat,gclon,gcrlat,gcrlon)
  deallocate(rlat_in,rlon_in)

  deallocate(region_dxi,region_dyi)
  deallocate(coeffx,coeffy)

end subroutine m_generate_anl_grid

end module m_mod_fv3_lola

subroutine m_rotate2deg(rlon_in,rlat_in,rlon_out,rlat_out,rlon0,rlat0,nx,ny)
!$$$  subprogram documentation block
!                .      .    .                                        .
! subprogram:    rotate2deg
!
!   prgmmr: parrish
!
!   Rotate right-handed spherical coordinate to new right-handed spherical
!   coordinate.  The coordinates are latitude (-90 to 90) and longitude.
!   Output for longitude is principle range of atan2d function ( -180 < rlon_out <= 180 )
!
! program history log:
!   2017-05-02  parrish
!
!  Method is as follows:
!  1.  define x,y,z coordinate system with origin at center of sphere,
!      x intersecting sphere at 0 deg N,  0 deg E,
!      y intersecting sphere at 0 deg N, 90 deg E,
!      z intersecting sphere at 90 deg N  (north pole).

!   4 steps:

!   1.  compute x,y,z from rlon_in, rlat_in

!   2.  rotate (x,y,z) about z axis by amount rlon0 -- (x,y,z) --> (xt,yt,zt)

!   3.  rotate (xt,yt,zt) about yt axis by amount rlat0 --- (xt,yt,zt) --> (xtt,ytt,ztt)

!   4.  compute rlon_out, rlat_out from xtt,ytt,ztt

!   This is the desired new orientation, where (0N, 0E) maps to point
!         (rlon0,rlat0) in original coordinate and the new equator is tangent to
!          the original latitude circle rlat0 at original longitude rlon0.
! attributes:
!   langauge: f90
!   machine:
!
!$$$ end documentation block


  use m_kinds, only: r_kind,i_kind
  use constants, only: deg2rad,rad2deg
  implicit none

  integer(i_kind), intent(in   ) :: nx,ny                 ! fv3 tile x- and y-dimensions
  real(r_kind),intent(in   ) :: rlon_in(nx,ny),rlat_in(nx,ny),rlon0,rlat0
  real(r_kind),intent(  out) :: rlon_out(nx,ny),rlat_out(nx,ny)

  real(r_kind) x,y,z, xt,yt,zt, xtt,ytt,ztt
  integer(i_kind) i,j

  do j=1,ny
     do i=1,nx
!   1.  compute x,y,z from rlon_in, rlat_in

        x=cos(rlat_in(i,j)*deg2rad)*cos(rlon_in(i,j)*deg2rad)
        y=cos(rlat_in(i,j)*deg2rad)*sin(rlon_in(i,j)*deg2rad)
        z=sin(rlat_in(i,j)*deg2rad)

!   2.  rotate (x,y,z) about z axis by amount rlon0 -- (x,y,z) --> (xt,yt,zt)

        xt= x*cos(rlon0*deg2rad)+y*sin(rlon0*deg2rad)
        yt=-x*sin(rlon0*deg2rad)+y*cos(rlon0*deg2rad)
        zt=z

!   3.  rotate (xt,yt,zt) about yt axis by amount rlat0 --- (xt,yt,zt) --> (xtt,ytt,ztt)

        xtt= xt*cos(rlat0*deg2rad)+zt*sin(rlat0*deg2rad)
        ytt= yt
        ztt=-xt*sin(rlat0*deg2rad)+zt*cos(rlat0*deg2rad)

!   4.  compute rlon_out, rlat_out from xtt,ytt,ztt

        rlat_out(i,j)=asin(ztt)*rad2deg
        rlon_out(i,j)=atan2(ytt,xtt)*rad2deg
     enddo
  enddo
end subroutine m_rotate2deg

subroutine m_unrotate2deg(rlon_in,rlat_in,rlon_out,rlat_out,rlon0,rlat0,nx,ny)
!$$$  subprogram documentation block
!                .      .    .                                        .
! subprogram:    unrotate2deg
!
!   prgmmr: parrish
!
! abstract:  inverse of rotate2deg.
!
! program history log:
!   2017-05-02  parrish

! attributes:
!   langauge: f90
!   machine:
!
!$$$ end documentation block

  use m_kinds, only: r_kind,i_kind
  use constants, only: deg2rad,rad2deg
  implicit none

  real(r_kind),intent(in   ) :: rlon_out(nx,ny),rlat_out(nx,ny),rlon0,rlat0
  integer(i_kind),intent(in   ) :: nx,ny
  real(r_kind),intent(  out) :: rlon_in(nx,ny),rlat_in(nx,ny)

  real(r_kind) x,y,z, xt,yt,zt, xtt,ytt,ztt
  integer(i_kind) i,j
  do j=1,ny
     do i=1,nx
        xtt=cos(rlat_out(i,j)*deg2rad)*cos(rlon_out(i,j)*deg2rad)
        ytt=cos(rlat_out(i,j)*deg2rad)*sin(rlon_out(i,j)*deg2rad)
        ztt=sin(rlat_out(i,j)*deg2rad)

        xt= xtt*cos(rlat0*deg2rad)-ztt*sin(rlat0*deg2rad)
        yt= ytt
        zt= xtt*sin(rlat0*deg2rad)+ztt*cos(rlat0*deg2rad)

        x= xt*cos(rlon0*deg2rad)-yt*sin(rlon0*deg2rad)
        y= xt*sin(rlon0*deg2rad)+yt*cos(rlon0*deg2rad)
        z= zt

        rlat_in(i,j)=asin(z)*rad2deg
        rlon_in(i,j)=atan2(y,x)*rad2deg
     enddo
  enddo

end subroutine m_unrotate2deg

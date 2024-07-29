module m_gsi

use m_kinds, only: i_kind,r_kind
use m_mpimod, only: npe,mype,mpi_character,gsi_mpi_comm_world
use m_mpimod, only: setworld
use m_mpimod, only: nxpe,nype
use gsimod, only: gsimain_initialize
use gridmod, only: lon2,lat2,lat1,lon1,nsig

use guess_grids, only: nfldsig
use guess_grids, only: ntguessig
use guess_grids, only: gsiguess_init
use guess_grids, only: gsiguess_set

use gsi_4dvar, only: nsubwin
use gsi_4dvar, only: lsqrtb

use bias_predictors, only: predictors,allocate_preds,deallocate_preds

use control_vectors, only: control_vector
use control_vectors, only: allocate_cv,deallocate_cv
use control_vectors, only: assignment(=)
use control_vectors, only: cvars3d
use control_vectors, only: prt_control_norms
use control_vectors, only: inquire_cv
use control_vectors, only: cvars2d, cvars3d

use state_vectors, only: allocate_state,deallocate_state

use hybrid_ensemble_parameters,only: l_hyb_ens
use hybrid_ensemble_isotropic, only: hybens_grid_setup, create_ensemble
use hybrid_ensemble_isotropic, only: bkerror_a_en
use hybrid_ensemble_parameters,only: ntlevs_ens

use ensctl2state_mod, only: ensctl2state_ad, ensctl2state
!use control2state_mod, only: control2state,control2state_ad

use mpeu_util, only: die
use mpeu_util, only: warn

use general_sub2grid_mod, only: sub2grid_info
!use general_sub2grid_mod, only: general_sub2grid_create_info
!use general_sub2grid_mod, only: general_sub2grid_destroy_info

use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: gsi_bundlegetpointer
use gsi_bundlemod, only: gsi_bundleprint
use gsi_bundlemod, only: assignment(=)
use gsi_bundlemod, only: self_add

use constants, only: zero,one
use constants, only: Pa_per_kPa
use constants, only: constoz
use m_berror_stats,only : berror_stats
use jfunc, only: nsclen,npclen,ntclen, set_pointer

implicit none

private
public gsi_init
public gsi_init_guess
public gsi_set_guess
public gsi_cv_space
public gsi_sv_space
public gsi_get_grid
public gsi_set_grid

interface gsi_init
  module procedure init_
end interface gsi_init

interface gsi_init_guess
  module procedure init_guess_
end interface gsi_init_guess

interface gsi_set_guess
  module procedure set_guess2_
  module procedure set_guess3_
end interface gsi_set_guess

interface gsi_get_grid
  module procedure get_hgrid_
end interface gsi_get_grid

interface gsi_set_grid
  module procedure set_vgrid_
end interface gsi_set_grid

interface gsi_cv_space
  module procedure be_cv_space0_
  module procedure be_cv_space1_
end interface gsi_cv_space

interface gsi_sv_space
  module procedure be_sv_space0_
  module procedure be_sv_space1_
end interface gsi_sv_space

logical,save :: gsi_initialized_ = .false.
logical,save :: gsi_iamset_ = .false.

character(len=*), parameter :: myname ="m_gsibec"
contains
  subroutine init_(cv,vgrid,bkgmock,nmlfile,befile,layout,&
                   jouter,inymd,inhms,&
                   comm)

  logical, intent(out) :: cv
  logical, optional, intent(in)  :: vgrid
  logical, optional, intent(out) :: bkgmock
  character(len=*),optional,intent(in) :: nmlfile
  character(len=*),optional,intent(in) :: befile
  integer,optional,intent(in) :: layout(2) ! 1=nx, 2=ny
  integer,optional,intent(out):: jouter
  integer,optional,intent(in) :: inymd(:),inhms(:)
  integer,optional,intent(in) :: comm
  integer :: jouter_def

  character(len=*), parameter :: myname_=myname//"init_"
  type(sub2grid_info) :: sg
  integer :: ier
  logical :: already_init_mpi

  jouter_def = 1
  if(present(jouter)) then
    jouter=jouter_def
  endif
  if(jouter_def > 1) then
     if (mype==0) call warn(myname_,': already initialized, skipping ...')
    return
  endif

  ier=0
  call mpi_initialized(already_init_mpi,ier)
  if(ier/=0) call die(myname,'mpi_initialized(), ieror =',ier)
  if(.not.already_init_mpi) then
     call mpi_init(ier)
     if(ier/=0) call die(myname,'mpi_init(), ier =',ier)
  endif
  call setworld(comm=comm)
  call mpi_comm_size(gsi_mpi_comm_world,npe,ier)
  call mpi_comm_rank(gsi_mpi_comm_world,mype,ier)
  if (present(layout)) then
     nxpe=layout(1)
     nype=layout(2)
     write(6,*)"nxpe,nype=",nxpe,nype
  endif
  if (present(befile)) then
     call befname_(befile,0)
  endif

  call gsimain_initialize

  nfldsig=1
  ntguessig=1
!  if (present(inymd) .and. present(inhms)) then
!      if (size(inymd)/=ntlevs_ens) then
!         print *, 'nymd,ntlevs ', size(inymd), ntlevs_ens
!         call die(myname,'inconsistent number of time slots ier =',99)
!      endif
!      allocate(nymd(ntlevs_ens),nhms(ntlevs_ens))
!      nymd = inymd
!      nhms = inhms
!      nfldsig=ntlevs_ens
!      ntguessig=(ntlevs_ens+1)/2
!  endif

  call set_(vgrid=vgrid)
  !if(l_hyb_ens) then
  !  call hybens_grid_setup()
  !  call create_ensemble
  !endif

!  call set_pointer_()  ! This is called in create_jfunc
!  call set_pointer  ! This is called in create_jfunc

! create subdomain/grid indexes 
! call general_sub2grid_create_info(sg,0,nlat,nlon,nsig,1,.false.)
! istart=sg%istart
! jstart=sg%jstart
! call general_sub2grid_destroy_info(sg)

  !cv = .true. 
  cv = .false. 
  if (present(bkgmock) ) then
    bkgmock = .false.
  endif
  gsi_initialized_=.true.

  end subroutine init_
!--------------------------------------------------------
  subroutine init_guess_
  call gsiguess_init(.false.)
! call gsiguess_bkgcov_init()  ! not where I want for this to be
  end subroutine init_guess_
!--------------------------------------------------------
  subroutine set_guess2_(varname,islot,var)
  character(len=*),intent(in) :: varname
  integer(i_kind),intent(in) :: islot
  real(r_kind),intent(in) :: var(:,:)
  call gsiguess_set(varname,islot,var)
  end subroutine set_guess2_
!--------------------------------------------------------
  subroutine set_guess3_(varname,islot,var)
  character(len=*),intent(in) :: varname
  integer(i_kind),intent(in) :: islot
  real(r_kind),intent(in) :: var(:,:,:)
  call gsiguess_set(varname,islot,var)
  end subroutine set_guess3_
!--------------------------------------------------------
  subroutine be_cv_space0_

  type(control_vector) :: gradx,grady

! apply B to vector: all in control space

! allocate vectors
  call allocate_cv(gradx)
  call allocate_cv(grady)
  gradx=zero
  grady=zero

  call set_silly_(gradx%step(1))
  call gsi2model_units_ad_(gradx%step(1))

  call bkerror(grady)
               
  !if (l_hyb_ens) then
  !   call bkerror_a_en(grady)
  !endif

  !if(bkgv_write_cv/='null') &
  !call write_bundle(grady%step(1),bkgv_write_cv)

  call gsi2model_units_(grady%step(1))

! clean up
  call deallocate_cv(gradx)
  call deallocate_cv(grady)

  end subroutine be_cv_space0_

  subroutine be_cv_space1_(gradx,internalcv,bypassbe)

  type(control_vector)  :: gradx
  logical,optional,intent(in) :: internalcv
  logical,optional,intent(in) :: bypassbe
  !type(gsi_enperts),optional,intent(in) :: epts

  type(control_vector) :: grady

  logical bypassbe_

  bypassbe_ = .false.
  if (present(bypassbe)) then
     if (bypassbe) bypassbe_ = .true.
  endif

! apply B to vector: all in control space
  if (present(internalcv)) then
     if(internalcv) call set_silly_(gradx%step(1))
  endif

! convert model units to gsi
  call gsi2model_units_ad_(gradx%step(1))

! allocate vectors
  call allocate_cv(grady)

  if (bypassbe_) then
     grady=gradx
  else
     grady=gradx
     call bkerror(grady) 
     !if (l_hyb_ens) then
     !   call ensemble_forward_model_ad(gradx%step(1),gradx%aens(1,:,:),1)
     !   call bkerror_a_en(grady)
     !   call ensemble_forward_model(grady%step(1),grady%aens(1,:,:),1)
     !endif
  endif

  !if(bkgv_write_cv/='null') &
  !call write_bundle(grady%step(1),bkgv_write_cv)

! return result in input vector
  gradx=grady

! convert units back to model units
  call gsi2model_units_(gradx%step(1))

! clean up
  call deallocate_cv(grady)

  end subroutine be_cv_space1_
!--------------------------------------------------------
  subroutine be_sv_space0_

  type(gsi_bundle), allocatable :: mval(:)
  integer ii

! start work space
  allocate(mval(nsubwin))
  do ii=1,nsubwin
      call allocate_state(mval(ii))
      mval(ii) = zero
  end do

  call be_sv_space1_(mval,internalsv=.true.)

  do ii=nsubwin,1,-1
      call deallocate_state(mval(ii))
  end do
  deallocate(mval)

  end subroutine be_sv_space0_
!--------------------------------------------------------
  subroutine be_sv_space1_(mval,internalsv,bypassbe)

  type(gsi_bundle) :: mval(nsubwin)
  logical,optional,intent(in) :: internalsv
  logical,optional,intent(in) :: bypassbe
  !type(gsi_enperts),optional,intent(in) :: epts

  character(len=*), parameter :: myname_ = myname//'*be_sv_space1_'
  type(gsi_bundle),allocatable :: eval(:)
  type(control_vector) :: gradx,grady
  type(predictors)     :: sbias
  logical bypassbe_
  integer ii,ier

  if (nsubwin/=1) then
     if(ier/=0) call die(myname,'cannot handle this nsubwin =',nsubwin)
  endif
! if (ntlevs_ens/=1) then
!    if(ier/=0) call die(myname,'cannot handle this ntlevs_ens =',ntlevs_ens)
! endif
  allocate(eval(ntlevs_ens))

  bypassbe_ = .false.
  if (present(bypassbe)) then
     if (bypassbe) bypassbe_ = .true.
  endif

! start work space
  !if (l_hyb_ens) then
  !   do ii=1,ntlevs_ens
  !     call allocate_state(eval(ii))
  !  end do
  !endif
  call allocate_preds(sbias)
  call allocate_cv(gradx)
  call allocate_cv(grady)
  gradx=zero
  grady=zero

! get test vector (mval)
  if (present(internalsv)) then
     if (internalsv) call set_silly_(mval(1))
  endif

! convert from model to gsi units
  call gsi2model_units_ad_(mval(1))

  !if (l_hyb_ens) then
  !   eval(1)=mval(1)
  !   call ensctl2state_ad(eval,mval(1),gradx)
  !endif
  call control2state_ad(mval,sbias,gradx)

! apply B to input (transformed) vector
  if (bypassbe_) then
    grady=gradx
  else
    grady=gradx
    call bkerror(grady) 
    !if (l_hyb_ens) then
    !   call bkerror_a_en(grady)
    !endif
  endif

  call control2state(grady,mval,sbias)
  !if (l_hyb_ens) then
  !   call ensctl2state(grady,mval(1),eval)
  !   mval(1)=eval(1)
  !end if

! if so write out fields from gsi (in GSI units)
!  if(bkgv_write_sv/='null') &
!  call write_bundle(mval(1),bkgv_write_sv)

! convert from gsi to model units
  call gsi2model_units_(mval(1))

! clean up work space
  call deallocate_cv(gradx)
  call deallocate_cv(grady)
  call deallocate_preds(sbias)
  !if (l_hyb_ens) then
  !   do ii=ntlevs_ens,1,-1
  !     call deallocate_state(eval(ii))
  !  end do
  !endif
  deallocate(eval)

  end subroutine be_sv_space1_
!--------------------------------------------------------
  subroutine befname_ (fname,root)
  implicit none
  character(len=*),intent(in) :: fname
  integer, intent(in) :: root
  character(len=*), parameter :: myname_ = myname//"*befname"
  integer ier,clen
  if(mype==root) then
    write(6,'(3a)') myname_, ": reading B error-coeffs from ", trim(fname)
    berror_stats = trim(fname)
  endif
  clen=len(berror_stats)
  call mpi_bcast(berror_stats,clen,mpi_character,root,gsi_mpi_comm_world,ier)
  end subroutine befname_
!--------------------------------------------------------
  subroutine set_(vgrid)

   use constants, only: pi,one,half,rearth
   use m_mpimod, only: mype
   use gridmod, only: nlon,nlat
   use gridmod, only: rlats,rlons,wgtlats
   use gridmod, only: coslon,sinlon
   use gridmod, only: rbs2
   use gridmod, only: sp_a
   use gridmod, only: create_grid_vars
   use gridmod, only: use_sp_eqspace
   use m_gridmod, only: gridmod_vgrid
   use compact_diffs, only: cdiff_created
   use compact_diffs, only: cdiff_initialized
   use compact_diffs, only: create_cdiff_coefs
   use compact_diffs, only: inisph
!  use mp_compact_diffs_mod1, only: init_mp_compact_diffs1
!  use compact_diffs, only: uv2vordiv
   implicit none
   logical,optional :: vgrid

   real(r_kind) :: dlat,dlon,pih
   integer i,j,i1,ifail

   if (gsi_iamset_ ) return

   call create_grid_vars()
   ifail=0
   if(.not.allocated(rlons)) ifail = 1
   if(.not.allocated(rlats)) ifail = 1
   if(ifail/=0) call die('init','dims not alloc', 99)
   call gengrid_vars
   if(present(vgrid)) then
     if(vgrid) call gridmod_vgrid(mype)
   endif
!  call init_mp_compact_diffs1(nsig+1,mype,.false.)
   gsi_iamset_ = .true.
  end subroutine set_
!--------------------------------------------------------
  subroutine get_hgrid_ (gsi_lats,gsi_lons) ! for now: redundant routine
  use constants, only: init_constants_derived, init_constants
  use constants, only: pi,one,two,half,rad2deg
  use mpeu_util, only: die
  use m_gsi_rfv3io_mod, only: gsi_rfv3io_get_grid_specs
  implicit none
   logical :: regional
   integer :: ierr
   real(r_kind),intent(inout) :: gsi_lats(:,:),gsi_lons(:,:)

   regional = .true.
!  Unfortunately, need to make sure basic constants are initialized
   call init_constants_derived
   call init_constants(regional) 
   call gsi_rfv3io_get_grid_specs(gsi_lats,gsi_lons,ierr)

  end subroutine get_hgrid_
!--------------------------------------------------------
  subroutine set_vgrid_(myid,akbk)
  use m_gridmod, only: gridmod_vgrid
  implicit none
  integer(i_kind),intent(in)  :: myid
  character(len=*),intent(in) :: akbk
  call gridmod_vgrid(myid,akbk)
  end subroutine set_vgrid_
!--------------------------------------------------------
  subroutine set_pointer_
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    set_pointer
!   prgmmr: treadon          org: np23                date: 2004-07-28
!
! abstract: Set length of control vector and other control 
!           vector constants
!
! program history log:
!   2004-07-28  treadon
!   2006-04-21  kleist - include pointers for more time tendency arrays
!   2008-12-04  todling - increase number of 3d fields from 6 to 8 
!   2009-09-16  parrish - add hybrid_ensemble connection in call to
!   setup_control_vectors
!   2010-03-01  zhu     - add nrf_levb and nrf_leve, generalize nval_levs
!                       - generalize vector starting points such as nvpsm, nst2,
!                       and others
!   2010-05-23  todling - remove pointers such as nvpsm, nst2, and others (intro
!   on 10/03/01)
!                       - move nrf_levb and nrf_leve to anberror where they are
!                       needed
!   2010-05-29  todling - generalized count for number of levels in state
!   variables
!   2013-10-22  todling - revisit level count in view of changes to bundle
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$
    use gridmod, only: latlon11,latlon1n,nsig,lat2,lon2
    use gridmod, only: nlat,nlon
    use gridmod, only: nnnn1o
    use state_vectors, only: ns2d,levels
    use constants, only : max_varname_length
    use bias_predictors, only: setup_predictors
    use control_vectors, only: nc2d,nc3d
    use control_vectors, only: setup_control_vectors
    use state_vectors, only: setup_state_vectors
    use hybrid_ensemble_parameters, only: l_hyb_ens,n_ens,generate_ens,grd_ens
    use jfunc, only: nval_lenz
    use jfunc, only: nclenz
    use radinfo, only: npred,jpch_rad
    use pcpinfo, only: npredp,npcptype
    implicit none

    character(len=*),parameter :: myname_=myname//'set_pointer_'
    integer(i_kind) n_ensz,nval_lenz_tot,nval_lenz_enz

    integer(i_kind) nvals_levs,nval_len
    integer(i_kind) nvals_len,nval_levs
    integer(i_kind) nclen,nrclen,nval2d

    if(lsqrtb) then
      call die(myname_,': cannot handle lsqrtb=.t.',999)
    endif

    nvals_levs=max(0,ns2d)+sum(levels)
    nvals_len=nvals_levs*latlon11

    nval_levs=max(0,nc3d)*nsig+max(0,nc2d)
    nval_len=nval_levs*latlon11
    !if(l_hyb_ens) then
    !   nval_len=nval_len+n_ens*grd_ens%nsig*grd_ens%latlon11
    !end if
    nsclen=npred*jpch_rad
    npclen=npredp*npcptype
    ntclen=0
    nclen=nsubwin*nval_len+nsclen+npclen+ntclen
    nrclen=nsclen+npclen+ntclen
  
    n_ensz=0
    nval_lenz_enz=0
    !if(l_hyb_ens.and.generate_ens) then
    !   call set_sqrt_2dsize_(nval2d)
    !   nval_lenz=nval2d*nnnn1o
    !   nval_lenz_tot=nval_lenz
    !   nclenz=nsubwin*nval_lenz_tot+nsclen+npclen+ntclen
    !else
       nval2d=latlon11
    !end if

    CALL setup_control_vectors(nsig,lat2,lon2,latlon11,latlon1n, &
                               nsclen,npclen,ntclen,nclen,nsubwin,&
                               nval_len,lsqrtb,n_ens, &
                               nval_lenz_enz)
    CALL setup_predictors(nrclen,nsclen,npclen,ntclen)
    CALL setup_state_vectors(latlon11,latlon1n,nvals_len,lat2,lon2,nsig)

  end subroutine set_pointer_
!--------------------------------------------------------

  subroutine set_silly_(bundle)
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  implicit none
  type(gsi_bundle) bundle
  character(len=*), parameter :: myname_ = myname//'*set_silly_'
  real(r_kind),pointer :: ptr3(:,:,:)=>NULL()
  real(r_kind),pointer :: ptr2(:,:)=>NULL()
  integer xloc, yloc, k,iset
  integer ier
  integer zex(4),zex072(4), zex127(4)
  real(r_kind) :: val
  character(len=2) :: var
  character(len=80):: ifname(1)
  character(len=80):: ofname
! logical :: fexist
!ifname = 'xinc.eta.nc4'
!ofname = 'outvec_bkgtest'
!inquire(file=ifname(1),exist=fexist)
!if(.not.fexist) then
!  call die ('main',': fishy', 99)
!endif
!
  if (mod(mype,3) /= 0) return
!
!            sfc  ~500  ~10   ~1
  zex072 = (/  1,   23,  48,  58 /)
  zex127 = (/  1,   53, 106, 116 /)
  iset=-1
  if (nsig==72) then
    zex=zex072
    iset=1
  endif
  if (nsig==127) then
    zex=zex127
    iset=1
  endif
  xloc=min(20,lat2)
  yloc=min(20,lon2)
  val=one
  var='sf'
  var='q'
  var='t'
  var='tv'
  if (iset<0) call die(myname_,'no input set',99)
  call gsi_bundlegetpointer(bundle,trim(var),ptr3,ier)
  if(ier==0) then
     if(var=='sf' .or. var=='vp') then
       val=val*1e-5
     endif
     do k=1,size(zex)
        ptr3(xloc,yloc,zex(k)) = val
     enddo
     if (mype==0) print *, myname_, ': var= ', trim(var)
     return
  endif
  if(var == 'tv') then
     call gsi_bundlegetpointer(bundle,trim(var),ptr3,ier)
     if(ier==0) then
        do k=1,size(zex)
           ptr3(xloc,yloc,zex(k)) = val
        enddo
        if (mype==0) print *, myname_, ': var= ', trim(var)
        return
     endif
  endif
  if (var == 'ps') then
     call gsi_bundlegetpointer(bundle,'ps',ptr2,ier)
     if(ier==0) then
        ptr2(xloc,yloc) = 100.
        if (mype==0) print *, myname_, ': var= ', 'ps(Pa)'
        return
     endif
  endif
  end subroutine set_silly_
!--------------------------------------------------------
  subroutine gsi2model_units_(bundle)
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  implicit none
  type(gsi_bundle) bundle
  real(r_kind),pointer :: ptr2(:,:)  =>NULL()
  real(r_kind),pointer :: ptr3(:,:,:)=>NULL()
  integer ier
  call gsi_bundlegetpointer(bundle,'ps',ptr2,ier)
  if(ier==0) then
     ptr2 = ptr2 * Pa_per_kPa
  endif
  call gsi_bundlegetpointer(bundle,'oz',ptr3,ier)
  if(ier==0) then
     ptr3 = ptr3 * constoz
  endif
  end subroutine gsi2model_units_
!--------------------------------------------------------
  subroutine gsi2model_units_ad_(bundle)
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  implicit none
  type(gsi_bundle) bundle
  real(r_kind),pointer :: ptr2(:,:)  =>NULL()
  real(r_kind),pointer :: ptr3(:,:,:)=>NULL()
  integer ier
  call gsi_bundlegetpointer(bundle,'ps',ptr2,ier)
  if(ier==0) then
     ptr2 = ptr2 * Pa_per_kPa
  endif
  call gsi_bundlegetpointer(bundle,'oz',ptr3,ier)
  if(ier==0) then
     ptr3 = ptr3 * constoz
  endif
  end subroutine gsi2model_units_ad_
end module m_gsi

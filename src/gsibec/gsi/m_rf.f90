module m_rf
use m_kinds, only: r_kind, i_kind
use constants, only: zero
use m_mpimod, only: mype

use berror, only: create_berror_vars_reg
use berror, only: destroy_berror_vars_reg
!use berror, only: final_rftable
use berror, only: fut2ps,cwcoveqqcov

use balmod, only: create_balance_vars_reg
use balmod, only: destroy_balance_vars_reg
use balmod, only: prebal_reg

use smooth_polcarf, only: destroy_smooth_polcas
use gridmod, only: nlon,nlat,lon2,lat2,lon2,nsig
use m_berror_stats_reg, only: berror_get_dims_reg
!use m_berror_stats, only: berror_init
!use m_berror_stats, only: berror_final
use jfunc, only: create_jfunc

use hybrid_ensemble_parameters, only: l_hyb_ens
use hybrid_ensemble_parameters, only: destroy_hybens_localization_parameters
use hybrid_ensemble_isotropic, only: load_ensemble
use hybrid_ensemble_isotropic, only: hybens_localization_setup
!use hybrid_ensemble_parameters, only: gsi_enperts
!use hybrid_ensemble_parameters, only: gsi_destroy_ensemble

use observermod, only: observer_init
use radiance_mod, only: radiance_obstype_init

use mpeu_util, only: getindex
use mpeu_util, only: die

use m_mpimod, only:  nvar_id
implicit none
private
public :: rf_set
!public :: rf_unset

interface rf_set; module procedure set_; end interface
!interface rf_unset; module procedure unset_; end interface

character(len=*), parameter :: myname = 'm_rf'
contains
  subroutine set_
  implicit none
!  type(gsi_enperts) :: epts
!  integer(i_kind), intent(in) :: nymd,nhms
 
  character(len=*), parameter :: mynmae_ = myname//'*set_'
  integer(i_kind) msig,mlat,mlon
  logical good
! Load background error arrays used by recursive filters
!  call berror_get_dims_reg(msig,mlat)
!  good=nlat==mlat.and.nsig==msig  !nlat=
!  if (.not. good) then
!    print *, nlat, mlat, nlon, mlon, nsig, msig
!    call die(myname,': bad dims',99)
!  endif
!  call berror_init(mlat,msig)
!  call radiance_obstype_init
!write(6,*)"Test5"
!  call observer_init()
  call create_jfunc
  call create_balance_vars_reg(mype)
  call create_berror_vars_reg
  call prebal_reg(cwcoveqqcov)
  call prewgt_reg(mype)
! If hybrid covariance
!  if(l_hyb_ens) then
!     call load_ensemble
!     call hybens_localization_setup
!  end if
  end subroutine set_
!  subroutine unset_(epts)
!  implicit none
!  type(gsi_enperts) :: epts
!  if (l_hyb_ens) then
!    call destroy_hybens_localization_parameters
!    call gsi_destroy_ensemble(epts)
!  endif
!  call destroy_smooth_polcas ! the set is called in prewgt - gsi typically has inconsistent set/unset
!  call destroy_berror_vars
!  call destroy_balance_vars
!  call final_rftable         ! out of place as consequence of gsi typically has inconsistent set/unset
!  call berror_final
!  call final_
!  end subroutine unset_
  subroutine final_
  end subroutine final_


end module m_rf

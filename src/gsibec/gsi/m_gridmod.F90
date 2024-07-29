module m_gridmod

! !USES:

  use m_kinds, only: i_byte,r_kind,r_single,i_kind
  use constants, only: kPa_per_Pa,Pa_per_kPa
  use gridmod, only: nsig,ak5,bk5,ck5,tref5
  implicit none

! set default to private
  private
! set subroutines to public
  public :: create_vgrid_vars
  public :: gridmod_vgrid

  interface gridmod_vgrid
     module procedure load_vert_coord_
  end interface

  character(len=*),parameter::myname='gridmod'
contains

  subroutine create_vgrid_vars
    implicit none
    if(.not.allocated(ak5)) &
      allocate(ak5(nsig+1),bk5(nsig+1),ck5(nsig+1),tref5(nsig))
  end subroutine create_vgrid_vars

  subroutine load_vert_coord_(mype,fname)
  use m_set_eta, only: set_eta
  use m_set_eta, only: set_eta_read
  ! ideally, these coordinates should be passed from JEDI
  implicit none
  integer(i_kind),intent(in) :: mype
  character(len=*),optional,intent(in) :: fname
  character(len=*),parameter:: myname_=myname//'*load_vert_coord_'
  integer ks,ifail,ier
  real(r_kind) :: ptop,pint
  ifail=0 
  if(.not.allocated(ak5)) ifail=1
  if(.not.allocated(bk5)) ifail=1
  if(ifail/=0) then
     call create_vgrid_vars()
     ifail=0
  endif
  ! Expect FV3 levels/orientation/units
  if (present(fname)) then
    if (trim(fname)=='/dev/null') then
!write(6,*)"Test4"
       call set_eta (nsig, ks, ptop, pint, ak5, bk5)
    else
!write(6,*)"Test5"
       call set_eta_read(trim(fname),ak5,bk5,ier,myid=mype)
!if(mype==0)write(6,*)"ak5=",ak5
!if(mype==0)write(6,*)"bk5=",bk5
    endif
  else
!write(6,*)"Test6"
    !call set_eta (nsig, ks, ptop, pint, ak5, bk5)
    call set_eta_read('akbk61.nc',ak5,bk5,ier,myid=mype)
!if(mype==0)write(6,*)"ak5=",ak5
!if(mype==0)write(6,*)"bk5=",bk5
  endif
  ! Reorient and adjust units for GSI
  ak5=kPa_per_Pa*ak5
  ak5=ak5(nsig+1:1:-1)
  bk5=bk5(nsig+1:1:-1)
  end subroutine load_vert_coord_

end module m_gridmod


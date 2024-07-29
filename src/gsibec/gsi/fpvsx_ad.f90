subroutine fpvsx_ad( t, es, t_ad, es_ad, adjoint )
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    fpvsx_ad     forward and adjoint model for saturation vapor pressure
!     prgmmr:    treadon     org: np23                date: 2003-12-18
!
! abstract:  This subroutine contains the forward and ajoint models for the
!            calculation of saturation vapor pressure.  
!
! program history log:
!   03-12-18  treadon - initial routine
!   04-06-14  treadon - reformat documenation
!
!   input argument list:
!     t       - temperature
!     t_ad    - partial derivative of vapor pressure with respect to temperature 
!     es_ad   - vapor pressure perturbation
!     adjoint - logical flag (.false.=forward model only, .true.=forward and ajoint)
!
!   output argument list:
!     es
!     t_ad    - partial derivative of vapor pressure with respect to temperature 
!     es_ad   - vapor pressure perturbation
!
! remarks:
!    The adjoint portion of this routine was generated by the 
!    Tangent linear and Adjoint Model Compiler,  TAMC 5.3.0
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
!==============================================
! all entries are defined explicitly
!==============================================
  use m_kinds, only: r_kind
  use constants, only: zero, one, tmix, xai, xbi, xa, xb, ttp, psatk
  implicit none

!==============================================
! define arguments
!==============================================
  logical     ,intent(in   ) :: adjoint
  real(r_kind),intent(inout) :: es_ad
  real(r_kind),intent(inout) :: t_ad
  real(r_kind),intent(  out) :: es
  real(r_kind),intent(in   ) :: t

!==============================================
! define local variables
!==============================================
  real(r_kind) tr_ad
  real(r_kind) w_ad
  real(r_kind) tr
  real(r_kind) w

!----------------------------------------------
! RESET LOCAL ADJOINT VARIABLES
!----------------------------------------------
  tr_ad = zero
  w_ad = zero
!----------------------------------------------
! ROUTINE BODY
!----------------------------------------------
!----------------------------------------------
! FUNCTION AND TAPE COMPUTATIONS
!----------------------------------------------

  tr = ttp/t
  if (t >= ttp) then
     es = psatk*tr**xa*exp(xb*(one-tr))
  else if (t < tmix) then
     es = psatk*tr**xai*exp(xbi*(one-tr))
  else
     w = (t-tmix)/(ttp-tmix)
     es = w*psatk*tr**xa*exp(xb*(one-tr))+(one-w)*psatk*tr**xai* &
          exp(xbi*(one-tr))
  endif
  if (.not.adjoint) return

!----------------------------------------------
! ADJOINT COMPUTATIONS
!----------------------------------------------
  if (t >= ttp) then
     tr_ad = tr_ad+es_ad*((-(psatk*tr**xa*xb*exp(xb*(one-tr))))+psatk*xa* &
          tr**(xa-one)*exp(xb*(one-tr)))
     es_ad = zero
  else if (t < tmix) then
     tr_ad = tr_ad+es_ad*((-(psatk*tr**xai*xbi*exp(xbi*(one-tr))))+psatk* &
          xai*tr**(xai-one)*exp(xbi*(one-tr)))
     es_ad = zero
  else
     tr_ad = tr_ad+es_ad*((-(w*psatk*tr**xa*xb*exp(xb*(one-tr))))+w* &
          psatk*xa*tr**(xa-one)*exp(xb*(one-tr))-(one-w)*psatk*tr**xai*xbi* &
          exp(xbi*(one-tr))+(one-w)*psatk*xai*tr**(xai-one)*exp(xbi*(one-tr)))
     w_ad = w_ad+es_ad*(psatk*tr**xa*exp(xb*(one-tr))-psatk*tr**xai* &
          exp(xbi*(one-tr)))
     es_ad = zero
     t_ad = t_ad+w_ad/(ttp-tmix)
     w_ad = zero
  endif
  t_ad = t_ad-tr_ad*(ttp/(t*t))
  tr_ad = zero
  
  return
end subroutine fpvsx_ad

subroutine fpvsx_tl( t, es, t_d, es_d )
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    fpvsx_tl     forward and tangent linear model for saturation vapor pressure
!     prgmmr:    kim     org: np23                date: 2012-02-16
!
! abstract:  This subroutine contains the forward and tangent linear models for the
!            calculation of saturation vapor pressure.  
!
! program history log:
!   2012-02-16  kim - initial routine based on Russ Treadon's fpvsx_ad subroutine
!
!$$$
!==============================================
! all entries are defined explicitly
!==============================================
  use m_kinds, only: r_kind
  use constants, only: zero, one, tmix, xai, xbi, xa, xb, ttp, psatk
  implicit none

!==============================================
! define arguments
!==============================================
  real(r_kind),intent(out) :: es_d
  real(r_kind),intent(in) :: t_d
  real(r_kind),intent(  out) :: es
  real(r_kind),intent(in   ) :: t

!==============================================
! define local variables
!==============================================
  real(r_kind) tr_d
  real(r_kind) w_d
  real(r_kind) tr
  real(r_kind) w

!----------------------------------------------
! RESET LOCAL ADJOINT VARIABLES
!----------------------------------------------
  tr_d = zero
  w_d = zero
!----------------------------------------------
! ROUTINE BODY
!----------------------------------------------
!----------------------------------------------
! FUNCTION AND TAPE COMPUTATIONS
!----------------------------------------------

  tr = ttp/t
  tr_d = -ttp*t_d/t**2
 
  if (t >= ttp) then
     es = psatk*tr**xa*exp(xb*(one-tr))
     es_d = xa*psatk*tr_d*tr**(xa-1)*exp(xb*(one-tr)) &
           -xb*tr_d*psatk*tr**xa*exp(xb*(one-tr))
  else if (t < tmix) then
     es = psatk*tr**xai*exp(xbi*(one-tr))
     es_d = psatk*xai*tr_d*tr**(xai-1)*exp(xbi*(one-tr)) &
           -xbi*tr_d*psatk*tr**xai*exp(xbi*(one-tr))
  else
     w = (t-tmix)/(ttp-tmix)
     w_d = t_d/(ttp-tmix)
     es = w*psatk*tr**xa*exp(xb*(one-tr))+(one-w)*psatk*tr**xai* &
          exp(xbi*(one-tr))
     es_d = w_d*psatk*tr**xa*exp(xb*(one-tr)) + w*psatk*xa*tr_d*tr**(xa-1)*exp(xb*(one-tr)) &
           -xb*tr_d*w*psatk*tr**xa*exp(xb*(one-tr)) &
           -w_d*psatk*tr**xai*exp(xbi*(one-tr)) &
           +(one-w)*psatk*xai*tr_d*tr**(xai-1)*exp(xbi*(one-tr)) &
           -xbi*tr_d*(one-w)*psatk*tr**xai*exp(xbi*(one-tr))
  endif
RETURN
END subroutine fpvsx_TL

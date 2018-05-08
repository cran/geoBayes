!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Commentary:
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine spcorr (dist,phi,kappa,n,icf)
  use covfun, only: create_spcor, covmat
  implicit none
  integer, intent(in) :: n, icf
  double precision, intent(in) :: phi, kappa
  double precision, intent(inout) :: dist(n)
  call create_spcor(icf,0)
  call covmat (dist,phi,kappa,n)
end subroutine spcorr

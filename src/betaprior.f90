!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Functions related to the prior for beta
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module betaprior
  implicit none
contains
  subroutine betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)
    implicit none
    integer, intent(in) :: n, p
    double precision, intent(in) :: betm0(p), betQ0(p,p), F(n,p), ssqdf
    logical, intent(out) :: lmxi
    double precision, intent(out) :: modeldfh, xi(n)
    integer i

    lmxi = .false.
    do i = 1, p
      lmxi = betQ0(i,i) .gt. 0d0
      if (lmxi) exit
    end do
    if (lmxi) then
      modeldfh = .5d0*(n + ssqdf)
      xi = matmul(F,betm0)
      lmxi = any(xi .ne. 0d0)
    else
      modeldfh = .5d0*(n - p + ssqdf)
    end if
  end subroutine betapriorz
end module betaprior

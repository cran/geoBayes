!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Functions related to the prior for beta
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module betaprior
  implicit none
contains
  subroutine betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf, &
     offset)
    implicit none
    integer, intent(in) :: n, p
    double precision, intent(in) :: betm0(p), betQ0(p,p), F(n,p), ssqdf, &
       offset(n) 
    logical, intent(out) :: lmxi
    double precision, intent(out) :: modeldfh, xi(n)
    lmxi = (betQ0(1,1) .gt. 0d0) ! Normal prior
    if (lmxi) then
      modeldfh = .5d0*(n + ssqdf)
      xi = matmul(F,betm0) + offset
      lmxi = any(xi .ne. 0d0)
    else ! Uniform prior
      modeldfh = .5d0*(n - p + ssqdf)
      xi = offset
      lmxi = any(xi .ne. 0d0)
    end if
  end subroutine betapriorz
end module betaprior

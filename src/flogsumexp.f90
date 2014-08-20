!!! flogsumexp.f90 --- 
!! 
!! Author: Evangelos Evangelou
!! Created: Thu, 19 Jun, 2014 13:34 (BST)
!! Last-Updated: Thu, 19 Jun, 2014 13:35 (BST)
!!     Update #: 1
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module flogsumexp
contains
  function logrsumexp (x,m,n)
    implicit none
    integer, intent(in) :: m, n
    double precision, intent(in) :: x(m,n)
    double precision logrsumexp(m)
    double precision xmm(m,n), se(m)
    logrsumexp = maxval(x,2)
    xmm = spread(logrsumexp,2,n)
    xmm = x - xmm
    xmm = exp(xmm)
    se = sum(xmm,2)
    logrsumexp = logrsumexp + log(se)
  end function logrsumexp

  function logcsumexp (x,m,n)
    implicit none
    integer, intent(in) :: m, n
    double precision, intent(in) :: x(m,n)
    double precision logcsumexp(n)
    double precision xmm(m,n), se(n)
    logcsumexp = maxval(x,1)
    xmm = spread(logcsumexp,1,m)
    xmm = x - xmm
    xmm = exp(xmm)
    se = sum(xmm,1)
    logcsumexp = logcsumexp + log(se)
  end function logcsumexp
end module flogsumexp

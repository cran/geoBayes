!!! flogsumexp.f90 --- 
!! 
!! Author: Evangelos Evangelou
!! Created: Thu, 19 Jun, 2014 13:34 (BST)
!! Last-Updated: Tue, 13 Jan, 2015 15:19 (GMT)
!!     Update #: 2
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module flogsumexp
contains
  pure function logsumexpv (x,n)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision logsumexpv
    double precision xmm(n), se
    logsumexpv = maxval(x)
    xmm = x - logsumexpv
    xmm = exp(xmm)
    se = sum(xmm)
    logsumexpv = logsumexpv + log(se)
  end function logsumexpv

  pure function logrsumexp (x,m,n)
    implicit none
    integer, intent(in) :: m, n
    double precision, intent(in) :: x(m,n)
    double precision logrsumexp(m)
    double precision xmm(m,n), se(m)
    integer i, j
    logrsumexp = maxval(x,2)
    do j = 1, n
      do i = 1, m
        xmm(i,j) = x(i,j) - logrsumexp(i)
      end do
    end do
    xmm = exp(xmm)
    se = sum(xmm,2)
    logrsumexp = logrsumexp + log(se)
  end function logrsumexp

  pure function logcsumexp (x,m,n)
    implicit none
    integer, intent(in) :: m, n
    double precision, intent(in) :: x(m,n)
    double precision logcsumexp(n)
    double precision xmm(m,n), se(n)
    integer i, j
    logcsumexp = maxval(x,1)
    do j = 1, n
      do i = 1, m
        xmm(i,j) = x(i,j) - logcsumexp(j)
      end do
    end do
    xmm = exp(xmm)
    se = sum(xmm,1)
    logcsumexp = logcsumexp + log(se)
  end function logcsumexp
end module flogsumexp

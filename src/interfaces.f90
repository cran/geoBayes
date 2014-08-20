!!! interfaces.f90 --- 
!! 
!! Author: Evangelos Evangelou
!! Created: Thu, 19 Jun, 2014 11:34 (BST)
!! Last-Updated: Thu, 19 Jun, 2014 11:34 (BST)
!!     Update #: 1
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Contains general-purpose functions
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module interfaces
  implicit none 
  public
  private logprobt0, logprobt1, logprobnorm0, logprobnorm1, &
     logproblogis0, logproblogis1, probinvs, probinvv, probinvm, &
     logpdfs, logpdfv, logpdfm
  private logprobs, logborps, logprobv, logborpv, logprobm, logborpm 
  private logbxcxs, logbxcxv, logbxcxm, invbxcxs, invbxcxv, invbxcxm
  double precision, parameter :: bigneg = -huge(1d0)

!!!!!!!!!!!!!!!!!!!!!!!!! Correlation functions !!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    elemental double precision function matern (d, kappa)
      double precision, intent(in) :: d, kappa
    end function matern
  end interface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RNG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    subroutine rngini ()
    end subroutine rngini
  end interface

  interface
    subroutine rngend ()
    end subroutine rngend
  end interface

  interface
    double precision function randnorm ()
    end function randnorm
  end interface

  interface
    double precision function randunif ()
    end function randunif
  end interface

  interface
    double precision function randt (d)
      double precision, intent(in) :: d
    end function randt
  end interface

  interface
    double precision function randgamma (shp)
      double precision, intent(in) :: shp
    end function randgamma
  end interface

  interface
    double precision function randchisq (df)
      double precision, intent(in) :: df
    end function randchisq
  end interface


!!!!!!!!!!!!!!!!!!!!!!!!!!! Auxliary functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    elemental double precision function flog1mexp (x)
      double precision, intent(in) :: x
    end function flog1mexp
  end interface

  interface
    elemental double precision function flogexpm1 (x)
      double precision, intent(in) :: x
    end function flogexpm1
  end interface

  interface
    elemental double precision function flog1p (x)
      double precision, intent(in) :: x
    end function flog1p
  end interface

  interface
    elemental double precision function fexpm1 (x)
      double precision, intent(in) :: x
    end function fexpm1
  end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!! Link functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface logprob
    module procedure logprobs, logprobv, logprobm
  end interface logprob

  interface logborp
    module procedure logborps, logborpv, logborpm
  end interface logborp

  interface probinv
    module procedure probinvs, probinvv, probinvm
  end interface probinv

  interface logpdf
    module procedure logpdfs, logpdfv, logpdfm
  end interface logpdf

  interface logbxcx
    module procedure logbxcxs, logbxcxv, logbxcxm
  end interface logbxcx

  interface invbxcx
    module procedure invbxcxs, invbxcxv, invbxcxm
  end interface invbxcx


!!!!!!!!!!!!!!!!!!!!!!!!!!!! Private functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! t distribution
  interface 
    elemental double precision function logprobt0 (q, d)
      double precision, intent(in) :: q, d
    end function logprobt0
  end interface

  interface 
    elemental double precision function logprobt1 (q, d)
      double precision, intent(in) :: q, d
    end function logprobt1
  end interface

  interface
    elemental double precision function quantt (p, d)
      double precision, intent(in) :: p, d
    end function quantt
  end interface

  interface
    elemental double precision function logpdft (x, d)
      double precision, intent(in) :: x, d
    end function logpdft
  end interface

!! normal distribution
  interface 
    elemental double precision function logprobnorm0 (q)
      double precision, intent(in) :: q
    end function logprobnorm0
  end interface

  interface 
    elemental double precision function logprobnorm1 (q)
      double precision, intent(in) :: q
    end function logprobnorm1
  end interface

  interface
    elemental double precision function quantnorm (p)
      double precision, intent(in) :: p
    end function quantnorm
  end interface

  interface
    elemental double precision function logpdfnorm (x)
      double precision, intent(in) :: x
    end function logpdfnorm
  end interface

!! logistic distribution
  interface 
    elemental double precision function logproblogis0 (q)
      double precision, intent(in) :: q
    end function logproblogis0
  end interface

  interface 
    elemental double precision function logproblogis1 (q)
      double precision, intent(in) :: q
    end function logproblogis1
  end interface

  interface
    elemental double precision function quantlogis (p)
      double precision, intent(in) :: p
    end function quantlogis
  end interface

  interface
    elemental double precision function logpdflogis (x)
      double precision, intent(in) :: x
    end function logpdflogis
  end interface

contains

  function geomean (y,n) result(g)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y(n)
    double precision g
    integer i, m
    g = 0d0
    m = 0
    do i = 1, n
      if (y(i) .gt. 0d0) then
        g = g + log(y(i))
        m = m + 1
      end if
    end do
    g = exp(g/dble(m))
  end function geomean

  function logprobs(q,d) result (p)
    double precision, intent(in) :: d, q
    double precision :: p
    if (d > 0d0) then
      p = logprobt1(q,d)
    else if (d < 0d0) then
      p = logproblogis1(q)
    else
      p = logprobnorm1(q)
    end if
  end function logprobs

  function logborps(q,d) result (p)
    double precision, intent(in) :: d, q
    double precision :: p
    if (d > 0d0) then
      p = logprobt0(q,d)
    else if (d < 0d0) then
      p = logproblogis0(q)
    else
      p = logprobnorm0(q)
    end if
  end function logborps

  function probinvs(p,d) result (q)
    double precision, intent(in) :: d, p
    double precision :: q
    if (d > 0d0) then
      q = quantt(p,d)
    else if (d < 0d0) then
      q = quantlogis(p)
    else
      q = quantnorm(p)
    end if
  end function probinvs

  function logpdfs(x,d) result(pdf)
    double precision, intent(in) :: x, d
    double precision :: pdf
    if (d > 0d0) then
      pdf = logpdft(x,d)
    else if (d < 0d0) then
      pdf = logpdflogis(x)
    else
      pdf = logpdfnorm(x)
    end if
  end function logpdfs


  function logprobv(q,d,n) result (p)
    integer, intent(in) :: n
    double precision, intent(in) :: d, q(n)
    double precision :: p(n)
    if (d > 0d0) then
      p = logprobt1(q,d)
    else if (d < 0d0) then
      p = logproblogis1(q)
    else
      p = logprobnorm1(q)
    end if
  end function logprobv

  function logborpv(q,d,n) result (p)
    integer, intent(in) :: n
    double precision, intent(in) :: d, q(n)
    double precision :: p(n)
    if (d > 0d0) then
      p = logprobt0(q,d)
    else if (d < 0d0) then
      p = logproblogis0(q)
    else
      p = logprobnorm0(q)
    end if
  end function logborpv

  function probinvv(p,d,n) result (q)
    integer, intent(in) :: n
    double precision, intent(in) :: d, p(n)
    double precision :: q(n)
    if (d > 0d0) then
      q = quantt(p,d)
    else if (d < 0d0) then
      q = quantlogis(p)
    else
      q = quantnorm(p)
    end if
  end function probinvv

  function logpdfv(x,d,n) result(pdf)
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), d
    double precision :: pdf(n)
    if (d > 0d0) then
      pdf = logpdft(x,d)
    else if (d < 0d0) then
      pdf = logpdflogis(x)
    else
      pdf = logpdfnorm(x)
    end if
  end function logpdfv

  function logprobm(q,d,n1,n2) result (p)
    integer, intent(in) :: n1, n2
    double precision, intent(in) :: d, q(n1,n2)
    double precision :: p(n1,n2)
    if (d > 0d0) then
      p = logprobt1(q,d)
    else if (d < 0d0) then
      p = logproblogis1(q)
    else
      p = logprobnorm1(q)
    end if
  end function logprobm

  function logborpm(q,d,n1,n2) result (p)
    integer, intent(in) :: n1, n2
    double precision, intent(in) :: d, q(n1,n2)
    double precision :: p(n1,n2)
    if (d > 0d0) then
      p = logprobt0(q,d)
    else if (d < 0d0) then
      p = logproblogis0(q)
    else
      p = logprobnorm0(q)
    end if
  end function logborpm

  function probinvm(p,d,n1,n2) result (q)
    integer, intent(in) :: n1, n2
    double precision, intent(in) :: d, p(n1,n2)
    double precision :: q(n1,n2)
    if (d > 0d0) then
      q = quantt(p,d)
    else if (d < 0d0) then
      q = quantlogis(p)
    else
      q = quantnorm(p)
    end if
  end function probinvm

  function logpdfm(x,d,n1,n2) result(pdf)
    integer, intent(in) :: n1,n2
    double precision, intent(in) :: x(n1,n2), d
    double precision :: pdf(n1,n2)
    if (d > 0d0) then
      pdf = logpdft(x,d)
    else if (d < 0d0) then
      pdf = logpdflogis(x)
    else
      pdf = logpdfnorm(x)
    end if
  end function logpdfm

  function logbxcxs(z,d) result (y)
    double precision, intent(in) :: d, z
    double precision y
    if (d /= 0d0) then
      y = d*z
      if (y > -1d0) then
        y = flog1p(y)/d
      else
        y = bigneg
      end if
    else
      y = z
    end if
  end function logbxcxs

  function logbxcxv(z,d,n) result (y)
    integer, intent(in) :: n
    double precision, intent(in) :: d, z(n)
    double precision y(n)
    if (d /= 0d0) then
      y = d*z
      where (y > -1d0)
        y = flog1p(y)/d
      elsewhere
        y = bigneg
      end where
    else
      y = z
    end if
  end function logbxcxv

  function logbxcxm(z,d,n1,n2) result (y)
    integer, intent(in) :: n1, n2
    double precision, intent(in) :: d, z(n1,n2)
    double precision y(n1,n2)
    if (d /= 0d0) then
      y = d*z
      where (y > -1d0)
        y = flog1p(y)/d
      elsewhere
        y = bigneg
      end where
    else
      y = z
    end if
  end function logbxcxm

  function invbxcxs(z,d,c,lext) result (y)
!! Inverse Box-Cox transformation
!! z = (y^d-1)/(d*c^(d-1)) <=> y = exp(log(1+z*d*c^(d-1))/d)
    double precision, intent(in) :: d, z, c
    logical, intent(in) :: lext
    double precision y
    if (d .eq. 0d0) then
      y = exp(z/c)
    else if (d .eq. .5d0) then
      y = 1d0 + .5d0*z/sqrt(c)
      if (y > 0d0) then
        y = y*y
      else if (lext .and. y < 0d0) then
        y = -y*y
      else
        y = 0d0
      end if
    else if (d .eq. 1d0) then
      y = 1d0 + z
    else if (d .eq. 2d0) then
      y = z*c
      y = y + y + 1d0
      if (y > 0d0) then
        y = sqrt(y)
      else if (lext .and. y < 0d0) then
        y = -sqrt(-y)
      else
        y = 0d0
      end if
    else if (d .gt. 0d0) then
      y = z*d*(c**(d-1d0)) + 1d0
      if (y > 0d0) then
        y = exp(log(y)/d)
      else if (lext .and. y < 0d0) then
        y = -exp(log(-y)/d)
      else
        y = 0d0
      end if
    end if
  end function invbxcxs

  function invbxcxv(z,d,n,c,lext) result (y)
    integer, intent(in) :: n
    double precision, intent(in) :: d, z(n), c
    logical, intent(in) :: lext
    double precision y(n)
    if (d .eq. 0d0) then
      y = exp(z/c)
    else if (d .eq. .5d0) then
      y = 1d0 + (.5d0/sqrt(c))*z
      where (y > 0d0)
        y = y*y
      elsewhere (lext .and. y < 0d0)
        y = -y*y
      elsewhere
        y = 0d0
      end where
    else if (d .eq. 1d0) then
      y = 1d0 + z
    else if (d .eq. 2d0) then
      y = z*c
      y = y + y + 1d0
      where (y > 0d0)
        y = sqrt(y)
      elsewhere (lext .and. y < 0d0)
        y = -sqrt(-y)
      elsewhere
        y = 0d0
      end where
    else if (d .gt. 0d0) then
      y = z*(d*c**(d-1d0)) + 1d0
      where (y > 0d0)
        y = exp(log(y)/d)
      elsewhere (lext .and. y < 0d0)
        y = -exp(log(-y)/d)
      elsewhere
        y = 0d0
      end where
    end if
  end function invbxcxv

  function invbxcxm(z,d,n1,n2,c,lext) result (y)
    integer, intent(in) :: n1, n2
    double precision, intent(in) :: d, z(n1,n2), c
    logical, intent(in) :: lext
    double precision y(n1,n2)
    if (d .eq. 0d0) then
      y = exp(z/c)
    else if (d .eq. .5d0) then
      y = 1d0 + (.5d0/sqrt(c))*z
      where (y > 0d0)
        y = y*y
      elsewhere (lext .and. y < 0d0)
        y = -y*y
      elsewhere
        y = 0d0
      end where
    else if (d .eq. 1d0) then
      y = 1d0 + z
    else if (d .eq. 2d0) then
      y = z*c
      y = y + y + 1d0
      where (y > 0d0)
        y = sqrt(y)
      elsewhere (lext .and. y < 0d0)
        y = -sqrt(-y)
      elsewhere
        y = 0d0
      end where
    else if (d .gt. 0d0) then
      y = z*(d*c**(d-1d0)) + 1d0
      where (y > 0d0)
        y = exp(log(y)/d)
      elsewhere (lext .and. y < 0d0)
        y = -exp(log(-y)/d)
      elsewhere
        y = 0d0
      end where
    end if
  end function invbxcxm
end module interfaces


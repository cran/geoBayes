!!! interfaces.f90 --- 
!! 
!! Author: Evangelos Evangelou
!! Created: Thu, 19 Jun, 2014 11:34 (BST)
!! Last-Updated: Fri, 23 Jan, 2015 15:07 (GMT)
!!     Update #: 11
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Contains general-purpose functions
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module interfaces
  implicit none 

!!!!!!!!!!!!!!!!!!!!!!!!!!! Correlation functions !!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    elemental double precision function matern (d, kappa)
      double precision, intent(in) :: d, kappa
    end function matern
  end interface


!!!!!!!!!!!!!!!!!!!!!!!!!!!! Auxliary functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    elemental double precision function flog1pexp (x)
      double precision, intent(in) :: x
    end function flog1pexp
  end interface

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RNG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    subroutine rngini ()
    end subroutine rngini
  end interface

  interface
    subroutine rngend ()
    end subroutine rngend
  end interface

  
!!!!!!!!!!!!!!!!!!!!!!!!!! Random sample functions !!!!!!!!!!!!!!!!!!!!!!!!!!

  interface
    double precision function randnorm ()
    end function randnorm
  end interface

  interface
    double precision function randunif ()
    end function randunif
  end interface

  interface
    pure double precision function randt (d)
      double precision, intent(in) :: d
    end function randt
  end interface

  interface
    pure double precision function randgamma (shp)
      double precision, intent(in) :: shp
    end function randgamma
  end interface

  interface
    pure double precision function randchisq (df)
      double precision, intent(in) :: df
    end function randchisq
  end interface


!!!!!!!!!!!!!!!!!!!!!!!!!!! Probability functions !!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    elemental double precision function logprobt (q, d)
      double precision, intent(in) :: q, d
    end function logprobt
  end interface
  
  interface
    elemental double precision function logborpt (q, d)
      double precision, intent(in) :: q, d
    end function logborpt
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
  
  interface
    elemental double precision function logprobnorm (q)
      double precision, intent(in) :: q
    end function logprobnorm
  end interface
  
  interface
    elemental double precision function logborpnorm (q)
      double precision, intent(in) :: q
    end function logborpnorm
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
  
  interface
    elemental double precision function logproblogis (q)
      double precision, intent(in) :: q
    end function logproblogis
  end interface
  
  interface
    elemental double precision function logborplogis (q)
      double precision, intent(in) :: q
    end function logborplogis
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
end module interfaces



! Local Variables:
! compile-command: "R CMD SHLIB"
! End:

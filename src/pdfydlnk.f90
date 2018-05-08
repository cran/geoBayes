!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Compute the derivative of log pdf(y|mu) w.r.t. mu
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pdfydlnk
contains
!!!!!!!!!!!!!! Compute the log-pdf of y for given parameter !!!!!!!!!!!!!!!
  elemental double precision function logpdfydlnk_gt (y1, y2, par)
    ! Transformed Gaussian
    ! y1 :: The average across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none 
    double precision, intent(in) :: y1, y2, par
    logpdfydlnk_gt = 2d0*y2*(par-y1)
  end function logpdfydlnk_gt

  elemental double precision function logpdfydlnk_ga (y1, y2, par)
    ! Gaussian
    ! y1 :: The total across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none 
    double precision, intent(in) :: y1, y2, par
    logpdfydlnk_ga = y1 - y2*par
  end function logpdfydlnk_ga

  elemental double precision function logpdfydlnk_bi (y1, y2, par)
    ! Binomial
    ! To get the asymmetric version, switch y1 and y2
    ! y1 :: The number of successes
    ! y2 :: The number of faillures
    ! par :: The logartithm of the probability of success
    use interfaces, only: fexpm1
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfydlnk_bi = y1 - y2/fexpm1(-par)
  end function logpdfydlnk_bi

  elemental double precision function logpdfydlnk_po (y1, y2, par)
    ! Poisson
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfydlnk_po = y1 - y2*exp(par)
  end function logpdfydlnk_po

  elemental double precision function logpdfydlnk_gm (y1, y2, par)
    ! Gamma
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfydlnk_gm = y1*exp(-par) - y2
  end function logpdfydlnk_gm
end module pdfydlnk

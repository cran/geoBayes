module modelfcns_pdfy
!! logpdfy is the log-pdf of y for given parameter value
!! logpdfydlnk is the derivative of the above w.r.t the parameter
!! logdffy is the difference between two such logpdfs at two different
!! parameter values
  implicit none 
  double precision, parameter :: bigpos = huge(1d0), bigneg = -bigpos, &
     smallpos = epsilon(1d0), smallneg = -smallpos
  private bigpos, bigneg, smallpos, smallneg

contains

!!!!!!!!!!!!!!!!!!!!!!!!!! Transformed Gaussian !!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure elemental double precision function logpdfy_gt (y1, y2, par)
    ! Transformed Gaussian
    ! y1 :: The average across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    double precision d
    d = y1 - par
    logpdfy_gt = y2*d*d
  end function logpdfy_gt

  pure elemental double precision function logdffy_gt (y1, y2, p1, p2)
    ! Transformed Gaussian
    ! y1 :: The average across all realisations
    ! y2 :: The number of replications
    ! p1, p2 :: The mean of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, p1, p2
    double precision d1, d2
    d1 = y1 - p1
    d2 = y1 - p2
    logdffy_gt = y2*(d1*d1 - d2*d2)
  end function logdffy_gt

  pure elemental double precision function logpdfydlnk_gt (y1, y2, par)
    ! Transformed Gaussian
    ! y1 :: The average across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfydlnk_gt = 2d0*y2*(par-y1)
  end function logpdfydlnk_gt

  pure elemental double precision function logpdfyhlnk_gt (y1, y2, par)
    ! Transformed Gaussian
    ! y1 :: The average across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfyhlnk_gt = 2d0*y2
  end function logpdfyhlnk_gt

  pure elemental double precision function logpdfy3lnk_gt (y1, y2, par)
    ! Gaussian
    ! y1 :: The total across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfy3lnk_gt = 0
  end function logpdfy3lnk_gt

  pure elemental double precision function mustart_gt (y1,y2)
    implicit none
    double precision, intent(in) :: y1, y2
    mustart_gt = y1
  end function mustart_gt

  pure elemental double precision function fcncum_gt (mu) result (w)
    ! Cumulant fcn as a fcn of mu
    implicit none
    double precision, intent(in) :: mu
    w = .5d0*mu*mu
  end function fcncum_gt

  pure elemental double precision function fcncumd2_gt (mu) result (w)
    ! Cumulant fcn as a fcn of mu, 2nd derivative
    implicit none
    double precision, intent(in) :: mu
    w = 1d0
  end function fcncumd2_gt

  pure elemental double precision function fcncumd3_gt (mu) result (w)
    ! Cumulant fcn as a fcn of mu, 3rd derivative
    implicit none
    double precision, intent(in) :: mu
    w = 0d0
  end function fcncumd3_gt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gaussian !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental double precision function logpdfy_ga (y1, y2, par)
    ! Gaussian
    ! y1 :: The total across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfy_ga = y1*par - .5d0*y2*par*par
  end function logpdfy_ga

  pure elemental double precision function logdffy_ga (y1, y2, p1, p2)
    ! Gaussian
    ! y1 :: The total across all realisations
    ! y2 :: The number of replications
    ! p1, p2 :: The mean of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, p1, p2
    logdffy_ga = y1*(p1 - p2) - .5d0*y2*(p1*p1 - p2*p2)
  end function logdffy_ga

  pure elemental double precision function logpdfydlnk_ga (y1, y2, par)
    ! Gaussian
    ! y1 :: The total across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfydlnk_ga = y1 - y2*par
  end function logpdfydlnk_ga

  pure elemental double precision function logpdfyhlnk_ga (y1, y2, par)
    ! Gaussian
    ! y1 :: The total across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfyhlnk_ga = -y2
  end function logpdfyhlnk_ga

  pure elemental double precision function logpdfy3lnk_ga (y1, y2, par)
    ! Gaussian
    ! y1 :: The total across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfy3lnk_ga = 0
  end function logpdfy3lnk_ga

  pure elemental double precision function mustart_ga (y1,y2)
    implicit none
    double precision, intent(in) :: y1, y2
    mustart_ga = y1/y2
  end function mustart_ga

  pure elemental double precision function fcncum_ga (mu) result (w)
    ! Cumulant fcn as a fcn of mu
    implicit none
    double precision, intent(in) :: mu
    w = .5d0*mu*mu
  end function fcncum_ga

  pure elemental double precision function fcncumd2_ga (mu) result (w)
    ! Cumulant fcn as a fcn of mu, 2nd derivative
    implicit none
    double precision, intent(in) :: mu
    w = 1d0
  end function fcncumd2_ga

  pure elemental double precision function fcncumd3_ga (mu) result (w)
    ! Cumulant fcn as a fcn of mu, 3rd derivative
    implicit none
    double precision, intent(in) :: mu
    w = 0d0
  end function fcncumd3_ga


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental double precision function logpdfy_bi (y1, y2, par)
    ! Binomial
    ! To get the asymmetric version, switch y1 and y2
    ! y1 :: The number of successes
    ! y2 :: The number of faillures
    ! par :: The logartithm of the probability of success
    use interfaces, only: flog1mexp
    implicit none
    double precision, intent(in) :: y1, y2, par
    if (y1 .eq. 0d0) then
      logpdfy_bi = y2*flog1mexp(par)
    else if (y2 .eq. 0d0) then
      logpdfy_bi = y1*par
    else 
      logpdfy_bi = y1*par + y2*flog1mexp(par)
    end if
  end function logpdfy_bi

  pure elemental double precision function logdffy_bi (y1, y2, p1, p2)
    ! Binomial
    ! To get the asymmetric version, switch y1 and y2
    ! y1 :: The number of successes
    ! y2 :: The number of faillures
    ! p1, p2 :: The logartithm of the probability of success
    use interfaces, only: flog1mexp
    implicit none
    double precision, intent(in) :: y1, y2, p1, p2
    if (y1 .eq. 0d0) then
      logdffy_bi = y2*(flog1mexp(p1) - flog1mexp(p2))
    else if (y2 .eq. 0d0) then
      logdffy_bi = y1*(p1 - p2)
    else
      logdffy_bi = y1*(p1 - p2) + y2*(flog1mexp(p1) - flog1mexp(p2))
    end if
  end function logdffy_bi

  pure elemental double precision function logpdfydlnk_bi (y1, y2, par)
    ! Binomial
    ! To get the asymmetric version, switch y1 and y2
    ! y1 :: The number of successes
    ! y2 :: The number of faillures
    ! par :: The logartithm of the probability of success
    use interfaces, only: fexpm1
    implicit none
    double precision, intent(in) :: y1, y2, par
    if (y2 .eq. 0d0) then
      logpdfydlnk_bi = y1
    else
      logpdfydlnk_bi = y1 - y2/fexpm1(-par)
    end if
  end function logpdfydlnk_bi

  pure elemental double precision function logpdfyhlnk_bi (y1, y2, par)
    ! Binomial
    ! To get the asymmetric version, switch y1 and y2
    ! y1 :: The number of successes
    ! y2 :: The number of faillures
    ! par :: The logartithm of the probability of success
    use interfaces, only: fexpm1
    implicit none
    double precision, intent(in) :: y1, y2, par
    if (y2 .eq. 0d0) then
      logpdfyhlnk_bi = 0d0
    else
      logpdfyhlnk_bi = 1d0/fexpm1(-par)
      logpdfyhlnk_bi = -y2*logpdfyhlnk_bi * (1d0+logpdfyhlnk_bi)
    end if
  end function logpdfyhlnk_bi

  pure elemental double precision function logpdfy3lnk_bi (y1, y2, par)
    ! Binomial
    ! To get the asymmetric version, switch y1 and y2
    ! y1 :: The number of successes
    ! y2 :: The number of faillures
    ! par :: The logartithm of the probability of success
    use interfaces, only: fexpm1
    implicit none
    double precision, intent(in) :: y1, y2, par
    if (y2 .eq. 0d0) then
      logpdfy3lnk_bi = 0d0
    else
      logpdfy3lnk_bi = 1d0/fexpm1(-par)
      logpdfy3lnk_bi = -y2*logpdfy3lnk_bi * (1d0+logpdfy3lnk_bi) &
         * (1d0+logpdfy3lnk_bi+logpdfy3lnk_bi)
    end if
  end function logpdfy3lnk_bi

  pure elemental double precision function mustart_bi (y1,y2)
    implicit none
    double precision, intent(in) :: y1, y2
    double precision t1, t2
    t1 = y1 + .5d0
    t2 = y1 + y2 + 1d0
    mustart_bi = log(t1) - log(t2)
  end function mustart_bi

  pure elemental double precision function fcncum_bi (mu) result (w)
    ! Cumulant fcn as a fcn of mu
    implicit none
    double precision, intent(in) :: mu
    w = -log(1-mu)
  end function fcncum_bi

  pure elemental double precision function fcncumd2_bi (mu) result (w)
    ! Cumulant fcn as a fcn of mu, 2nd derivative
    implicit none
    double precision, intent(in) :: mu
    w = mu*(1d0-mu)
  end function fcncumd2_bi

  pure elemental double precision function fcncumd3_bi (mu) result (w)
    ! Cumulant fcn as a fcn of mu, 3rd derivative
    implicit none
    double precision, intent(in) :: mu
    w = mu*(1d0-mu)*(1d0-mu-mu)
  end function fcncumd3_bi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Poisson !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental double precision function logpdfy_po (y1, y2, par)
    ! Poisson
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    if (y1 .eq. 0d0) then
      logpdfy_po = -y2*exp(par)
    else
      logpdfy_po = y1*par - y2*exp(par)
    end if
  end function logpdfy_po

  pure elemental double precision function logdffy_po (y1, y2, p1, p2)
    ! Poisson
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! p1, p2 :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, p1, p2
    if (y1 .eq. 0d0) then
      logdffy_po = -y2*(exp(p1) - exp(p2))
    else
      logdffy_po = y1*(p1 - p2) - y2*(exp(p1) - exp(p2))
    end if
  end function logdffy_po

  pure elemental double precision function logpdfydlnk_po (y1, y2, par)
    ! Poisson
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfydlnk_po = y1 - y2*exp(par)
  end function logpdfydlnk_po

  pure elemental double precision function logpdfyhlnk_po (y1, y2, par)
    ! Poisson
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfyhlnk_po = -y2*exp(par)
  end function logpdfyhlnk_po

  pure elemental double precision function logpdfy3lnk_po (y1, y2, par)
    ! Poisson
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfy3lnk_po = -y2*exp(par)
  end function logpdfy3lnk_po

  pure elemental double precision function mustart_po (y1,y2)
    implicit none
    double precision, intent(in) :: y1, y2
    double precision t1, t2
    t1 = y1 + .5d0
    t2 = y2 + 1d0
    mustart_po = log(t1) - log(t2)
  end function mustart_po

  pure elemental double precision function fcncum_po (mu) result (w)
    ! Cumulant fcn as a fcn of mu
    implicit none
    double precision, intent(in) :: mu
    w = mu
  end function fcncum_po

  pure elemental double precision function fcncumd2_po (mu) result (w)
    ! Cumulant fcn as a fcn of mu, 2nd derivative
    implicit none
    double precision, intent(in) :: mu
    w = mu
  end function fcncumd2_po

  pure elemental double precision function fcncumd3_po (mu) result (w)
    ! Cumulant fcn as a fcn of mu, 3rd derivative
    implicit none
    double precision, intent(in) :: mu
    w = mu
  end function fcncumd3_po


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gamma !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure elemental double precision function logpdfy_gm (y1, y2, par)
    ! Gamma
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfy_gm = -y1*exp(-par) - y2*par
  end function logpdfy_gm

  pure elemental double precision function logdffy_gm (y1, y2, p1, p2)
    ! Gamma
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! p1, p2 :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, p1, p2
    logdffy_gm = -y1*(exp(-p1) - exp(-p2)) - y2*(p1 - p2)
  end function logdffy_gm

  pure elemental double precision function logpdfydlnk_gm (y1, y2, par)
    ! Gamma
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfydlnk_gm = y1*exp(-par) - y2
  end function logpdfydlnk_gm

  pure elemental double precision function logpdfyhlnk_gm (y1, y2, par)
    ! Gamma
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfyhlnk_gm = -y1*exp(-par)
  end function logpdfyhlnk_gm

  pure elemental double precision function logpdfy3lnk_gm (y1, y2, par)
    ! Gamma
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfy3lnk_gm = y1*exp(-par)
  end function logpdfy3lnk_gm

  pure elemental double precision function mustart_gm (y1,y2)
    implicit none
    double precision, intent(in) :: y1, y2
    mustart_gm = log(y1) - log(y2)
  end function mustart_gm

  pure elemental double precision function fcncum_gm (mu) result (w)
    ! Cumulant fcn as a fcn of mu
    implicit none
    double precision, intent(in) :: mu
    w = log(mu)
  end function fcncum_gm

  pure elemental double precision function fcncumd2_gm (mu) result (w)
    ! Cumulant fcn as a fcn of mu, 2nd derivative
    implicit none
    double precision, intent(in) :: mu
    w = mu*mu
  end function fcncumd2_gm

  pure elemental double precision function fcncumd3_gm (mu) result (w)
    ! Cumulant fcn as a fcn of mu, 3rd derivative
    implicit none
    double precision, intent(in) :: mu
    w = mu*mu*mu
    w = w + w
  end function fcncumd3_gm

end module modelfcns_pdfy

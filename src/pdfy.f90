!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Computes the conditional log-pdf of y
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pdfy
contains
!!!!!!!!!!!!!! Compute the log-pdf of y for given parameter !!!!!!!!!!!!!!!
  elemental double precision function logpdfy_gt (y1, y2, par)
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

  elemental double precision function logpdfy_ga (y1, y2, par)
    ! Gaussian
    ! y1 :: The total across all realisations
    ! y2 :: The number of replications
    ! par :: The mean of a single realisation
    implicit none 
    double precision, intent(in) :: y1, y2, par
    logpdfy_ga = y1*par - .5d0*y2*par*par
  end function logpdfy_ga

  elemental double precision function logpdfy_bi (y1, y2, par)
    ! Binomial
    ! To get the asymmetric version, switch y1 and y2
    ! y1 :: The number of successes
    ! y2 :: The number of faillures
    ! par :: The logartithm of the probability of success
    use interfaces, only: flog1mexp
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfy_bi = y1*par + y2*flog1mexp(par)
  end function logpdfy_bi

  elemental double precision function logpdfy_po (y1, y2, par)
    ! Poisson
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfy_po = y1*par - y2*exp(par)
  end function logpdfy_po

  elemental double precision function logpdfy_gm (y1, y2, par)
    ! Gamma
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! par :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, par
    logpdfy_gm = -y1*exp(-par) - y2*par
  end function logpdfy_gm

  
!!!!!!!!! Compute the difference log pdf of y wrt two parameters !!!!!!!!!!
  elemental double precision function logdffy_gt (y1, y2, p1, p2)
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

  elemental double precision function logdffy_ga (y1, y2, p1, p2)
    ! Gaussian
    ! y1 :: The total across all realisations
    ! y2 :: The number of replications
    ! p1, p2 :: The mean of a single realisation
    implicit none 
    double precision, intent(in) :: y1, y2, p1, p2
    logdffy_ga = y1*(p1 - p2) - .5d0*y2*(p1*p1 - p2*p2)
  end function logdffy_ga

  elemental double precision function logdffy_bi (y1, y2, p1, p2)
    ! Binomial
    ! To get the asymmetric version, switch y1 and y2
    ! y1 :: The number of successes
    ! y2 :: The number of faillures
    ! p1, p2 :: The logartithm of the probability of success
    use interfaces, only: flog1mexp
    implicit none
    double precision, intent(in) :: y1, y2, p1, p2
    logdffy_bi = y1*(p1 - p2) + y2*(flog1mexp(p1) - flog1mexp(p2))
  end function logdffy_bi

  elemental double precision function logdffy_po (y1, y2, p1, p2)
    ! Poisson
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! p1, p2 :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, p1, p2
    logdffy_po = y1*(p1 - p2) - y2*(exp(p1) - exp(p2))
  end function logdffy_po

  elemental double precision function logdffy_gm (y1, y2, p1, p2)
    ! Gamma
    ! y1 :: The total number of observations
    ! y2 :: The number of realisations
    ! p1, p2 :: The log(mean) of a single realisation
    implicit none
    double precision, intent(in) :: y1, y2, p1, p2
    logdffy_gm = -y1*(exp(-p1) - exp(-p2)) - y2*(p1 - p2)
  end function logdffy_gm
end module pdfy




  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Computes the log-pdf of z
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pdfz
contains
  function logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, modeldfh, xi(n)
    double precision logpdfz
    double precision Upsz(n), zUz, zmxi(n)
    if (lmxi) then
      zmxi = z - xi
      call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
      zUz = dot_product(zmxi,Upsz) + ssqdfsc
    else
      call dsymv ('u',n,1d0,Ups,n,z,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
      zUz = dot_product(z,Upsz) + ssqdfsc
    end if
    logpdfz = ldh_Ups - modeldfh*log(zUz)
  end function logpdfz
end module pdfz




  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Computes the conditional log-pdf of y|z
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module condyz
contains
  function condyz_gt (n, y1, y2, z, nu, tsqdfsc, respdfh)
    use pdfy, only: logpdfy_gt
    use linkfcn, only: invlink_ga
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsqdfsc, respdfh
    double precision condyz_gt
    integer i
    double precision mu, lfy
    lfy = tsqdfsc
    do i = 1, n
      mu = invlink_ga(z(i),nu)
      lfy = lfy + logpdfy_gt(y1(i), y2(i), mu)
    end do
    condyz_gt = -respdfh*log(lfy)
  end function condyz_gt

  function condyz_ga (n, y1, y2, z, nu, tsq)
    use pdfy, only: logpdfy_ga
    use linkfcn, only: invlink_ga
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision condyz_ga
    integer i
    double precision mu, lfy
    lfy = 0d0
    do i = 1, n
      mu = invlink_ga(z(i),nu)
      lfy = lfy + logpdfy_ga(y1(i), y2(i), mu)
    end do
    condyz_ga = lfy/tsq
  end function condyz_ga

  function condyz_bi (n, y1, y2, z, nu, tsq)
    use pdfy, only: logpdfy_bi
    use linkfcn, only: invlink_bi
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision condyz_bi
    integer i
    double precision mu, lfy
    lfy = 0d0
    do i = 1, n
      mu = invlink_bi(z(i),nu)
      lfy = lfy + logpdfy_bi(y1(i), y2(i), mu)
    end do
    condyz_bi = lfy/tsq
  end function condyz_bi

  function condyz_po (n, y1, y2, z, nu, tsq)
    use pdfy, only: logpdfy_po
    use linkfcn, only: invlink_po
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision condyz_po
    integer i
    double precision mu, lfy
    lfy = 0d0
    do i = 1, n
      mu = invlink_po(z(i),nu)
      lfy = lfy + logpdfy_po(y1(i), y2(i), mu)
    end do
    condyz_po = lfy/tsq
  end function condyz_po

  function condyz_gm (n, y1, y2, z, nu, tsq)
    use pdfy, only: logpdfy_gm
    use linkfcn, only: invlink_gm
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision condyz_gm
    integer i
    double precision mu, lfy
    lfy = 0d0
    do i = 1, n
      mu = invlink_gm(z(i),nu)
      lfy = lfy + logpdfy_gm(y1(i), y2(i), mu)
    end do
    condyz_gm = lfy/tsq
  end function condyz_gm

  function condyz_ba (n, y1, y2, z, nu, tsq)
    use pdfy, only: logpdfy_bi
    use linkfcn, only: invlink_ba
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision condyz_ba
    integer i
    double precision mu, lfy
    lfy = 0d0
    do i = 1, n
      mu = invlink_ba(z(i),nu)
      lfy = lfy + logpdfy_bi(y2(i), y1(i), mu)
    end do
    condyz_ba = lfy/tsq
  end function condyz_ba

  function condyz_bd (n, y1, y2, z, nu, tsq)
    use pdfy, only: logpdfy_bi
    use linkfcn, only: invlink_bd
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision condyz_bd
    integer i
    double precision mu, lfy
    lfy = 0d0
    do i = 1, n
      mu = invlink_bd(z(i),nu)
      lfy = lfy + logpdfy_bi(y2(i), y1(i), mu)
    end do
    condyz_bd = lfy/tsq
  end function condyz_bd
end module condyz




  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Computes the joint log-pdf of (y,z)
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module jointyz
contains
  function jointyz_gt(n, z, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
    use condyz, only: condyz_gt
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, tsqdfsc, modeldfh, respdfh, nu, xi(n)
    double precision jointyz_gt
    double precision lfz, lfy
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condyz_gt(n, y, l, z, nu, tsqdfsc, respdfh)
    jointyz_gt = lfz + lfy
  end function jointyz_gt

  function jointyz_ga (n, z, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condyz, only: condyz_ga
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, tsq, nu, modeldfh, xi(n)
    double precision jointyz_ga
    double precision lfz, lfy
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condyz_ga(n, y, l, z, nu, tsq)
    jointyz_ga = lfz + lfy
  end function jointyz_ga

  function jointyz_bi (n, z, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condyz, only: condyz_bi
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, tsq, nu, modeldfh, xi(n)
    double precision jointyz_bi
    double precision lfz, lfy
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condyz_bi(n, y, l, z, nu, tsq)
    jointyz_bi = lfz + lfy
  end function jointyz_bi

  function jointyz_po (n, z, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condyz, only: condyz_po
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, tsq, nu, modeldfh, xi(n)
    double precision jointyz_po
    double precision lfz, lfy
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condyz_po(n, y, l, z, nu, tsq)
    jointyz_po = lfz + lfy
  end function jointyz_po

  function jointyz_gm (n, z, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condyz, only: condyz_gm
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, tsq, nu, modeldfh, xi(n)
    double precision jointyz_gm
    double precision lfz, lfy
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condyz_gm(n, y, l, z, nu, tsq)
    jointyz_gm = lfz + lfy
  end function jointyz_gm

  function jointyz_ba (n, z, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condyz, only: condyz_ba
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, tsq, nu, modeldfh, xi(n)
    double precision jointyz_ba
    double precision lfz, lfy
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condyz_ba(n, y, l, z, nu, tsq)
    jointyz_ba = lfz + lfy
  end function jointyz_ba

  function jointyz_bd (n, z, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condyz, only: condyz_bd
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, tsq, nu, modeldfh, xi(n)
    double precision jointyz_bd
    double precision lfz, lfy
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condyz_bd(n, y, l, z, nu, tsq)
    jointyz_bd = lfz + lfy
  end function jointyz_bd
end module jointyz






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Computes the log-pdf of mu
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pdfmu
  implicit none
  double precision, parameter :: bigneg = -huge(1d0)
  private bigneg
contains
  function logpdfmu_ga (n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    use linkfcn, only: flink_ga
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, modeldfh
    double precision logpdfmu_ga
    integer i
    double precision z(n), logjac, lfz
    ! Linear predictor
    do i = 1, n
      z(i) = flink_ga(mu(i), nu)
    end do
    ! Jacobian
    if (nu .gt. 0d0) then
      logjac = 0d0
      do i = 1, n
        logjac = logjac + log(abs(mu(i)))
      end do
      logjac = (nu - 1d0)*logjac
    else if (all(mu .gt. 0d0)) then
      logjac = 0d0
      do i = 1, n
        logjac = logjac + log(mu(i))
      end do
      logjac = (nu - 1d0)*logjac
    else
      logjac = bigneg
    end if
    ! log-likelihood for z
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    ! Put all together
    logpdfmu_ga = lfz + logjac
  end function logpdfmu_ga

  function logpdfmu_bi (n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    use linkfcn, only: flink_bi
    use interfaces, only: flog1pexp, flog1p ! logpdft, logpdfnorm, logpdflogis
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, modeldfh
    double precision logpdfmu_bi
    integer i
    double precision z(n), logjac, lfz, tmp
    double precision, parameter :: lgits = 0.6458d0
    ! Linear predictor
    do i = 1, n
      z(i) = flink_bi(mu(i), nu)
    end do
    ! Jacobian
    logjac = 0d0
    if (nu .gt. 0d0) then ! pdf t
      do i = 1, n
        tmp = z(i)*z(i)/nu
        logjac = logjac + flog1p(tmp)
      end do
      logjac = .5d0*logjac*(nu + 1d0)
    else if (nu .lt. 0d0) then ! pdf logistic 
      do i = 1, n
        tmp = -z(i)/lgits
        logjac = logjac - tmp + 2d0*flog1pexp(tmp)
        ! - logpdflogis(z(i))
      end do
    else ! pdf normal
      do i = 1, n
        logjac = logjac + .5*z(i)*z(i) ! - logpdfnorm(z(i))
      end do
    end if
    ! log-likelihood for z
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    ! Put all together
    logpdfmu_bi = lfz + logjac
  end function logpdfmu_bi

  function logpdfmu_po (n, tht, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, &
     modeldfh)
    ! tht is log(mean)
    use linkfcn, only: flink_po
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: tht(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, modeldfh
    double precision logpdfmu_po
    integer i
    double precision z(n), logjac, lfz
    ! Linear predictor
    do i = 1, n
      z(i) = flink_po(tht(i), nu)
    end do
    ! Jacobian
    logjac = 0d0
    if (nu .ne. 1d0) then
      do i = 1, n
        logjac = logjac + tht(i)
      end do
      logjac = (nu - 1d0)*logjac
    end if
    ! log-likelihood for z
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    ! Put all together
    logpdfmu_po = lfz + logjac
  end function logpdfmu_po

  function logpdfmu_gm (n, tht, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, &
     modeldfh)
    ! tht is log(mean)
    use linkfcn, only: flink_gm
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: tht(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, modeldfh
    double precision logpdfmu_gm
    integer i
    double precision z(n), logjac, lfz
    ! Linear predictor
    do i = 1, n
      z(i) = flink_gm(tht(i), nu)
    end do
    ! Jacobian
    logjac = 0d0
    if (nu .ne. 1d0) then
      do i = 1, n
        logjac = logjac + tht(i)
      end do
      logjac = (nu - 1d0)*logjac
    end if
    ! log-likelihood for z
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    ! Put all together
    logpdfmu_gm = lfz + logjac
  end function logpdfmu_gm

  function logpdfmu_ba (n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    use linkfcn, only: flink_ba
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, modeldfh
    double precision logpdfmu_ba
    integer i
    double precision z(n), logjac, lfz, logjac1, logjac2
    ! Linear predictor
    do i = 1, n
      z(i) = flink_ba(mu(i), nu)
    end do
    ! Jacobian
    logjac1 = 0d0
    logjac2 = 0d0
    do i = 1, n
      logjac1 = logjac1 + log(-mu(i))
      logjac2 = logjac2 - mu(i)
    end do
    logjac = (nu - 1d0)*logjac1 + logjac2
    ! log-likelihood for z
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    ! Put all together
    logpdfmu_ba = lfz + logjac
  end function logpdfmu_ba

  function logpdfmu_bd (n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    use linkfcn, only: flink_bd
    use pdfz
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, modeldfh
    double precision logpdfmu_bd
    integer i
    double precision z(n), logjac, lfz, logjac1, logjac2
    ! Linear predictor
    do i = 1, n
      z(i) = flink_bd(mu(i), nu)
    end do
    ! Jacobian
    logjac1 = 0d0
    logjac2 = 0d0
    do i = 1, n
      logjac1 = logjac1 + log(-mu(i))
      logjac2 = logjac2 - mu(i)
    end do
    logjac = (nu - 1d0)*logjac1 + logjac2
    ! log-likelihood for z
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    ! Put all together
    logpdfmu_bd = lfz + logjac
  end function logpdfmu_bd
end module pdfmu







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Computes the conditional log-pdf of y|mu
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module condymu
contains
  function condymu_gt (n, y1, y2, mu, tsqdfsc, respdfh)
    use pdfy, only: logpdfy_gt
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqdfsc, &
       respdfh
    double precision condymu_gt
    integer i
    double precision lfy
    lfy = tsqdfsc
    do i = 1, n
      lfy = lfy + logpdfy_gt(y1(i), y2(i), mu(i))
    end do
    condymu_gt = -respdfh*log(lfy)
  end function condymu_gt

  function condymu_ga (n, y1, y2, mu, tsq)
    use pdfy, only: logpdfy_ga
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsq
    double precision condymu_ga
    integer i
    double precision lfy
    lfy = 0d0
    do i = 1, n
      lfy = lfy + logpdfy_ga(y1(i), y2(i), mu(i))
    end do
    condymu_ga = lfy/tsq
  end function condymu_ga

  function condymu_bi (n, y1, y2, mu, tsq)
    use pdfy, only: logpdfy_bi
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsq
    double precision condymu_bi
    integer i
    double precision lfy
    lfy = 0d0
    do i = 1, n
      lfy = lfy + logpdfy_bi(y1(i), y2(i), mu(i))
    end do
    condymu_bi = lfy/tsq
  end function condymu_bi

  function condymu_po (n, y1, y2, mu, tsq)
    use pdfy, only: logpdfy_po
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsq
    double precision condymu_po
    integer i
    double precision lfy
    lfy = 0d0
    do i = 1, n
      lfy = lfy + logpdfy_po(y1(i), y2(i), mu(i))
    end do
    condymu_po = lfy/tsq
  end function condymu_po

  function condymu_gm (n, y1, y2, mu, tsq)
    use pdfy, only: logpdfy_gm
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsq
    double precision condymu_gm
    integer i
    double precision lfy
    lfy = 0d0
    do i = 1, n
      lfy = lfy + logpdfy_gm(y1(i), y2(i), mu(i))
    end do
    condymu_gm = lfy/tsq
  end function condymu_gm

  function condymu_ba (n, y1, y2, mu, tsq)
    use pdfy, only: logpdfy_bi
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsq
    double precision condymu_ba
    integer i
    double precision lfy
    lfy = 0d0
    do i = 1, n
      lfy = lfy + logpdfy_bi(y2(i), y1(i), mu(i))
    end do
    condymu_ba = lfy/tsq
  end function condymu_ba

  function condymu_bd (n, y1, y2, mu, tsq)
    use pdfy, only: logpdfy_bi
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsq
    double precision condymu_bd
    integer i
    double precision lfy
    lfy = 0d0
    do i = 1, n
      lfy = lfy + logpdfy_bi(y2(i), y1(i), mu(i))
    end do
    condymu_bd = lfy/tsq
  end function condymu_bd
end module condymu







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Computes the joint log-pdf of (y,mu)
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module jointymu
  implicit none
  double precision, parameter :: bigneg = -huge(1d0)
  private bigneg
contains
  function jointymu_gt (n, mu, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
    use condymu, only: condymu_gt
    use pdfmu, only: logpdfmu_ga
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, tsqdfsc, modeldfh, respdfh
    double precision jointymu_gt
    double precision lfmu, lfy
    lfmu = logpdfmu_ga(n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condymu_gt(n, y, l, mu, tsqdfsc, respdfh)
    jointymu_gt = lfmu + lfy
  end function jointymu_gt

  function jointymu_ga (n, mu, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condymu, only: condymu_ga
    use pdfmu, only: logpdfmu_ga
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, tsq, modeldfh
    double precision jointymu_ga
    double precision lfmu, lfy
    lfmu = logpdfmu_ga(n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condymu_ga(n, y, l, mu, tsq)
    jointymu_ga = lfmu + lfy
  end function jointymu_ga

  function jointymu_bi (n, mu, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condymu, only: condymu_bi
    use pdfmu, only: logpdfmu_bi
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, tsq, modeldfh
    double precision jointymu_bi
    double precision lfmu, lfy
    lfmu = logpdfmu_bi(n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condymu_bi(n, y, l, mu, tsq)
    jointymu_bi = lfmu + lfy
  end function jointymu_bi

  function jointymu_po (n, mu, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condymu, only: condymu_po
    use pdfmu, only: logpdfmu_po
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, tsq, modeldfh
    double precision jointymu_po
    double precision lfmu, lfy
    lfmu = logpdfmu_po(n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condymu_po(n, y, l, mu, tsq)
    jointymu_po = lfmu + lfy
  end function jointymu_po

  function jointymu_gm (n, mu, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condymu, only: condymu_gm
    use pdfmu, only: logpdfmu_gm
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, tsq, modeldfh
    double precision jointymu_gm
    double precision lfmu, lfy
    lfmu = logpdfmu_gm(n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condymu_gm(n, y, l, mu, tsq)
    jointymu_gm = lfmu + lfy
  end function jointymu_gm

  function jointymu_ba (n, mu, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condymu, only: condymu_ba
    use pdfmu, only: logpdfmu_ba
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, tsq, modeldfh
    double precision jointymu_ba
    double precision lfmu, lfy
    lfmu = logpdfmu_ba(n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condymu_ba(n, y, l, mu, tsq)
    jointymu_ba = lfmu + lfy
  end function jointymu_ba

  function jointymu_bd (n, mu, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condymu, only: condymu_bd
    use pdfmu, only: logpdfmu_bd
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, tsq, modeldfh
    double precision jointymu_bd
    double precision lfmu, lfy
    lfmu = logpdfmu_bd(n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, &
       modeldfh) 
    lfy = condymu_bd(n, y, l, mu, tsq)
    jointymu_bd = lfmu + lfy
  end function jointymu_bd
end module jointymu

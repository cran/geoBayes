module cor_fcns
  double precision, parameter :: bigneg = -huge(1d0)
  private matern, logbesselk_dk, logbesselk_hk, bigneg

!!!!!!!!!!!!!!!!!!!!!!!!!! Correlation functions !!!!!!!!!!!!!!!!!!!!!!!!!!
contains

!!! Matern
  elemental double precision function matern (h,k) result (c)
    use interfaces, only: fgamma, fbesselk
    implicit none
    double precision, intent(in) :: h, k
    if (h .eq. 0d0) then
      c = 1d0
    else
      c = 2d0/fgamma(k) * (.5d0*h)**k * fbesselk(h,k)
    end if
  end function matern

  elemental double precision function logbesselk_dk (h,k) result (c)
    ! Numerical derivatives. See
    ! http://www.math.ubc.ca/~jfeng/CHBE553/Example7/Formulae.pdf
    use interfaces, only: fbesselkexp
    implicit none
    double precision, intent(in) :: h, k
    double precision eps
    double precision kk
    eps = sqrt(epsilon(1d0))
    kk = abs(k)
    if (kk .eq. 0d0) then
      c = 0d0
    else if (kk .gt. eps) then
      c = .5d0*(log(fbesselkexp(h,kk+eps)) - log(fbesselkexp(h,kk-eps)))
    else
      c = .5d0*(-log(fbesselkexp(h,kk+eps+eps)) &
         + 4d0*log(fbesselkexp(h,kk+eps)) - 3d0*log(fbesselkexp(h,kk)))
    end if
    c = c/eps
    if (k .lt. 0d0) c = -c
  end function logbesselk_dk

  elemental double precision function logbesselk_hk (h,k) result (c)
    ! Numerical derivatives. See
    ! http://www.math.ubc.ca/~jfeng/CHBE553/Example7/Formulae.pdf
    use interfaces, only: fbesselkexp
    implicit none
    double precision, intent(in) :: h, k
    double precision eps
    double precision kk
    eps = sqrt(epsilon(1d0))
    kk = abs(k)
    if (kk .eq. 0d0) then
      c = 0d0
    else if (kk .gt. eps) then
      c = log(fbesselkexp(h,kk+eps)) + log(fbesselkexp(h,kk-eps)) &
         - 2d0*log(fbesselkexp(h,kk))
    else
      c = 4d0*log(fbesselkexp(h,kk+eps+eps)) + 2d0*log(fbesselkexp(h,kk)) &
         - 5d0*log(fbesselkexp(h,kk+eps)) - log(fbesselkexp(h,kk+eps*3d0))
    end if
    c = c/(eps*eps)
  end function logbesselk_hk

  pure double precision function cor_matern (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 1d0
    else if (h .gt. 0d0) then
      if (kappa .eq. .5d0) then
        c = exp(-h)
      else if (kappa .eq. 1.5d0) then
        c = exp(-h)*(1d0 + h)
      else if (kappa .eq. 2.5d0) then
        c = exp(-h)*(1d0 + h + h*h/3d0)
      else if (kappa .gt. 0d0) then
        c = matern(h,kappa)
      end if
    else
      c = bigneg
    end if
  end function cor_matern

  pure double precision function cor_dh_matern (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      if (kappa .eq. .5d0) then
        c = -exp(-h)
      else if (kappa .eq. 1.5d0) then
        c = -h*exp(-h)
      else if (kappa .eq. 2.5d0) then
        c = -exp(-h)*(h + h*h)/3d0
      else if (kappa .gt. 0d0) then
        c = matern_dh(h,kappa)
      end if
    else
      c = bigneg
    end if
  contains
    elemental double precision function matern_dh (h,k) result (c)
      use interfaces, only: fgamma, fbesselk
      implicit none
      double precision, intent(in) :: h, k
      if (h .eq. 0d0) then
        c = 0d0
      else
        c = -2d0/fgamma(k) * (.5d0*h)**k * fbesselk(h,k-1d0)
      end if
    end function matern_dh
  end function cor_dh_matern

  pure double precision function cor_hh_matern (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      if (kappa .eq. .5d0) then
        c = exp(-h)
      else if (kappa .eq. 1.5d0) then
        c = (h-1d0)*exp(-h)
      else if (kappa .eq. 2.5d0) then
        c = exp(-h)*(-1d0 - h + h*h)/3d0
      else if (kappa .gt. 0d0) then
        c = matern_hh(h,kappa)
      end if
    else
      c = bigneg
    end if
  contains
    elemental double precision function matern_hh (h,k) result (c)
      use interfaces, only: fgamma, fbesselk
      implicit none
      double precision, intent(in) :: h, k
      if (h .eq. 0d0) then
        c = 0d0
      else
        c = 1d0/fgamma(k) * (.5d0*h)**(k-1) * &
           (h*fbesselk(h,k-2d0) - fbesselk(h,k-1d0))
      end if
    end function matern_hh
  end function cor_hh_matern

  pure double precision function cor_dk_matern (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      c = matern_dk(h,kappa)
    end if
  contains
    elemental double precision function matern_dk (h,k) result (c)
      use interfaces, only: fdigamma
      double precision, intent(in) :: h, k
      double precision d1K, t
      if (h .eq. 0d0) then
        c = 0d0
      else
        d1K = logbesselk_dk(h,k)
        t = fdigamma(k) - log(.5d0*h)
        c = matern(h,k)*(d1K - t)
      end if
    end function matern_dk
  end function cor_dk_matern

  pure double precision function cor_dhdk_matern (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      c = matern_dhdk(h,kappa)
    end if
  contains
    elemental double precision function matern_dhdk (h,k) result (c)
      use interfaces, only: fbesselk, fgamma, fdigamma
      double precision, intent(in) :: h, k
      double precision d1K, t
      if (h .eq. 0d0) then
        c = 0d0
      else
        d1K = logbesselk_dk(h,k-1d0)
        t = fdigamma(k) - log(.5d0*h)
        c = -2d0/fgamma(k) * (.5d0*h)**k * fbesselk(h,k-1d0) * (d1K - t)
      end if
    end function matern_dhdk
  end function cor_dhdk_matern

  pure double precision function cor_hk_matern (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      c = matern_hk(h,kappa)
    end if
  contains
    elemental double precision function matern_hk (h,k) result (c)
      use interfaces, only: fdigamma, ftrigamma
      double precision, intent(in) :: h, k
      double precision d1K, d2K, t
      if (h .eq. 0d0) then
        c = 0d0
      else
        d1K = logbesselk_dk(h,k)
        d2K = logbesselk_hk(h,k) + d1K*d1K
        t = fdigamma(k) - log(.5d0*h)
        c = matern(h,k)*(d2K - (t+t)*d1K + t*t - ftrigamma(k))
      end if
    end function matern_hk
  end function cor_hk_matern


!!! Spherical
  pure double precision function cor_spher (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 1d0
    else if (h .ge. 1d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      c = 1d0 - 1.5d0*h + .5d0*h*h*h
    else
      c = bigneg
    end if
  end function cor_spher

  pure double precision function cor_dh_spher (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = -1.5d0
    else if (h .gt. 1d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      c = -1.5d0 + 1.5d0*h*h
    else
      c = bigneg
    end if
  end function cor_dh_spher

  pure double precision function cor_dk_spher (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .ge. 0d0) then
      c = 0d0
    else
      c = bigneg
    end if
  end function cor_dk_spher

  pure double precision function cor_hh_spher (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 1d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      c = 3d0*h
    else
      c = bigneg
    end if
  end function cor_hh_spher

  pure double precision function cor_dhdk_spher (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .ge. 0d0) then
      c = 0d0
    else
      c = bigneg
    end if
  end function cor_dhdk_spher

  pure double precision function cor_hk_spher (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .ge. 0d0) then
      c = 0d0
    else
      c = bigneg
    end if
  end function cor_hk_spher

!!! Power-exponential
  pure double precision function cor_powexp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 1d0
    else if (h .gt. 0d0) then
      if (kappa .eq. 1d0) then
        c = exp(-h)
      else if (kappa .eq. 2d0) then
        c = exp(-h*h)
      else if (kappa .gt. 0d0 .and. kappa .le. 2d0) then
        c = exp(-(h**kappa))
      else
        c = bigneg
      end if
    else
      c = bigneg
    end if
  end function cor_powexp

  pure double precision function cor_dh_powexp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      if (kappa .eq. 1d0) then
        c = -exp(-h)
      else if (kappa .eq. 2d0) then
        c = -2d0*h*exp(-h*h)
      else if (kappa .gt. 0d0 .and. kappa .le. 2d0) then
        c = -kappa*(h**(kappa-1d0))*exp(-(h**kappa))
      else
        c = bigneg
      end if
    else
      c = bigneg
    end if
  end function cor_dh_powexp

  pure double precision function cor_dk_powexp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    double precision hk
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      if (kappa .eq. 1d0) then
        c = -h*log(h)*exp(-h)
      else if (kappa .eq. 2d0) then
        c = -h*h*log(h)*exp(-h*h)
      else if (kappa .gt. 0d0 .and. kappa .le. 2d0) then
        hk = h**kappa
        c = -hk*log(h)*exp(-hk)
      else
        c = bigneg
      end if
    else
      c = bigneg
    end if
  end function cor_dk_powexp

  pure double precision function cor_dhdk_powexp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    double precision lh, hk1, hk0, hk1lh, h2k1lh
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      lh = log(h)
      if (kappa .eq. 1d0) then
        c = (h*lh - lh - 1d0)*exp(-h)
      else if (kappa .eq. 2d0) then
        hk0 = h*h
        hk1 = hk0*h
        c = (2d0*lh*hk1 - 2d0*lh*h - h)*exp(-hk0)
      else if (kappa .gt. 0d0 .and. kappa .le. 2d0) then
        hk1 = h**(kappa -1d0)
        hk0 = hk1*h
        hk1lh = kappa*lh*hk1
        h2k1lh = hk1lh*hk0
        c = (h2k1lh - hk1lh - hk1)*exp(-hk0)
      else
        c = bigneg
      end if
    else
      c = bigneg
    end if
  end function cor_dhdk_powexp

  pure double precision function cor_hh_powexp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    double precision hk2, hk1, hk0
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      if (kappa .eq. 1d0) then
        c = exp(-h)
      else if (kappa .eq. 2d0) then
        c = 2d0*(2d0*h*h-1d0)*exp(-h*h)
      else if (kappa .gt. 0d0 .and. kappa .le. 2d0) then
        hk2 = h**(kappa - 2d0)
        hk1 = hk2*h
        hk0 = hk1*h
        c = (kappa*kappa*hk1*hk1 - kappa*(kappa - 1d0)*hk2)*exp(-hk0)
      else
        c = bigneg
      end if
    else
      c = bigneg
    end if
  end function cor_hh_powexp

  pure double precision function cor_hk_powexp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    double precision hk, lh2
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      lh2 = log(h)
      lh2 = lh2*lh2
      if (kappa .eq. 1d0) then
        c = lh2*h*(h-1d0)*exp(-h)
      else if (kappa .eq. 2d0) then
        hk = h*h
        c = lh2*hk*(hk-1d0)*exp(-hk)
      else if (kappa .gt. 0d0 .and. kappa .le. 2d0) then
        hk = h**kappa
        c = lh2*hk*(hk-1d0)*exp(-hk)
      else
        c = bigneg
      end if
    else
      c = bigneg
    end if
  end function cor_hk_powexp

!!! Exponential
  pure double precision function cor_exp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 1d0
    else if (h .gt. 0d0) then
      c = exp(-h)
    else
      c = bigneg
    end if
  end function cor_exp

  pure double precision function cor_dh_exp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      c = -exp(-h)
    else
      c = bigneg
    end if
  end function cor_dh_exp

  pure double precision function cor_dk_exp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .ge. 0d0) then
      c = 0d0
    else
      c = bigneg
    end if
  end function cor_dk_exp

  pure double precision function cor_dhdk_exp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .ge. 0d0) then
      c = 0d0
    else
      c = bigneg
    end if
  end function cor_dhdk_exp

  pure double precision function cor_hh_exp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      c = exp(-h)
    else
      c = bigneg
    end if
  end function cor_hh_exp

  pure double precision function cor_hk_exp (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .ge. 0d0) then
      c = 0d0
    else
      c = bigneg
    end if
  end function cor_hk_exp

!!! Gaussian
  pure double precision function cor_gaussian (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 1d0
    else if (h .gt. 0d0) then
      c = exp(-h*h)
    else
      c = bigneg
    end if
  end function cor_gaussian

  pure double precision function cor_dh_gaussian (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      c = -2d0*h*exp(-h*h)
    else
      c = bigneg
    end if
  end function cor_dh_gaussian

  pure double precision function cor_dk_gaussian (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .ge. 0d0) then
      c = 0d0
    else
      c = bigneg
    end if
  end function cor_dk_gaussian

  pure double precision function cor_dhdk_gaussian (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .ge. 0d0) then
      c = 0d0
    else
      c = bigneg
    end if
  end function cor_dhdk_gaussian

  pure double precision function cor_hh_gaussian (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .eq. 0d0) then
      c = 0d0
    else if (h .gt. 0d0) then
      c = 2d0*(2d0*h*h-1d0)*exp(-h*h)
    else
      c = bigneg
    end if
  end function cor_hh_gaussian

  pure double precision function cor_hk_gaussian (h,kappa) result (c)
    implicit none
    double precision, intent(in) :: h, kappa
    if (h .ge. 0d0) then
      c = 0d0
    else
      c = bigneg
    end if
  end function cor_hk_gaussian
end module cor_fcns

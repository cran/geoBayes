!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module linkfcn
  implicit none
  double precision, parameter :: bigpos = huge(1d0), bigneg = -bigpos, &
     smallpos = epsilon(1d0), smallneg = -smallpos
  private bigpos, bigneg, smallpos, smallneg
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link function
  elemental function flink_bi (mu, d) result (z)
    use interfaces, only: quantt, quantnorm, quantlogis
    implicit none
    double precision, intent(in) :: mu, d
    double precision z
    if (d .gt. 0d0) then
      z = quantt(mu, d)
    else if (d .lt. 0d0) then
      z = quantlogis(mu)
    else
      z = quantnorm(mu)
    end if
  end function flink_bi

!! Inverse link function (log scale)
  elemental function invlink_bi (z, d) result (y)
    ! Binomial symmetric (robit) link fcn
    use interfaces, only: logprobt, logprobnorm, logproblogis
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    if (d .gt. 0d0) then
      y = logprobt(z, d)
    else if (d .lt. 0d0) then
      y = logproblogis(z)
    else
      y = logprobnorm(z)
    end if
  end function invlink_bi

!! 1st derivative w.r.t. z
  pure function invlink1_bi (z, d) result (y)
    use interfaces, only: logprobt, logprobnorm, logpdft, logpdfnorm
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision f, ff
    if (d .gt. 0d0) then
      ff = logprobt(z, d)
      f = logpdft(z, d)
      y = exp(f - ff)
    else if (d .lt. 0d0) then
      f = 0.6458
      y = 1d0/(f*(1d0 + exp(z/f)))
    else
      ff = logprobnorm(z)
      f = logpdfnorm(z)
      y = exp(f - ff)      
    end if
  end function invlink1_bi

!! 2nd derivative w.r.t. z
  pure function invlink2_bi (z, d) result (y)
    use interfaces, only: logprobt, logprobnorm, logpdft, logpdfnorm
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision f, ff
    if (d .gt. 0d0) then
      ff = logprobt(z, d)
      f = logpdft(z, d)
      y = f/ff
      y = -y*(z+y)
    else if (d .lt. 0d0) then
      f = 0.6458
      ff = f + f
      y = z/ff
      y = 1d0/(cosh(y)*ff)
      y = -y*y
    else
      ff = logprobnorm(z)
      f = logpdfnorm(z)
      y = f/ff
      y = -y*((d+1d0)*z/(d+z*z) + y)
    end if
  end function invlink2_bi

!! 3rd derivative w.r.t. z
  pure function invlink3_bi (z, d) result (y)
    use interfaces, only: logprobt, logprobnorm, logpdft, logpdfnorm
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision f, ff
    if (d .gt. 0d0) then
      f = invlink1_bi(z, d)
      ff = invlink2_bi(z, d)
      y = d + z*z
      y = (d+1d0)*(d + d - y)/(y*y)*f + (d+1d0)*z/y*ff + 2d0*f*ff
      y = -y
    else if (d .lt. 0d0) then
      f = 0.6458
      y = z/f
      ff = .5d0*y
      y = 1d0/(f*sinh(y))
      ff = sinh(ff)
      ff = ff*ff
      y = y*y*y*ff*ff
      y = y + y
    else
      f = invlink1_bi(z, d)
      ff = invlink2_bi(z, d)
      y = f + z*ff + 2d0*f*ff
      y = -y
    end if
  end function invlink3_bi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GEV binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link function
!! mu = 1 - exp{-max(0,1+d*z)**(1/d)} , if d /= 0
!!    = 1 - exp(-exp(z))              , if d == 0
  elemental function flink_ba (mu, d) result (z)
    ! mu is log(1-mu)
    implicit none
    double precision, intent(in) :: mu, d
    double precision z
      z = -mu
      if (d .eq. 0d0) then
        z = log(z)
      else if (d .eq. .5d0) then
        z = -1d0 + sqrt(z)
        z = z + z
      else if (d .eq. 1d0) then
        z = z - 1d0
      else if (d .eq. 2d0) then
        z = .5d0*(z*z - 1d0)
      else if (d .eq. -.5d0) then
        z = -1d0 + 1d0/sqrt(z)
        z = z + z
      else if (d .eq. -1d0) then
        z = 1d0 - 1d0/z
      else if (d .eq. -2d0) then
        z = .5d0 - 1d0/(z*z)
      else
        z = (-1d0 + z**d)/d
      end if
  end function flink_ba

!! Inverse link fcn. Returns log(1-mu)
  elemental function invlink_ba (z, d) result (y)
    ! Binomial asymmetric (GEV) link fcn
    ! mu = 1 - exp{-max(0,1+d*z)**(1/d)} , if d /= 0
    !    = 1 - exp(-exp(z))              , if d == 0
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    if (d .eq. 0d0) then
      y = -exp(z)
    else if (d .eq. 0.5d0) then
      y = 1d0 + 0.5d0*z
      if (y .gt. 0d0) then
        y = -y*y
      else
        y = smallneg
      end if
    else if (d .eq. -0.5d0) then
      y = 1d0 - 0.5d0*z
      if (y .gt. 0d0) then
        y = -1d0/(y*y)
      else
        y = bigneg
      end if
    else if (d .eq. 1d0) then
      y = 1d0 + z
      if (y .gt. 0d0) then
        y = -y
      else
        y = smallneg
      end if
    else if (d .eq. -1d0) then
      y = 1d0 - z
      if (y .gt. 0d0) then
        y = -1d0/y
      else
        y = bigneg
      end if
    else if (d .eq. 2d0) then
      y = 1d0 + z + z
      if (y .gt. 0d0) then
        y = -sqrt(y)
      else
        y = smallneg
      end if
    else if (d .eq. -2d0) then
      y = 1d0 - z - z
      if (y .gt. 0d0) then
        y = -1d0/sqrt(y)
      else
        y = bigneg
      end if
    else if (d .gt. 0d0) then 
      y = 1d0 + d*z
      if (y .gt. 0d0) then
        y = -y**(1d0/d)
      else
        y = smallneg
      end if
    else
      y = 1d0 + d*z
      if (y .gt. 0d0) then
        y = -y**(1d0/d)
      else
        y = bigneg
      end if
    end if
  end function invlink_ba

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GEVD binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link function
!! mu = 1 - exp{-max(0, 1+d*z)**(-1/d)} , if d /= 0
!!    = 1 - exp(-exp(-z))               , if d == 0
  elemental function flink_bd (mu, d) result (z)
    ! mu is log (1-mu)
    implicit none
    double precision, intent(in) :: mu, d
    double precision z
    z = -mu
    if (d .eq. 0d0) then
      z = -log(z)
    else ! d .ne. 0d0
      z = (-1d0 + z**(-d))/d
    end if
  end function flink_bd
  
!! Inverse link fcn. Returns log(1-mu)
  elemental function invlink_bd (z, d) result (y)
    ! Binomial asymmetric (GEV) decreasing link fcn
    ! p = 1 - exp{-max(0, 1+d*z)**(-1/d)}
    ! y = log(1-p) = -max(0, 1+d*z)**(-1/d)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    if (d .eq. 0d0) then
      y = -exp(-z)
    else
      y = 1d0 + d*z
      if (y .gt. 0d0) then
        y = -y**(-1d0/d)
      else if (d .gt. 0d0) then
        y = bigneg
      else
        y = smallneg
      end if
    end if
  end function invlink_bd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Poisson !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link fcn
  elemental function flink_po (mu, d) result (z)
    ! mu is log(mu)
    use interfaces, only: fexpm1
    implicit none
    double precision, intent(in) :: mu, d
    double precision z
    if (d .eq. 0d0) then
      z = mu
    else if (d .eq. 1d0) then
      z = fexpm1(mu)
    else if (d .eq. -1d0) then
      z = -fexpm1(-mu)
    else
      z = fexpm1(mu*d)/d
    end if
  end function flink_po
  
!! Inverse link fcn, returns log(mu)
  elemental function invlink_po (z,d) result (y)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    if (d .eq. 0d0) then
      y = z
    else 
      y = d*z
      if (y .gt. -1d0) then
        y = flog1p(y)/d
      else
        y = bigneg
      end if
    end if
  end function invlink_po

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gamma !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link fcn
  elemental function flink_gm (mu, d) result (z)
    implicit none
    double precision, intent(in) :: mu, d
    double precision z
    z = flink_po(mu, d)
  end function flink_gm
  
!! Inverse link fcn, returns log(mu)
  elemental function invlink_gm (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = invlink_po(z, d)
  end function invlink_gm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gaussian !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Linf fcn
  elemental function flink_ga (mu, d) result (z)
    implicit none
    double precision, intent(in) :: mu, d
    double precision z
    if (d .eq. 0d0) then
      if (mu > 0d0) then
        z = log(mu)
      else
        z = bigneg
      end if
    else if (d .gt. 0d0) then
      if (d .eq. 2d0) then
        z = .5d0*(mu*mu - 1d0)
      else if (d .eq. 1d0) then
        z = mu - 1d0
      else if (d .eq. .5d0) then
        z = sign(sqrt(abs(mu)), mu) - 1d0
        z = z + z
      else
        z = (sign(abs(mu)**d, mu) - 1d0)/d
      end if
    else ! d .lt. 0d0
      if (mu .gt. 0d0) then
        if (d .eq. -0.5d0) then
          z = -(1d0/sqrt(mu) - 1d0)
          z = z + z
        else if (d .eq. -1d0) then
          z = -(1d0/mu - 1d0)
        else if (d .eq. -2d0) then
          z = -.5d0*(1d0/(mu*mu) - 1d0)
        else
          z = (mu**d - 1d0)/d
        end if
      else ! mu .le. 0d0
        z = bigneg
      end if
    end if
  end function flink_ga
  
!! Inverse link fcn 
  elemental function invlink_ga (z,d) result(y)
    ! Inverse Box-Cox transformation. Using extended if d > 0.
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    if (d .eq. 0d0) then
      y = exp(z)
    else if (d .gt. 0d0) then
      y = d*z + 1d0
      if (d .eq. 2d0) then
        y = sign(sqrt(abs(y)),y)
      else if (d .eq. .5d0) then
        y = sign(y*y,y)
      else if (d .ne. 1d0) then
        y = sign(abs(y)**(1d0/d),y)
      end if
    else ! d .lt. 0d0
      y = d*z + 1d0
      if (y .gt. 0d0) then
        y = 1d0/y
        if (d .eq. -2d0) then
          y = sqrt(y)
        else if (d .eq. -.5d0) then
          y = y*y
        else if (d .ne. -1d0) then
          y = y**(1d0/(-d))
        end if
      else
        y = bigpos
      end if
    end if
  end function invlink_ga
  
  
!!!!!!!!!!!!!!!!!!! Wallace transformation for binomial !!!!!!!!!!!!!!!!!!!
!!! See function u4 in
!!!  Wallace, D. L. (1959). Bounds on normal approximations to Student's
!!!  and the chi-square distributions. The Annals of Mathematical
!!!  Statistics, 1121-1130. 
!! Link fcn
  elemental function flink_bw (mu, d) result(z)
    use interfaces, only: quantnorm, fexpm1
    implicit none
    double precision, intent(in) :: mu, d
    double precision :: z
    double precision t, e
    t = 8d0*d
    e = quantnorm(mu)*(t + 3d0)/(t + 1d0)
    t = sqrt(d*fexpm1(e*e/d))
    z = sign(t, e)
  end function flink_bw

!! Inverse link fcn (log scale)
  elemental function invlink_bw (z,d) result (y)
    use interfaces, only: logprobnorm, flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision x, t
    t = 8d0*d
    x = sqrt(d*flog1p(z*z/d))*(t + 1d0)/(t + 3d0)
    t = sign(x, z)
    y = logprobnorm(t)
  end function invlink_bw
end module linkfcn


! Local Variables:
! compile-command: "R CMD SHLIB"
! End:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module linkfcns
  implicit none
  double precision, parameter :: bigpos = huge(1d0), bigneg = -bigpos, &
     smallpos = epsilon(1d0), smallneg = -smallpos
  private bigpos, bigneg, smallpos, smallneg
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link function
  elemental function flink_bi (w, d) result (z)
    use interfaces, only: quantt, quantnorm, quantlogis
    implicit none
    double precision, intent(in) :: w, d
    double precision z
    if (d .gt. 0d0) then
      z = quantt(w, d)
    else if (d .lt. 0d0) then
      z = quantlogis(w)
    else
      z = quantnorm(w)
    end if
  end function flink_bi

!! Inverse link function (log scale)
  elemental function invlink_bi (z, d) result (w)
    ! Binomial symmetric (robit) link fcn
    use interfaces, only: logprobt, logprobnorm, logproblogis
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .gt. 0d0) then
      w = logprobt(z, d)
    else if (d .lt. 0d0) then
      w = logproblogis(z)
    else
      w = logprobnorm(z)
    end if
  end function invlink_bi

!! Derivative of inverse link function w.r.t. parameter
  function invlinkdnu_bi (z, d) result (y)
    ! Binomial symmetric (robit) link fcn
    use interfaces, only: logprobt
    use tcdfder
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    if (d .gt. 0d0) then
      y = tcdfdnu(z, d)
    else 
      y = 0d0
    end if
  end function invlinkdnu_bi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GEV binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link function
!! mu = 1 - exp{-max(0,1+d*z)**(1/d)} , if d /= 0
!!    = 1 - exp(-exp(z))              , if d == 0
  elemental function flink_ba (w, d) result (z)
    ! w is log(1-mu)
    implicit none
    double precision, intent(in) :: w, d
    double precision z
      z = -w
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
  elemental function invlink_ba (z, d) result (w)
    ! Binomial asymmetric (GEV) link fcn
    ! mu = 1 - exp{-max(0,1+d*z)**(1/d)} , if d /= 0
    !    = 1 - exp(-exp(z))              , if d == 0
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = -exp(z)
    else if (d .eq. 0.5d0) then
      w = 1d0 + 0.5d0*z
      if (w .gt. 0d0) then
        w = -w*w
      else
        w = smallneg
      end if
    else if (d .eq. -0.5d0) then
      w = 1d0 - 0.5d0*z
      if (w .gt. 0d0) then
        w = -1d0/(w*w)
      else
        w = bigneg
      end if
    else if (d .eq. 1d0) then
      w = 1d0 + z
      if (w .gt. 0d0) then
        w = -w
      else
        w = smallneg
      end if
    else if (d .eq. -1d0) then
      w = 1d0 - z
      if (w .gt. 0d0) then
        w = -1d0/w
      else
        w = bigneg
      end if
    else if (d .eq. 2d0) then
      w = 1d0 + z + z
      if (w .gt. 0d0) then
        w = -sqrt(w)
      else
        w = smallneg
      end if
    else if (d .eq. -2d0) then
      w = 1d0 - z - z
      if (w .gt. 0d0) then
        w = -1d0/sqrt(w)
      else
        w = bigneg
      end if
    else if (d .gt. 0d0) then 
      w = 1d0 + d*z
      if (w .gt. 0d0) then
        w = -w**(1d0/d)
      else
        w = smallneg
      end if
    else
      w = 1d0 + d*z
      if (w .gt. 0d0) then
        w = -w**(1d0/d)
      else
        w = bigneg
      end if
    end if
  end function invlink_ba

!! Inverse link fcn derivative w.r.t. nu
  elemental function invlinkdnu_ba (z, d) result (y)
    ! Binomial asymmetric (GEV) link fcn
    ! mu = 1 - exp{-max(0,1+d*z)**(1/d)} , if d /= 0
    !    = 1 - exp(-exp(z))              , if d == 0
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    if (d .eq. 0d0) then
      y = .5d0*exp(z)*z*z
    else if (d .eq. 0.5d0) then
      y = 2d0 + z
      if (y .gt. 0d0) then
        y = y*(y*log(.5d0*y) - z)
      else
        y = 0d0
      end if
    else if (d .eq. -0.5d0) then
      y = 2d0 - z
      if (y .gt. 0d0) then
        y = 16d0 * (y*log(.5d0*y) + z)/(y*y*y)
      else
        y = bigpos
      end if
    else if (d .eq. 1d0) then
      y = 1d0 + z
      if (y .gt. 0d0) then
        y = y*log(y) - z
      else
        y = -z
      end if
    else if (d .eq. -1d0) then
      y = 1d0 - z
      if (y .gt. 0d0) then
        y = (y*log(y) + z)/(y*y)
      else
        y = bigpos
      end if
    else if (d .eq. 2d0) then
      y = 1d0 + z + z
      if (y .gt. 0d0) then
        y = (y*log(y) - z - z)/(4d0*sqrt(y))
      else
        y = bigpos
      end if
    else if (d .eq. -2d0) then
      y = 1d0 - z - z
      if (y .gt. 0d0) then
        y = (y*log(y) + z + z)/(4d0*y*sqrt(y))
      else
        y = bigpos
      end if
    else
      y = 1d0 + d*z
      if (y .gt. 0d0) then
        y = (y*log(y) - d*z)/(d*d*y**(1d0-1d0/d))
      else
        y = bigpos
      end if
    end if
  end function invlinkdnu_ba

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GEVD binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link function
!! mu = 1 - exp{-max(0, 1+d*z)**(-1/d)} , if d /= 0
!!    = 1 - exp(-exp(-z))               , if d == 0
  elemental function flink_bd (w, d) result (z)
    ! w is log (1-mu)
    implicit none
    double precision, intent(in) :: w, d
    double precision z
    z = -w
    if (d .eq. 0d0) then
      z = -log(z)
    else ! d .ne. 0d0
      z = (-1d0 + z**(-d))/d
    end if
  end function flink_bd
  
!! Inverse link fcn. Returns log(1-mu)
  elemental function invlink_bd (z, d) result (w)
    ! Binomial asymmetric (GEV) decreasing link fcn
    ! p = 1 - exp{-max(0, 1+d*z)**(-1/d)}
    ! w = log(1-p) = -max(0, 1+d*z)**(-1/d)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = -exp(-z)
    else
      w = 1d0 + d*z
      if (w .gt. 0d0) then
        w = -w**(-1d0/d)
      else if (d .gt. 0d0) then
        w = bigneg
      else
        w = smallneg
      end if
    end if
  end function invlink_bd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Box-Cox !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link fcn
  elemental function flink_boxcox (w, d) result (z)
    ! z = (mu^d - 1)/d
    ! w is log(mu)
    use interfaces, only: fexpm1
    implicit none
    double precision, intent(in) :: w, d
    double precision z
    if (d .eq. 0d0) then
      z = w
    else if (d .eq. 1d0) then
      z = fexpm1(w)
    else if (d .eq. -1d0) then
      z = -fexpm1(-w)
    else
      z = fexpm1(w*d)/d
    end if
  end function flink_boxcox
  
!! Inverse link fcn, returns log(mu)
  elemental function invlink_boxcox (z,d) result (w)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = z
    else 
      w = d*z
      if (w .gt. -1d0) then
        w = flog1p(w)/d
      else
        w = bigneg
      end if
    end if
  end function invlink_boxcox

!!   elemental function invlinkdnu_boxcox (z,d) result (y)
!!     implicit none
!!     double precision, intent(in) :: z, d
!!     double precision y
!!     if (d .eq. 0d0) then
!!       y = -.5d0*z*z
!!     else 
!!       y = 1d0 + d*z
!!       if (y .gt. 0d0) then
!!         y = z/(d*y) - log(y)/(d*d)
!!       else
!!         y = bigneg
!!       end if
!!     end if
!!   end function invlinkdnu_boxcox

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Poisson !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link fcn
  elemental function flink_po (w, d) result (z)
    ! w is log(mu)
    ! z = sign(w)*expm1(d*|w|)/d for d > 0
    use interfaces, only: fexpm1
    implicit none
    double precision, intent(in) :: w, d
    double precision z
    if (d .eq. 0d0) then
      z = w
    else if (d .eq. 1d0) then
      z = fexpm1(abs(w))
      z = sign(z,w)
    else if (d .gt. 0d0) then
      z = fexpm1(d*abs(w))/d
      z = sign(z,w)
    else ! Use regular Box-Cox
      z = fexpm1(w*d)/d
    end if
  end function flink_po

  elemental function invlink_po (z,d) result (w)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = z
    else if (d .eq. 1d0) then
      w = flog1p(abs(z))
      w = sign(w,z)
    else if (d .gt. 0d0) then
      w = flog1p(d*abs(z))/d
      w = sign(w,z)
    else ! Use regular Box-Cox
      w = d*z
      if (w .gt. -1d0) then
        w = flog1p(w)/d
      else
        w = bigneg
      end if
    end if
  end function invlink_po

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gamma !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link fcn
  elemental function flink_gm (w, d) result (z)
    ! z = (mu^d - 1)/d
    ! w is log(mu)
    implicit none
    double precision, intent(in) :: w, d
    double precision z
    z = flink_po(w, d)
  end function flink_gm
  
!! Inverse link fcn, returns log(mu)
  elemental function invlink_gm (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = invlink_po(z, d)
  end function invlink_gm

!!   elemental function invlinkdnu_gm (z,d) result (y)
!!     implicit none
!!     double precision, intent(in) :: z, d
!!     double precision y
!!     y = invlinkdnu_po(z, d)
!!   end function invlinkdnu_gm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gaussian !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Linf fcn
  elemental function flink_ga (w, d) result (z)
    implicit none
    double precision, intent(in) :: w, d
    double precision z
    if (d .eq. 0d0) then
      if (w > 0d0) then
        z = log(w)
      else
        z = bigneg
      end if
    else if (d .gt. 0d0) then
      if (d .eq. 2d0) then
        z = .5d0*(w*w - 1d0)
      else if (d .eq. 1d0) then
        z = w - 1d0
      else if (d .eq. .5d0) then
        z = sign(sqrt(abs(w)), w) - 1d0
        z = z + z
      else
        z = (sign(abs(w)**d, w) - 1d0)/d
      end if
    else ! d .lt. 0d0
      if (w .gt. 0d0) then
        if (d .eq. -0.5d0) then
          z = -(1d0/sqrt(w) - 1d0)
          z = z + z
        else if (d .eq. -1d0) then
          z = -(1d0/w - 1d0)
        else if (d .eq. -2d0) then
          z = -.5d0*(1d0/(w*w) - 1d0)
        else
          z = (w**d - 1d0)/d
        end if
      else ! w .le. 0d0
        z = bigneg
      end if
    end if
  end function flink_ga
  
!! Inverse link fcn 
  elemental function invlink_ga (z,d) result(w)
    ! Inverse Box-Cox transformation. Using extended if d > 0.
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = exp(z)
    else if (d .gt. 0d0) then
      w = d*z + 1d0
      if (d .eq. 2d0) then
        w = sign(sqrt(abs(w)),w)
      else if (d .eq. .5d0) then
        w = sign(w*w,w)
      else if (d .ne. 1d0) then
        w = sign(abs(w)**(1d0/d),w)
      end if
    else ! d .lt. 0d0
      w = d*z + 1d0
      if (w .gt. 0d0) then
        w = 1d0/w
        if (d .eq. -2d0) then
          w = sqrt(w)
        else if (d .eq. -.5d0) then
          w = w*w
        else if (d .ne. -1d0) then
          w = w**(1d0/(-d))
        end if
      else
        w = bigpos
      end if
    end if
  end function invlink_ga
  
  
!!!!!!!!!!!!!!!!!!! Wallace transformation for binomial !!!!!!!!!!!!!!!!!!!
!!! See function u4 in
!!!  Wallace, D. L. (1959). Bounds on normal approximations to Student's
!!!  and the chi-square distributions. The Annals of Mathematical
!!!  Statistics, 1121-1130.
!! Link fcn
!! mu = PHI[sign(z) * c(nu) * sqrt(nu*log(1+z*z/nu))]
!! where c(nu) = (8*nu+1)/(8*nu+3)
!! w = log(mu)
  elemental function flink_bw (w, d) result(z)
    use interfaces, only: quantnorm, fexpm1
    implicit none
    double precision, intent(in) :: w, d
    double precision :: z
    double precision t, e
    t = 8d0*d
    e = quantnorm(w)*(t + 3d0)/(t + 1d0)
    t = sqrt(d*fexpm1(e*e/d))
    z = sign(t, e)
  end function flink_bw

!! Inverse link fcn (log scale)
  elemental function invlink_bw (z,d) result (w)
    use interfaces, only: logprobnorm, flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision x, t
    t = 8d0*d
    x = sqrt(d*flog1p(z*z/d))*(t + 1d0)/(t + 3d0)
    t = sign(x, z)
    w = logprobnorm(t)
  end function invlink_bw
end module linkfcns


! Local Variables:
! compile-command: "R CMD SHLIB"
! End:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Commentary:
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module modelfcns_link
  implicit none
  double precision, parameter :: bigpos = huge(1d0), bigneg = -bigpos, &
     smallpos = epsilon(1d0), smallneg = -smallpos
  private bigpos, bigneg, smallpos, smallneg
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! binomial robit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link function
  pure function flink_robit (w,d) result (z)
    use interfaces, only: quantt
    implicit none
    double precision, intent(in) :: w, d
    double precision z
    z = quantt(w, d)
  end function flink_robit

!! Inverse link function (log scale)
  pure function invlink_robit (z,d) result (w)
    ! Binomial symmetric (robit) link fcn
    use interfaces, only: logprobt
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = logprobt(z, d)
  end function invlink_robit

  pure function invlinkdz_robit (z,d) result (y)
    use interfaces, only: logprobt, logpdft
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision f, ff
    ff = logprobt(z, d)
    f = logpdft(z, d)
    y = exp(f - ff)
  end function invlinkdz_robit

  pure function invlinkhz_robit (z,d) result (y)
    ! If g_d(z) is the pdf of the t_d distribution, then we need
    ! d2w = d(log g_d(z)) * d1w - (d1w)^2
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision dwdz, dlogtdz
    dwdz = invlinkdz_robit(z,d)
    dlogtdz = -(d+1d0)*z/(d+z*z)
    y = dlogtdz*dwdz - dwdz*dwdz
  end function invlinkhz_robit

  pure function invlink3z_robit (z,d) result (y)
    ! If g_d(z) is the pdf of the t_d distribution, then we need
    ! d3w = d2(log g_d(z)) * d1w - (d1w)*(d2w) + (d2w)^2/(d1w)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision dw, hw, hlogt, zsq, dpzsq
    dw = invlinkdz_robit(z,d)
    hw = invlinkhz_robit(z,d)
    zsq = z*z
    dpzsq = d + zsq
    hlogt = -(d+1d0)/dpzsq*(1d0 - 2*zsq/dpzsq)
    y = hlogt*dw - dw*hw + hw*hw/dw
  end function invlink3z_robit

  pure function invlinkdn_robit (z,d) result (y)
    ! First derivative w.r.t. d
    double precision, intent(in) :: z, d
    double precision y
    double precision, parameter :: eps = sqrt(epsilon(1d0))
    y = 8d0*invlink_robit(z,d+eps) - 8d0*invlink_robit(z,d-eps) &
       - invlink_robit(z,d+eps*2d0) + invlink_robit(z,d-eps*2d0)
    y = y/(12d0*eps)
  end function invlinkdn_robit

  pure function invlinkhn_robit (z,d) result (y)
    ! First derivative w.r.t. d
    double precision, intent(in) :: z, d
    double precision y
    double precision, parameter :: eps = sqrt(sqrt(epsilon(1d0)))
    y = 16d0*invlink_robit(z,d+eps) + 16d0*invlink_robit(z,d-eps) &
       - invlink_robit(z,d+eps*2d0) - invlink_robit(z,d-eps*2d0) &
       - 30d0*invlink_robit(z,d)
    y = y/(12d0*eps*eps)
  end function invlinkhn_robit

  pure function invlinkdzdn_robit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision dz, dlogF, dlogt
    dz = invlinkdz_robit(z,d)
    dlogF = invlinkdn_robit(z,d)
    dlogt = dnu_logpdft(z,d)
    y = (dlogt - dlogF)*dz
  contains
    pure function dnu_logpdft (z,d)
      use interfaces, only: fdigamma, flog1p
      implicit none
      double precision, intent(in) :: z, d
      double precision dnu_logpdft
      double precision zzd
      zzd = z*z/d
      dnu_logpdft = -.5d0*flog1p(zzd) + .5d0*(d+1)/d*zzd/(zzd+1) &
         - .5d0/d -.5d0*fdigamma(.5d0*d) + .5d0*fdigamma(.5d0*(d+1))
    end function dnu_logpdft
  end function invlinkdzdn_robit

  pure function invlinkhzdn_robit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision dz, hz, dzdn, ddlogt
    dz = invlinkdz_robit(z,d)
    hz = invlinkhz_robit(z,d)
    dzdn = invlinkdzdn_robit(z,d)
    ddlogt = dzdn_logpdft(z,d)
    y = ddlogt*dz + hz*dzdn/dz - dzdn*dz
  contains
    pure function dzdn_logpdft (z,d)
      implicit none
      double precision, intent(in) :: z, d
      double precision dzdn_logpdft
      double precision zzd
      zzd = z*z+d
      !!dzdn_logpdft = (d+1)*z/(d*zzd) - z/zzd - (d+1)*z*z*z/(d*zzd*zzd)
      dzdn_logpdft = -z/zzd*(1d0 - (d+1)/d*(1d0 - z*z/zzd))
    end function dzdn_logpdft
  end function invlinkhzdn_robit

  pure function invlinkdzhn_robit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision dz, hlogt, hlogF, dlogg
    dz = invlinkdz_robit(z,d)
    hlogt = hnu_logpdft(z,d)
    hlogF = invlinkhn_robit(z,d)
    dlogg = invlinkdzdn_robit(z,d)
    dlogg = dlogg/dz
    y = dz*(dlogg*dlogg + hlogt - hlogF)
  contains
    pure function hnu_logpdft (z,d)
      use interfaces, only: ftrigamma, flog1p
      implicit none
      double precision, intent(in) :: z, d
      double precision hnu_logpdft
      double precision zzd
      zzd = z*z/d
      zzd = zzd/(1d0+zzd)
      hnu_logpdft = zzd/d - zzd*(d+1)/(d*d)*(1d0 - .5d0*zzd) + .5/(d*d) &
         - .25d0*ftrigamma(.5d0*d) + .25d0*ftrigamma(.5d0*(d+1))
    end function hnu_logpdft
  end function invlinkdzhn_robit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! binomial logit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link function
  pure function flink_logit (w,d) result (z)
    use interfaces, only: quantlogis
    implicit none
    double precision, intent(in) :: w, d
    double precision z
    z = quantlogis(w)
  end function flink_logit

!! Inverse link function (log scale)
  pure function invlink_logit (z,d) result (w)
    ! Binomial symmetric (robit) link fcn
    use interfaces, only: logproblogis
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = logproblogis(z)
  end function invlink_logit

  pure function invlinkdz_logit (z,d) result (y)
    use interfaces, only: logproblogis
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    ! y = 1d0/(1d0 + exp(z))
    y = logproblogis(z)
    y = exp(y)
    y = 1d0 - y
  end function invlinkdz_logit

  pure function invlinkhz_logit (z,d) result (y)
    use interfaces, only: logproblogis
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = logproblogis(z)
    y = exp(y)
    y = y*(y-1d0)
  end function invlinkhz_logit

  pure function invlink3z_logit (z,d) result (y)
    use interfaces, only: logproblogis
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = logproblogis(z)
    y = exp(y)
    y = y*(1d0-y)*(y+y-1d0)
  end function invlink3z_logit

  pure function invlinkdn_logit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = 0d0
  end function invlinkdn_logit

  pure function invlinkhn_logit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = 0d0
  end function invlinkhn_logit

  pure function invlinkdzdn_logit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = 0d0
  end function invlinkdzdn_logit

  pure function invlinkhzdn_logit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = 0d0
  end function invlinkhzdn_logit

  pure function invlinkdzhn_logit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = 0d0
  end function invlinkdzhn_logit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! binomial probit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link function
  pure function flink_probit (w,d) result (z)
    use interfaces, only: quantnorm
    implicit none
    double precision, intent(in) :: w, d
    double precision z
    z = quantnorm(w)
  end function flink_probit

!! Inverse link function (log scale)
  pure function invlink_probit (z,d) result (w)
    ! Binomial symmetric (robit) link fcn
    use interfaces, only: logprobnorm
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = logprobnorm(z)
  end function invlink_probit

  pure function invlinkdz_probit (z,d) result (y)
    use interfaces, only: logprobnorm, logpdfnorm
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision f, ff
    ff = logprobnorm(z)
    f = logpdfnorm(z)
    y = exp(f - ff)
  end function invlinkdz_probit

  pure function invlinkhz_probit (z,d) result (y)
    ! If g(z) is the pdf of the N(0,1) distribution, then we need
    ! d2w = d(log g(z)) * d1w - (d1w)^2
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision dwdz
    dwdz = invlinkdz_probit(z,d)
    y = -z*dwdz - dwdz*dwdz
  end function invlinkhz_probit

  pure function invlink3z_probit (z,d) result (y)
    ! If g(z) is the pdf of the N(0,1) distribution, then we need
    ! d3w = h(log g(z)) * dw - dw * hw + hw*hw/dw
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision dw, hw
    dw = invlinkdz_probit(z,d)
    hw = invlinkhz_probit(z,d)
    y = -dw - dw*hw + hw*hw/dw
  end function invlink3z_probit

  pure function invlinkdn_probit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = 0d0
  end function invlinkdn_probit

  pure function invlinkhn_probit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = 0d0
  end function invlinkhn_probit

  pure function invlinkdzdn_probit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = 0d0
  end function invlinkdzdn_probit

  pure function invlinkhzdn_probit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = 0d0
  end function invlinkhzdn_probit

  pure function invlinkdzhn_probit (z,d) result (y)
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = 0d0
  end function invlinkdzhn_probit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Box-Cox !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link fcn
  pure function flink_boxcox (w,d) result (z)
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
  pure function invlink_boxcox (z,d) result (w)
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
      else if (d .gt. 0d0) then
        w = bigneg
      else
        w = bigpos
      end if
    end if
  end function invlink_boxcox

  pure function invlinkdz_boxcox (z,d) result (w)
    ! w = log(mu)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = 1d0
    else
      w = d*z
      if (w .gt. -1d0) then
        w = 1d0/(1d0+w)
      else
        w = 0d0
      end if
    end if
  end function invlinkdz_boxcox

  pure function invlinkhz_boxcox (z,d) result (w)
    ! w = log(mu)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = 0d0
    else
      w = d*z
      if (w .gt. -1d0) then
        w = 1d0+w
        w = -d/(w*w)
      else
        w = 0d0
      end if
    end if
  end function invlinkhz_boxcox

  pure function invlink3z_boxcox (z,d) result (w)
    ! w = log(mu)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = 0d0
    else
      w = d*z
      if (w .gt. -1d0) then
        w = 1d0+w
        w = 2*d*d/(w*w*w)
      else
        w = 0d0
      end if
    end if
  end function invlink3z_boxcox

  pure function invlinkdn_boxcox (z,d) result (w)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = -.5d0*z*z
    else
      w = d*z
      if (w .gt. -1d0) then
        w = (w/(1d0+w) - flog1p(w))/(d*d)
      else
        w = 0d0
      end if
    end if
  end function invlinkdn_boxcox

  pure function invlinkhn_boxcox (z,d) result (w)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision zdp1
    if (d .eq. 0d0) then
      w = 2d0*z*z*z/3d0
    else
      w = d*z
      if (w .gt. -1d0) then
        zdp1 = w + 1d0
        w = -w*(3d0*zdp1 - 1d0) + 2d0*zdp1*zdp1*flog1p(w)
        w = w/(d*d*d*zdp1*zdp1)
      else
        w = 0d0
      end if
    end if
  end function invlinkhn_boxcox

  pure function invlinkdzdn_boxcox (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = -z
    else
      w = d*z
      if (w .gt. -1d0) then
        w = 1d0 + w
        w = -z/(w*w)
      else
        w = 0d0
      end if
    end if
  end function invlinkdzdn_boxcox

  pure function invlinkhzdn_boxcox (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w, dz
    if (d .eq. 0d0) then
      w = -1
    else
      dz = d*z
      if (dz .gt. -1d0) then
        w = 1d0 + dz
        w = (dz-1)/(w*w*w)
      else
        w = 0d0
      end if
    end if
  end function invlinkhzdn_boxcox

  pure function invlinkdzhn_boxcox (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = 2*z*z
    else
      w = d*z
      if (w .gt. -1d0) then
        w = 1d0 + w
        w = 2*z*z/(w*w*w)
      else
        w = 0d0
      end if
    end if
  end function invlinkdzhn_boxcox


!!!!!!!!!!!!!!!!!!!!!!!!!!!! Modified Box-Cox !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Link fcn
  pure function flink_modbc (w, d) result (z)
    ! w is log(mu)
    ! z = sign(w)*expm1(d*|w|)/d for d > 0
    use interfaces, only: fexpm1
    implicit none
    double precision, intent(in) :: w, d
    double precision z
    if (d .eq. 0d0) then
      z = w
    else
      z = fexpm1(abs(d*w))/abs(d)
      z = sign(z,w)
!     else if (d .eq. 1d0) then
!       z = fexpm1(abs(w))
!       z = sign(z,w)
!     else if (d .gt. 0d0) then
!       z = fexpm1(d*abs(w))/d
!       z = sign(z,w)
!     else ! Use regular Box-Cox
!       z = fexpm1(w*d)/d
    end if
  end function flink_modbc

  pure function invlink_modbc (z,d) result (w)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = z
    else
      w = flog1p(abs(d*z))/abs(d)
      w = sign(w,z)
!     else if (d .eq. 1d0) then
!       w = flog1p(abs(z))
!       w = sign(w,z)
!     else if (d .gt. 0d0) then
!       w = flog1p(d*abs(z))/d
!       w = sign(w,z)
!     else ! Use regular Box-Cox
!       w = d*z
!       if (w .gt. -1d0) then
!         w = flog1p(w)/d
!       else
!         w = bigneg
!       end if
    end if
  end function invlink_modbc

  pure function invlinkdz_modbc (z,d) result (w)
    ! w = log(mu)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = 1d0
    else
      w = 1d0/(1d0 + abs(d*z))
!     else if (d .gt. 0d0) then
!       w = 1d0 + d*abs(z)
!       w = 1d0/w
!     else ! Use regular Box-Cox
!       w = d*z
!       if (w .gt. -1d0) then
!         w = 1d0/(1d0+w)
!       else
!         w = 0d0
!       end if
    end if
  end function invlinkdz_modbc

  pure function invlinkhz_modbc (z,d) result (w)
    ! w = log(mu)
    implicit none
    double precision, intent(in) :: z, d
    double precision w, y
    if (d .eq. 0d0) then
      w = 0d0
    else
      w = d*z
      y = 1d0 + abs(w)
      if (w .lt. 0d0) then
        w = d/(y*y)
      else
        w = -d/(y*y)
      end if
    end if
  end function invlinkhz_modbc

  pure function invlink3z_modbc (z,d) result (w)
    ! w = log(mu)
    implicit none
    double precision, intent(in) :: z, d
    double precision w, y
    if (d .eq. 0d0) then
      w = 0d0
    else
      w = d*z
      y = 1d0 + abs(w)
      w = 2*d*d/(y*y*y)
    end if
  end function invlink3z_modbc

  pure function invlinkdn_modbc (z,d) result (w)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      ! The derivative at 0 is undefined
      w = 0d0
      ! w = -.5d0*z*z
    else
      w = abs(z*d)
      w = (w/(1d0+w) - flog1p(w))/(d*d)
      if (d .lt. 0d0) w = -w
    end if
    if (z .lt. 0d0) w = -w
  end function invlinkdn_modbc

  pure function invlinkhn_modbc (z,d) result (w)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision zdp1
    if (d .eq. 0d0) then
      w = 2d0*z*z*z/3d0
    else
      w = abs(d*z)
      zdp1 = w + 1d0
      w = -w*(3d0*zdp1 - 1d0) + 2d0*zdp1*zdp1*flog1p(w)
      w = w/(d*d*d*zdp1*zdp1)
      if (d .lt. 0d0) w = -w
      if (z .lt. 0d0) w = -w
    end if
  end function invlinkhn_modbc

  pure function invlinkdzdn_modbc (z,d) result (w)
    ! w = log(mu)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      ! The derivative at 0 is undefined
      w = 0d0
      !w = -abs(z)
    else
      w = 1d0 + abs(d*z)
      w = -abs(z)/(w*w)
      if (d .lt. 0d0) w = -w
    end if
  end function invlinkdzdn_modbc

  pure function invlinkhzdn_modbc (z,d) result (w)
    ! w = log(mu)
    implicit none
    double precision, intent(in) :: z, d
    double precision w, adz
    if (d .eq. 0d0) then
      ! The derivative at 0 is undefined
      w = 0d0
      !w = -1
    else
      adz = abs(d*z)
      w = 1d0 + adz
      w = (adz - 1)/(w*w*w)
      if (d .lt. 0d0) w = -w
    end if
    if (z .lt. 0d0) w = -w
  end function invlinkhzdn_modbc

  pure function invlinkdzhn_modbc (z,d) result (w)
    ! w = log(mu)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = 2*z*z
    else
      w = 1d0 + abs(d*z)
      w = 2*z*z/(w*w*w)
    end if
  end function invlinkdzhn_modbc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gaussian !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Linf fcn
  pure function flink_ga (w, d) result (z)
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
  pure function invlink_ga (z,d) result(w)
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

  pure function invlinkdz_ga (z,d) result(w)
    ! Inverse Box-Cox transformation. Using extended if d > 0.
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = exp(z)
    else if (d .eq. 1d0) then
      w = 1d0
    else if (d .gt. 0d0) then
      w = abs(d*z + 1d0)
      if (d .eq. 2d0) then
        w = 1d0/sqrt(w)
      else if (d .ne. .5d0) then
        w = w**(1d0/d - 1d0)
      end if
    else ! d .lt. 0d0
      w = d*z + 1d0
      if (w .gt. 0d0) then
        if (d .eq. -1d0) then
          w = 1d0/(w*w)
        else if (d .eq. -2d0) then
          w = 1d0/(w*sqrt(w))
        else if (d .eq. -.5d0) then
          w = 1d0/(w*w*w)
        else
          w = w**(1d0/d - 1d0)
        end if
      else
        w = 0d0
      end if
    end if
  end function invlinkdz_ga

  pure function invlinkhz_ga (z,d) result(w)
    ! Inverse Box-Cox transformation. Using extended if d > 0.
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = exp(z)
    else if (d .eq. 1d0) then
      w = 0d0
    else if (d .eq. .5d0) then
      w = .5d0
    else if (d .gt. 0d0) then
      w = d*z + 1d0
      if (w .ge. 0d0) then
        w = (1-d) * w**(1d0/d - 2d0)
      else
        w = (d-1) * (-w)**(1d0/d - 2d0)
      end if
    else ! d .lt. 0d0
      w = d*z + 1d0
      if (w .gt. 0d0) then
        w = (1-d) * w**(1d0/d - 2d0)
      else
        w = 0d0
      end if
    end if
  end function invlinkhz_ga

  pure function invlink3z_ga (z,d) result(w)
    ! Inverse Box-Cox transformation. Using extended if d > 0.
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = exp(z)
    else if (d .eq. 1d0) then
      w = 0d0
    else if (d .eq. .5d0) then
      w = 0d0
    else if (d .gt. 0d0) then
      w = d*z + 1d0
      w = (1-d)*(1-2*d) * abs(w)**(1d0/d - 3d0)
    else ! d .lt. 0d0
      w = d*z + 1d0
      if (w .gt. 0d0) then
        w = (1-d)*(1-2*d) * (w)**(1d0/d - 3d0)
      else
        w = 0d0
      end if
    end if
  end function invlink3z_ga

  pure function invlinkdn_ga (z,d) result(w)
    ! Inverse Box-Cox transformation. Using extended if d > 0.
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision zdp1
    if (d .eq. 0d0) then
      w = -.5d0*z*z*exp(z)
    else if (d .gt. 0d0) then
      w = d*z
      zdp1 = w + 1d0
      w = abs(zdp1)**(1d0/d-1d0)*(w - zdp1*log(abs(zdp1)))/(d*d)
    else ! d .lt. 0d0
      w = d*z
      zdp1 = w + 1d0
      if (zdp1 .gt. 0d0) then
        w = zdp1**(1d0/d-1d0)*(w - zdp1*log(zdp1))/(d*d)
      else
        w = 0d0
      end if
    end if
  end function invlinkdn_ga

  pure function invlinkhn_ga (z,d) result(w)
    ! Inverse Box-Cox transformation. Using extended if d > 0.
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision zdp1, zoz, t
    if (d .eq. 0d0) then
      w = z*z*z*(8d0 + 3d0*z)*exp(z)/12d0
    else if (d .gt. 0d0) then
      w = d*z
      zdp1 = w + 1d0
      zoz = w/zdp1
      t = zoz - log(abs(zdp1))
      w = abs(zdp1)**(1d0/d)*((-zoz*zoz - t - t) + t*t/d)/(d*d*d)
      if (zdp1 .lt. 0d0) w = -w
    else ! d .lt. 0d0
      w = d*z
      zdp1 = w + 1d0
      if (zdp1 .gt. 0d0) then
      zdp1 = w + 1d0
      zoz = w/zdp1
      t = zoz - log(zdp1)
      w = zdp1**(1d0/d)*((-zoz*zoz - t - t) + t*t/d)/(d*d*d)
      else
        w = 0d0
      end if
    end if
  end function invlinkhn_ga

  pure function invlinkdzdn_ga (z,d) result(w)
    ! Inverse Box-Cox transformation. Using extended if d > 0.
    implicit none
    double precision, intent(in) :: z, d
    double precision w, omod
    if (d .eq. 0d0) then
      w = -.5d0*exp(z)*z*(z+2)
    else if (d .eq. 1d0) then
      w = -log(abs(z+1d0))
    else if (d .gt. 0d0) then
      w = d*z + 1d0
      omod = 1d0 - 1d0/d
      if (w .lt. 0d0) then
        w = -w
        w = omod*z*(w**(-1d0-omod)) - (w**(-omod))*log(w)/(d*d)
      else
        w = -omod*z*(w**(-1d0-omod)) - (w**(-omod))*log(w)/(d*d)
      end if
    else ! d .lt. 0d0
      w = d*z + 1d0
      omod = 1d0 - 1d0/d
      if (w .gt. 0d0) then
        w = -omod*z*(w**(-1d0-omod)) - (w**(-omod))*log(w)/(d*d)
      else
        w = 0d0
      end if
    end if
  end function invlinkdzdn_ga

  pure function invlinkhzdn_ga (z,d) result(w)
    ! Inverse Box-Cox transformation. Using extended if d > 0.
    implicit none
    double precision, intent(in) :: z, d
    double precision w, omod
    if (d .eq. 0d0) then
      w = -.5d0*exp(z)*(2+4*z+z*z)
    else if (d .eq. 1d0) then
      w = -1d0/(1d0+z)
    else if (d .gt. 0d0) then
      w = d*z + 1d0
      omod = 1d0 - 1d0/d
      if (w .lt. 0d0) then
        w = -w
        w = (w**(-1d0-omod))*(1d0 + omod*d*((1+omod)*z/w - log(w)/(d*d)))
      else
        w = -(w**(-1d0-omod))*(1d0 - omod*d*((1+omod)*z/w + log(w)/(d*d)))
      end if
    else ! d .lt. 0d0
      w = d*z + 1d0
      omod = 1d0 - 1d0/d
      if (w .gt. 0d0) then
        w = -(w**(-1d0-omod))*(1d0 - omod*d*((1+omod)*z/w + log(w)/(d*d)))
      else
        w = 0d0
      end if
    end if
  end function invlinkhzdn_ga

  pure function invlinkdzhn_ga (z,d) result(w)
    ! Inverse Box-Cox transformation. Using extended if d > 0.
    implicit none
    double precision, intent(in) :: z, d
    double precision w, omod, zo1pzd
    if (d .eq. 0d0) then
      w = exp(z)*z*z*(24+20*z + 3*z*z)/12d0
    else if (d .eq. 1d0) then
      w = log(abs(z+1d0))
      w = w*w + 2*w - 2*z/(1+z)
    else if (d .gt. 0d0) then
      w = abs(d*z + 1d0)
      omod = 1d0 - 1d0/d
      zo1pzd = z/w
      w = w**(-omod)*(omod*zo1pzd*zo1pzd - 2*zo1pzd/(d*d) &
         + 2*log(w)/(d*d*d) + (omod*zo1pzd + log(w)/(d*d))**2)
    else ! d .lt. 0d0
      w = d*z + 1d0
      omod = 1d0 - 1d0/d
      zo1pzd = z/w
      if (w .gt. 0d0) then
      w = w**(-omod)*(omod*zo1pzd*zo1pzd - 2*zo1pzd/(d*d) &
         + 2*log(w)/(d*d*d) + (omod*zo1pzd + log(w)/(d*d))**2)
      else
        w = 0d0
      end if
    end if
  end function invlinkdzhn_ga


!!!!!!!!!!!!!!!!!!! Wallace transformation for binomial !!!!!!!!!!!!!!!!!!!
!!! See function u4 in
!!!  Wallace, D. L. (1959). Bounds on normal approximations to Student's
!!!  and the chi-square distributions. The Annals of Mathematical
!!!  Statistics, 1121-1130.
!! Link fcn
!! mu = PHI[sign(z) * c(nu) * sqrt(nu*log(1+z*z/nu))]
!! where c(nu) = (8*nu+1)/(8*nu+3)
!! w = log(mu)
  pure function flink_wallace (w, d) result(z)
    use interfaces, only: quantnorm, fexpm1
    implicit none
    double precision, intent(in) :: w, d
    double precision :: z
    double precision t, e
    t = 8d0*d
    e = quantnorm(w)*(t + 3d0)/(t + 1d0)
    t = sqrt(d*fexpm1(e*e/d))
    z = sign(t, e)
  end function flink_wallace

!! Inverse link fcn (log scale)
  pure function invlink_wallace (z,d) result (w)
    use interfaces, only: logprobnorm, flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision x, t
    t = 8d0*d
    x = sqrt(d*flog1p(z*z/d))*(t + 1d0)/(t + 3d0)
    t = sign(x, z)
    w = logprobnorm(t)
  end function invlink_wallace

  pure function invlinkdz_wallace (z,d) result (w)
    use interfaces, only: logprobnorm, logpdfnorm, flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision t, c, f, u, dzt, zt
    c = 8d0*d
    c = (c + 1d0)/(c + 3d0)
    t = z*z/d
    u = sqrt(d*flog1p(t))
    if (u .eq. 0d0) then
      dzt = c
      zt = 0d0
    else
      dzt = c*abs(z)/u/(1+t)
      zt = c*u
      if (z .lt. 0d0) zt = -zt
    end if
    w = logprobnorm(zt)
    f = logpdfnorm(zt)
    w = exp(f-w)*dzt
  end function invlinkdz_wallace

  pure function invlinkhz_wallace (z,d) result (w)
    use interfaces, only: logprobnorm, logpdfnorm, flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision t, c, f, u, dzt, zt, hzt
    c = 8d0*d
    c = (c + 1d0)/(c + 3d0)
    t = z*z/d
    u = sqrt(d*flog1p(t))
    t = t+1d0
    if (u .eq. 0d0) then
      dzt = c
      zt = 0d0
      hzt = 0d0
    else
      zt = c*u
      dzt = c*abs(z)/u/t
      hzt = t*u
      u = u*u
      hzt = c*((1-z*z/d)*u-z*z)/(t*u*hzt)
      if (z .lt. 0d0) then
        zt = -zt
        hzt = -hzt
      end if
    end if
    w = logprobnorm(zt)
    f = logpdfnorm(zt)
    w = exp(f-w)
    dzt = dzt*dzt
    w = w*(hzt - zt*dzt) - w*w*dzt
  end function invlinkhz_wallace

  pure function invlink3z_wallace (z,d) result (w)
    use interfaces, only: logprobnorm, logpdfnorm, flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision f, w, dw, hw, dzt, hzt, szt, zt
    zt = wallace_zeta(z,d)
    dzt = wallace_zeta_dz(z,d)
    hzt = wallace_zeta_hz(z,d)
    szt = wallace_zeta_3z(z,d)
    dw = invlinkdz_wallace(z,d)
    hw = invlinkhz_wallace(z,d)
    w = logprobnorm(zt)
    f = logpdfnorm(zt)
    w = exp(f-w)
    w = (zt*zt-1)*w*dzt*dzt*dzt -3*zt*w*dzt*hzt + w*szt - dw*dw*dw - 3*dw*hw
  end function invlink3z_wallace

  pure function invlinkdn_wallace (z,d) result (w)
    use interfaces, only: logprobnorm, logpdfnorm, flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision zzd, c, cd, lzzd, u, zt, dzt, phip, phid
    zzd = z*z/d
    c = 8*d
    cd = c + 3
    c = (c+1)/cd
    lzzd = flog1p(zzd)
    u = sqrt(d*lzzd)
    if (u .eq. 0d0) then
      zt = 0d0
      dzt = 0d0
    else
      zt = c*u
      dzt = 16*u/(cd*cd) + .5d0*c*(lzzd - zzd/(1+zzd))/u
      if (z .lt. 0d0) then
        zt = -zt
        dzt = -dzt
      end if
    end if
    phip = logprobnorm(zt)
    phid = logpdfnorm(zt)
    w = exp(phid - phip)*dzt
  end function invlinkdn_wallace

  pure function invlinkhn_wallace (z,d) result (w)
    use interfaces, only: logprobnorm, logpdfnorm, flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision zzd, c, cd, lzzd, u, zt, dzt, hzt, phip, phid, zo
    zzd = z*z/d
    zo = zzd/(1+zzd)
    c = 8*d
    cd = c + 3
    c = (c+1)/cd
    lzzd = flog1p(zzd)
    u = sqrt(d*lzzd)
    if (u .eq. 0d0) then
      zt = 0d0
      dzt = 0d0
      hzt = 0d0
    else
      zt = c*u
      dzt = 16*u/(cd*cd) + .5d0*c*(lzzd - zo)/u
      hzt = -.25d0*c*zo*zo/(u*u*u) &
         + .5d0*((3/d+64*d)*zo/(1+zzd)-32*zo*zo)/(u*cd*cd) &
         - .25d0*(9-72*d+960*d*d+512*d*d*d)*u/(d*d*cd*cd*cd)
      if (z .lt. 0d0) then
        zt = -zt
        dzt = -dzt
        hzt = -hzt
      end if
    end if
    phip = logprobnorm(zt)
    phid = logpdfnorm(zt)
    w = exp(phid-phip)
    dzt = dzt*dzt
    w = w*(hzt - zt*dzt) - w*w*dzt
  end function invlinkhn_wallace

  pure function invlinkdzdn_wallace (z,d) result (w)
    use interfaces, only: logprobnorm, logpdfnorm
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision zt, ztdz, ztdn, ztdzdn, phid, phip
    zt = wallace_zeta(z,d)
    ztdz = wallace_zeta_dz(z,d)
    ztdn = wallace_zeta_dn(z,d)
    ztdzdn = wallace_zeta_dzdn(z,d)
    phid = logpdfnorm(zt)
    phip = logprobnorm(zt)
    w = exp(phid - phip)
    w = w*ztdzdn - (zt+w)*w*ztdz*ztdn
  end function invlinkdzdn_wallace

  pure function invlinkhzdn_wallace (z,d) result (w)
    use interfaces, only: logprobnorm, logpdfnorm
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision zt, ztdz, ztdn, zthz, ztdzdn, zthzdn, phid, phip
    double precision wdz, wdn, whz, wdzdn
    zt = wallace_zeta(z,d)
    ztdz = wallace_zeta_dz(z,d)
    zthz = wallace_zeta_hz(z,d)
    ztdn = wallace_zeta_dn(z,d)
    ztdzdn = wallace_zeta_dzdn(z,d)
    zthzdn = wallace_zeta_hzdn(z,d)
    wdz = invlinkdz_wallace(z,d)
    whz = invlinkhz_wallace(z,d)
    wdn = invlinkdn_wallace(z,d)
    wdzdn = invlinkdzdn_wallace(z,d)
    phid = logpdfnorm(zt)
    phip = logprobnorm(zt)
    w = exp(phid - phip)
    w = (zt*zt-1)*w*ztdz*ztdz*ztdn - 2*zt*w*ztdz*ztdzdn - zt*w*zthz*ztdn &
       + w*zthzdn - wdz*wdz*wdn - 2*wdz*wdzdn - whz*wdn
  end function invlinkhzdn_wallace

  pure function invlinkdzhn_wallace (z,d) result (w)
    use interfaces, only: logprobnorm, logpdfnorm
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision zt, ztdz, ztdn, ztdzdn, zthn, ztdzhn, phid, phip
    double precision dw1, dw2, dw3
    zt = wallace_zeta(z,d)
    ztdz = wallace_zeta_dz(z,d)
    ztdn = wallace_zeta_dn(z,d)
    ztdzdn = wallace_zeta_dzdn(z,d)
    zthn = wallace_zeta_hn(z,d)
    ztdzhn = wallace_zeta_dzhn(z,d)
    phid = logpdfnorm(zt)
    phip = logprobnorm(zt)
    dw1 = exp(phid - phip)
    dw2 = -zt*dw1 - dw1*dw1
    dw3 = -dw1 - zt*dw2 - 2*dw1*dw2
    w = dw3*ztdz*ztdn*ztdn &
       + dw2*(ztdzdn*ztdn*2 + ztdz*zthn) +dw1*ztdzhn
  end function invlinkdzhn_wallace

  pure double precision function wallace_zeta (z,d) result (zt)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision zzd, c, c1, c3, lgz
    c = 8*d
    c1 = 1 + c
    c3 = 3 + c
    zzd = z*z/d
    lgz = flog1p(zzd)
    lgz = d*lgz
    zt = c1/c3*sqrt(lgz)
    if (z .lt. 0d0) then
      zt = -zt
    end if
  end function wallace_zeta

  pure double precision function wallace_zeta_dz (z,d) result (zt)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision zzd, c, c1, c3, lgz
    c = 8*d
    c1 = 1 + c
    c3 = 3 + c
    zzd = z*z/d
    lgz = flog1p(zzd)
    lgz = d*lgz
    zzd = 1 + zzd
    zt = z*c1/(c3*zzd*sqrt(lgz))
    if (z .lt. 0d0) then
      zt = -zt
    end if
  end function wallace_zeta_dz

  pure double precision function wallace_zeta_hz (z,d) result (zt)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision zzd, c, c1, c3, lgz, num, dum
    c = 8*d
    c1 = 1 + c
    c3 = 3 + c
    zzd = z*z/d
    lgz = flog1p(zzd)
    num = c1*(d*lgz - d*zzd*(1+lgz))
    lgz = d*lgz
    zzd = 1 + zzd
    dum = c3*zzd*zzd*lgz*sqrt(lgz)
    zt = num/dum
    if (z .lt. 0d0) then
      zt = -zt
    end if
  end function wallace_zeta_hz

  pure double precision function wallace_zeta_3z (z,d) result (zt)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision zzd, c, lgz, zlgz, num, dum, zsq
    c = 8*d
    zsq = z*z
    zzd = zsq/d
    lgz = flog1p(zzd)
    zlgz = z*lgz
    num = -d*(1+c)*(-3*z*z*z + 3*(d-zsq)*zlgz + 2*(3*d - zsq)*zlgz*lgz)
    zzd = zzd+1
    lgz = d*lgz
    dum = (3+c)*d*zzd*zzd*zzd*lgz*lgz*sqrt(lgz)
    zt = num/dum
    if (z .lt. 0d0) then
      zt = -zt
    end if
  end function wallace_zeta_3z

  pure double precision function wallace_zeta_dn (z,d) result (zt)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision zzd, c, c1, c3, lgz, num, dum
    c = 8*d
    c1 = 1 + c
    c3 = 3 + c
    zzd = z*z/d
    lgz = flog1p(zzd)
    num = -c1*c3*zzd + (3+8*c+c*c)*(1+zzd)*lgz
    lgz = d*lgz
    zzd = 1 + zzd
    dum = 2*c3*c3*zzd*sqrt(lgz)
    zt = num/dum
    if (z .lt. 0d0) then
      zt = -zt
    end if
  end function wallace_zeta_dn

  pure double precision function wallace_zeta_hn (z,d) result (zt)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision zzd, c, c1, c3, cc, lgz, lgzd, lgzr
    c = 8*d
    c1 = 1 + c
    c3 = 3 + c
    cc = c1/c3
    zzd = z*z/d
    lgz = flog1p(zzd)
    zzd = zzd/(1+zzd)
    lgzd = sqrt(d*lgz)
    lgzr = lgz/lgzd
    zt = -.25d0*cc*lgzr*lgzr/lgzd + (8/c3)*(1 - cc)*lgzr &
       - 128/(c3*c3)*(1 - cc)*lgzd &
       - .5d0*zzd*zzd*cc/lgzd*(1d0/d + .5d0/(lgzd*lgzd)) &
       + zzd/lgzd*(.5d0*cc*lgzr/lgzd - 8*(1-cc)/c3)
    if (z .lt. 0d0) then
      zt = -zt
    end if
  end function wallace_zeta_hn

  pure double precision function wallace_zeta_dzdn (z,d) result (zt)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision zzd, c, c1, c3, lgz, num, dum
    c = 8*d
    c1 = 1 + c
    c3 = 3 + c
    zzd = z*z/d
    lgz = flog1p(zzd)
    num = c1*c3*z*zzd - z*(3+c*c)*lgz + z*(3 + 8*c + c*c)*zzd*lgz
    lgz = d*lgz
    zzd = 1 + zzd
    dum = 2*c3*c3*zzd*zzd*lgz*sqrt(lgz)
    zt = num/dum
    if (z .lt. 0d0) then
      zt = -zt
    end if
  end function wallace_zeta_dzdn

  pure double precision function wallace_zeta_hzdn (z,d) result (zt)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision zsq, dsq, zzd, c3, lgz, num, dum
    c3 = 3 + 8*d
    dsq = d*d
    zsq = z*z
    zzd = z*z/d
    lgz = flog1p(zzd)
    num = -3*(3 + 32*d + 64*dsq)*zsq*zsq &
       + 2*zsq*(128*dsq*d - 3*zsq + dsq*(48 - 64*zsq) + d*(6 - 48*zsq))*lgz &
       - (64*dsq*dsq - 384*dsq*d*zsq + 3*zsq*zsq + 2*d*zsq*(-9 + 32*zsq) &
       + dsq*(3 - 192*zsq + 64*zsq*zsq))*lgz*lgz
    lgz = d*lgz
    zzd = 1+zzd
    dum = 2*c3*c3*zzd*zzd*zzd*d*lgz*lgz*sqrt(lgz)
    zt = num/dum
    if (z .lt. 0d0) then
      zt = -zt
    end if
  end function wallace_zeta_hzdn

  pure double precision function wallace_zeta_dzhn (z,d) result (zt)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision zzd, c, c1, c3, cc, lgz, lgzd, lgzd2
    double precision t1, t2, t3
    c = 8*d
    c1 = 1 + c
    c3 = 3 + c
    cc = c1/c3
    zzd = z*z/d
    lgz = flog1p(zzd)
    zzd = zzd/(1+zzd)
    lgzd2 = d*lgz
    lgzd = sqrt(lgzd2)
    t1 = zzd/(z*lgzd)*(.75d0*cc/d - 16/(c3*c3) - 256/(c3*c3*c3)*d)
    t2 = zzd*zzd/(z*lgzd)*(-1.5d0*cc/lgzd2 + 16*d/(c3*c3*lgzd2) &
       - 3*cc/d + 32/(c3*c3))
    t3 = zzd*zzd*zzd*cc/(z*lgzd)*(1.5d0/lgzd2 + 2/d + .75d0*d/(lgzd2*lgzd2))
    zt =  t1 + t2 + t3
    if (z .lt. 0d0) then
      zt = -zt
    end if
  end function wallace_zeta_dzhn


!!!!!!!!!!!!!!!!!!!!!!!!!! modified GEV binomial !!!!!!!!!!!!!!!!!!!!!!!!!!
!! mu = exp{-max(0,1+nu*z)^{-1/nu}} if nu neq 0
!!    = exp{-exp(z)}                if nu == 0
  pure function flink_modgev (w,d) result(z)
    implicit none
    double precision, intent(in) :: w, d
    double precision :: z
    z = -flink_modbc(log(-w),-d)
  end function flink_modgev

  pure function invlink_modgev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = -exp(invlink_modbc(-z,-d))
  end function invlink_modgev

  pure function invlinkdz_modgev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = exp(invlink_modbc(-z,-d))*invlinkdz_modbc(-z,-d)
  end function invlinkdz_modgev

  pure function invlinkhz_modgev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w, y, hz
    y = exp(invlink_modbc(-z,-d))
    w = invlinkdz_modbc(-z,-d)
    hz = invlinkhz_modbc(-z,-d)
    w = -(w*w + hz)*y
  end function invlinkhz_modgev

  pure function invlink3z_modgev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision dz, y, hz, w
    y = exp(invlink_modbc(-z,-d))
    dz = invlinkdz_modbc(-z,-d)
    hz = invlinkhz_modbc(-z,-d)
    w = invlink3z_modbc(-z,-d)
    w = (dz*dz*dz + 3*dz*hz + w)*y
  end function invlink3z_modgev

  pure function invlinkdn_modgev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision dg
    w = exp(invlink_modbc(-z,-d))
    dg = invlinkdn_modbc(-z,-d)
    w = dg*w
  end function invlinkdn_modgev

  pure function invlinkhn_modgev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision dg, d2g
    w = exp(invlink_modbc(-z,-d))
    dg = invlinkdn_modbc(-z,-d)
    d2g = invlinkhn_modbc(-z,-d)
    w = -(dg*dg + d2g)*w
  end function invlinkhn_modgev

  pure function invlinkdzdn_modgev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision dz, dn, dzdn
    w = exp(invlink_modbc(-z,-d))
    dz = invlinkdz_modbc(-z,-d)
    dn = invlinkdn_modbc(-z,-d)
    dzdn = invlinkdzdn_modbc(-z,-d)
    w = -(dz*dn + dzdn)*w
  end function invlinkdzdn_modgev

  pure function invlinkhzdn_modgev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision dz, dn, hz, dzdn, hzdn
    w = exp(invlink_modbc(-z,-d))
    dz = invlinkdz_modbc(-z,-d)
    hz = invlinkhz_modbc(-z,-d)
    dn = invlinkdn_modbc(-z,-d)
    dzdn = invlinkdzdn_modbc(-z,-d)
    hzdn = invlinkhzdn_modbc(-z,-d)
    w = (dz*dz*dn + 2*dz*dzdn + hz*dn +hzdn)*w
  end function invlinkhzdn_modgev

  pure function invlinkdzhn_modgev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision dz, dn, hn, dzdn, dzhn
    w = exp(invlink_modbc(-z,-d))
    dz = invlinkdz_modbc(-z,-d)
    dn = invlinkdn_modbc(-z,-d)
    hn = invlinkhn_modbc(-z,-d)
    dzdn = invlinkdzdn_modbc(-z,-d)
    dzhn = invlinkdzhn_modbc(-z,-d)
    w = (dz*dn*dn + dz*hn + 2*dzdn*dn + dzhn)*w
  end function invlinkdzhn_modgev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GEV binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mu = exp{-max(0,1+nu*z)^{-1/nu}} if nu neq 0
!!    = exp{-exp(z)}                if nu == 0
  pure function flink_gev (w,d) result(z)
    implicit none
    double precision, intent(in) :: w, d
    double precision :: z
    z = -flink_boxcox(log(-w),-d)
  end function flink_gev

  pure function invlink_gev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = -exp(invlink_boxcox(-z,-d))
  end function invlink_gev

  pure function invlinkdz_gev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = exp(invlink_boxcox(-z,-d))*invlinkdz_boxcox(-z,-d)
  end function invlinkdz_gev

  pure function invlinkhz_gev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w, y, hz
    y = exp(invlink_boxcox(-z,-d))
    w = invlinkdz_boxcox(-z,-d)
    hz = invlinkhz_boxcox(-z,-d)
    w = -(w*w + hz)*y
  end function invlinkhz_gev

  pure function invlink3z_gev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision dz, y, hz, w
    y = exp(invlink_boxcox(-z,-d))
    dz = invlinkdz_boxcox(-z,-d)
    hz = invlinkhz_boxcox(-z,-d)
    w = invlink3z_boxcox(-z,-d)
    w = (dz*dz*dz + 3*dz*hz + w)*y
  end function invlink3z_gev

  pure function invlinkdn_gev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision dg
    w = exp(invlink_boxcox(-z,-d))
    dg = invlinkdn_boxcox(-z,-d)
    w = dg*w
  end function invlinkdn_gev

  pure function invlinkhn_gev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision dg, d2g
    w = exp(invlink_boxcox(-z,-d))
    dg = invlinkdn_boxcox(-z,-d)
    d2g = invlinkhn_boxcox(-z,-d)
    w = -(dg*dg + d2g)*w
  end function invlinkhn_gev

  pure function invlinkdzdn_gev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision dz, dn, dzdn
    w = exp(invlink_boxcox(-z,-d))
    dz = invlinkdz_boxcox(-z,-d)
    dn = invlinkdn_boxcox(-z,-d)
    dzdn = invlinkdzdn_boxcox(-z,-d)
    w = -(dz*dn + dzdn)*w
  end function invlinkdzdn_gev

  pure function invlinkhzdn_gev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision dz, dn, hz, dzdn, hzdn
    w = exp(invlink_boxcox(-z,-d))
    dz = invlinkdz_boxcox(-z,-d)
    hz = invlinkhz_boxcox(-z,-d)
    dn = invlinkdn_boxcox(-z,-d)
    dzdn = invlinkdzdn_boxcox(-z,-d)
    hzdn = invlinkhzdn_boxcox(-z,-d)
    w = (dz*dz*dn + 2*dz*dzdn + hz*dn +hzdn)*w
  end function invlinkhzdn_gev

  pure function invlinkdzhn_gev (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision dz, dn, hn, dzdn, dzhn
    w = exp(invlink_boxcox(-z,-d))
    dz = invlinkdz_boxcox(-z,-d)
    dn = invlinkdn_boxcox(-z,-d)
    hn = invlinkhn_boxcox(-z,-d)
    dzdn = invlinkdzdn_boxcox(-z,-d)
    dzhn = invlinkdzhn_boxcox(-z,-d)
    w = (dz*dn*dn + dz*hn + 2*dzdn*dn + dzhn)*w
  end function invlinkdzhn_gev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Other functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function identity (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = z
  end function identity

  pure function constant (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = 1d0
  end function constant

  pure function zero (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = 0d0
  end function zero

end module modelfcns_link

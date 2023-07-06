!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Derivatives of inverse link fcn w.r.t. first arg.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module linkdz

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 1st derivative w.r.t. z
  pure function invlinkdz_bi (z, d) result (y)
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
      y = 1d0/(1d0 + exp(z))
    else
      ff = logprobnorm(z)
      f = logpdfnorm(z)
      y = exp(f - ff)      
    end if
  end function invlinkdz_bi

!! 2nd derivative w.r.t. z
  pure function invlinkddz_bi (z, d) result (y)
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
      y = z*.5d0
      y = .5d0/cosh(y)
      y = -y*y
    else
      ff = logprobnorm(z)
      f = logpdfnorm(z)
      y = f/ff
      y = -y*((d+1d0)*z/(d+z*z) + y)
    end if
  end function invlinkddz_bi

!! 3rd derivative w.r.t. z
  pure function invlinkdddz_bi (z, d) result (y)
    use interfaces, only: logprobt, logprobnorm, logpdft, logpdfnorm
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    double precision f, ff
    if (d .gt. 0d0) then
      f = invlinkdz_bi(z, d)
      ff = invlinkddz_bi(z, d)
      y = d + z*z
      y = (d+1d0)*(d + d - y)/(y*y)*f + (d+1d0)*z/y*ff + 2d0*f*ff
      y = -y
    else if (d .lt. 0d0) then
      y = z
      ff = .5d0*y
      y = 1d0/(sinh(y))
      ff = sinh(ff)
      ff = ff*ff
      y = y*y*y*ff*ff
      y = y + y
    else
      f = invlinkdz_bi(z, d)
      ff = invlinkddz_bi(z, d)
      y = f + z*ff + 2d0*f*ff
      y = -y
    end if
  end function invlinkdddz_bi
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GEV Binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure function invlinkdz_ba (z,d) result(w)
    ! Binomial asymmetric (GEV) link fcn
    ! mu = 1 - exp{-max(0,1+d*z)**(1/d)} , if d /= 0
    !    = 1 - exp(-exp(z))              , if d == 0
    ! w = log(1-mu)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = -exp(z)
    else if (d .eq. 0.5d0) then
      w = 1d0 + 0.5d0*z
      if (w .gt. 0d0) then
        w = -w
      else
        w = 0d0
      end if
    else if (d .eq. -0.5d0) then
      w = 1d0 - 0.5d0*z
      if (w .gt. 0d0) then
        w = -1d0/(w*w*w)
      else
        w = 0d0
      end if
    else if (d .eq. 1d0) then
      w = 1d0 + z
      if (w .gt. 0d0) then
        w = -1d0
      else
        w = 0d0
      end if
    else if (d .eq. -1d0) then
      w = 1d0 - z
      if (w .gt. 0d0) then
        w = -1d0/(w*w)
      else
        w = 0d0
      end if
    else if (d .eq. 2d0) then
      w = 1d0 + z + z
      if (w .gt. 0d0) then
        w = -1d0/sqrt(w)
      else
        w = 0d0
      end if
    else if (d .eq. -2d0) then
      w = 1d0 - z - z
      if (w .gt. 0d0) then
        w = -1d0/(w*sqrt(w))
      else
        w = 0d0
      end if
    else
      w = 1d0 + d*z
      if (w .gt. 0d0) then
        w = -w**(1d0/d - 1d0)
      else
        w = 0d0
      end if
    end if
  end function invlinkdz_ba

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GEVD Binomial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure elemental function invlinkdz_bd (z, d) result (w)
    ! Binomial asymmetric (GEV) decreasing link fcn
    ! p = 1 - exp{-max(0, 1+d*z)**(-1/d)}
    ! w = log(1-p) = -max(0, 1+d*z)**(-1/d)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = exp(-z)
    else
      w = 1d0 + d*z
      if (w .gt. 0d0) then
        if (d .eq. -1d0) then
          w = 1d0
        else
          w = w**(-1d0/d - 1d0)
        end if
      else
        w = 0d0
      end if
    end if
  end function invlinkdz_bd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Box-Cox !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure elemental function invlinkdz_boxcox (z,d) result (w)
    ! w = log(mu)
    use interfaces, only: flog1p
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Poisson !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure elemental function invlinkdz_po (z,d) result (w)
    ! w = log(mu)
    use interfaces, only: flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    if (d .eq. 0d0) then
      w = 1d0
    else if (d .gt. 0d0) then
      w = 1d0 + d*abs(z)
      w = 1d0/w
    else ! Use regular Box-Cox
      w = d*z
      if (w .gt. -1d0) then
        w = 1d0/(1d0+w)
      else
        w = 0d0
      end if
    end if
  end function invlinkdz_po

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gamma !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure elemental function invlinkdz_gm (z,d) result (w)
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    w = invlinkdz_po(z, d)
  end function invlinkdz_gm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gaussian !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure elemental function invlinkdz_ga (z,d) result(w)
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

!!!!!!!!!!!!!!!!!!! Wallace transformation for binomial !!!!!!!!!!!!!!!!!!!
  pure elemental function invlinkdz_bw (z,d) result (w)
    use interfaces, only: logprobnorm, logpdfnorm, flog1p
    implicit none
    double precision, intent(in) :: z, d
    double precision w
    double precision x, t, c, f, u
    c = 8d0*d
    c = (c + 1d0)/(c + 3d0)
    t = z*z/d
    u = sqrt(d*flog1p(t))
    u = sign(u,z)
    x = c*u
    w = logprobnorm(x)
    f = logpdfnorm(x)
    w = exp(f-w)*c*z/((1d0+t)*u)
  end function invlinkdz_bw
end module linkdz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module linkfcn
  implicit none
  double precision, parameter :: bigpos = huge(1d0), bigneg = -bigpos
  private bigpos, bigneg
contains
  function invlink_bi (z,d) result (y)
    use interfaces, only: logprob
    implicit none
    double precision, intent(in) :: z, d
    double precision y
    y = logprob(z,d)
  end function invlink_bi

  function invlink_po (z,d) result (y)
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

  function invlink_gm (z,d) result (y)
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
    y = -exp(-y)
  end function invlink_gm

  function invlink_ga (z,d) result(y)
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
end module linkfcn

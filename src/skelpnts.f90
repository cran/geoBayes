!!! skelpnts.f90 --- 
!! 
!! Author: Evangelos Evangelou
!! Created: Thu, 25 Sep, 2014 16:19 (BST)
!! Last-Updated: Wed, 29 Oct, 2014 16:15 (GMT)
!!     Update #: 33
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: These routines are used for finding appropriate skeleton
!!! points. 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine jointyzfc_bi (fc, gr, hm, z, y1, y2, Ups, nu, ssq, &
   tsq, xi, n)
  ! fc: log-likelihood
  ! gr: score
  ! hm: Fisher info matrix
  use interfaces, only: fexpm1, flog1mexp
  use linkfcn
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: z(n), y1(n), y2(n), Ups(n, n), nu, ssq, &
     tsq, xi(n)
  double precision, intent(out) :: fc, gr(n), hm(n, n)
  double precision p0, p1, p2, zmxi(n), tmp
  integer i

  !! Gaussian likelihood
  zmxi = z - xi
  call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,gr,1) ! gr = Ups*(z-xi)
  gr = gr/ssq
  fc = -0.5d0*dot_product(zmxi,gr) ! WARNING No determinant term

  !! Exponential likelihood
  do i = 1, n
    p0 = invlink_bi(z(i), nu)
    p1 = invlink1_bi(z(i), nu)
    p2 = invlink2_bi(z(i), nu)
    tmp = 1d0/fexpm1(-p0)
    fc = fc + (y1(i)*p0 + y2(i)*flog1mexp(p0))/tsq
    gr(i) = gr(i) + p1*(y1(i) - y2(i)*tmp)/tsq
    tmp = p2*(y1(i) - y2(i)*tmp)/tsq &
       - p1*p1*y2(i)*tmp*(1d0 + tmp)/tsq
    hm(:i, i) = Ups(:i, i)/ssq
    hm(i, i) = hm(i, i) + tmp
  end do
end subroutine jointyzfc_bi



subroutine gaussaprx (mean, prec, fcyz, y1, y2, Ups, nu, xi, ssq, &
   tsq, n, ifam)
  ! Gaussian approximation to f(z|y;phi,theta)
  ! On input, mean is the starting value for z
  implicit none
  integer, intent(in) :: n, ifam
  double precision, intent(in) :: y1(n), y2(n), Ups(n, n), nu, xi(n), ssq, &
     tsq 
  double precision, intent(inout) :: mean(n)
  double precision, intent(out) :: prec(n, n), fcyz
  double precision gr(n)
  integer opnb(n), opipr, iflag, maxiter, i
  double precision oplb(n), opub(n), opfac, oppgt
  parameter ( maxiter = 200, opfac = 1d7, oppgt = 0d0 )

  opnb = 0
  iflag = 0
  select case (ifam)
  case (2) ! Binomial
    do i = 1, maxiter
      call jointyzfc_bi (fcyz, gr, prec, mean, y1, y2, Ups, nu, ssq, tsq, &
         xi, n) 
      fcyz = -fcyz
      gr = -gr
      call lbfgsbmod (n, mean, oplb, opub, opnb, fcyz, gr, opipr, &
         opfac, oppgt, iflag)
      if (iflag .eq. 0) then
        exit
      else if (iflag .lt. 0) then
        call rwarn ('The optimisation for the Gaussian approximation didn''&
           &t converge')
        exit
      end if
    end do
  end select
  if (iflag .gt. 0) then
    call rwarn ('The optimisation for the Gaussian approximation needs more&
       & iterations')
  end if
  fcyz = -fcyz
  ! TODO Where to put the log determinant of precision?
end subroutine gaussaprx



subroutine poster (fff, y1, y2, Ups, nu, xi, ssq, &
   tsq, n, ifam)
  ! Gives the approximation to the log posterior of the nuissance
  ! parameters f(psi|y,theta) propto f(y,z,psi,theta)/f(z|y,psi,theta);
  ! psi = ssq here
  implicit none
  integer, intent(in) :: n, ifam
  double precision, intent(in) :: y1(n), y2(n), Ups(n, n), nu, xi(n), ssq, &
     tsq
  double precision :: precision(n, n), fcyz
  double precision z(n)
  select case (ifam)
  case (2) ! Binomial
    z = log((y1 + .5d0)/(y2 + .5d0)) ! TODO Better StVal using link
    call jointyzfc_bi (fcyz, gr, hm, z, y1, y2, Ups, nu, ssq, tsq, xi, n)
  end select
  
end subroutine poster




! Local Variables:
! compile-command: "gfortran -c -fpic -Wunused-parameter -Wall \
!   -pedantic -o skelpnts.o skelpnts.f90"
! End:

!!! skelpnts.f90 ---
!!
!! Author: Evangelos Evangelou
!! Created: Thu, 25 Sep, 2014 16:19 (BST)
!! Last-Updated: Wed, 1 Jun, 2016 13:26 (BST)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Commentary: These routines are used for finding appropriate skeleton
!!! points.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gaussaprx (mean, prec, fcyz, y1, y2, Ups, ldh_Ups, &
   nu, xi, lmxi, ssq, tsq, tsqdf, n, ifam)
  ! Gaussian approximation to f(z|y;phi,ssq)
  ! mean is the Gaussian mean of the approximation
  ! prec is the precision matrix
  ! fcyz is the value at the maximum
  ! Rest are input to the joint pdf
  use lbfgsbmod
  use pdfdz, only: logcondyzdz_gt, logcondyzhs_gt
  use linkfcns, only: flink_ga
  use modelfcns, only: mustart, flink_sp => flink, &
     logpdfzdz, logcondyzdz_sp => logcondyzdz, logcondyzhs_sp => logcondyzhs
  implicit none
  integer, intent(in) :: n, ifam
  logical, intent(in) :: lmxi
  double precision, intent(in) :: y1(n), y2(n), Ups(n, n), ldh_Ups, &
     nu, xi(n), ssq, tsq, tsqdf
  double precision, intent(out) :: mean(n), prec(n, n), fcyz
  double precision gr(n), gr1(n), gr2(n), fc1, fc2, hs(n)
  double precision tsqdfsc, respdf
  integer opnb(n), opipr, iflag, maxiter, i
  double precision oplb(n), opub(n), opfac, oppgt
  parameter ( maxiter = 200, opfac = 1d7, oppgt = 0d0 )

  opnb = 0
  iflag = 0
  select case (ifam)
  case (0) ! Gaussian transformed
    tsqdfsc = tsq*tsqdf
    respdf = n + tsqdf
    mean = flink_ga(y1,nu) ! Starting value
    do i = 1, maxiter
      call logpdfzdz (fc1, gr1, mean, Ups, ldh_Ups, xi, lmxi, ssq, n)
      call logcondyzdz_gt (fc2, gr2, nu, y1, y2, mean, n, tsqdfsc)
      fcyz = -fc1 - respdf*fc2
      gr = -gr1 - respdf*gr2
      call lbfgsb (n, mean, oplb, opub, opnb, fcyz, gr, opipr, &
         opfac, oppgt, iflag)
      if (iflag .eq. 0) then
        exit
      else if (iflag .lt. 0) then
        call rwarn ('The optimisation for the Gaussian approximation didn''&
           &t converge')
        exit
      end if
    end do
    if (iflag .gt. 0) then
      call rwarn ('The optimisation for the Gaussian approximation needs more&
         & iterations')
    end if
    call logcondyzhs_gt (hs, nu, y1, y2, mean, n, tsqdfsc)
    hs = respdf*hs
  case default
    gr = mustart(y1,y2)
    mean = flink_sp(gr,nu) ! Starting value
    do i = 1, maxiter
      call logpdfzdz (fc1, gr1, mean, Ups, ldh_Ups, xi, lmxi, ssq, n)
      call logcondyzdz_sp (fc2, gr2, nu, y1, y2, mean, n, tsq)
      fcyz = -fc1 - fc2
      gr = -gr1 - gr2
      call lbfgsb (n, mean, oplb, opub, opnb, fcyz, gr, opipr, &
         opfac, oppgt, iflag)
      if (iflag .eq. 0) then
        exit
      else if (iflag .lt. 0) then
        call rwarn ('The optimisation for the Gaussian approximation didn''&
           &t converge')
        exit
      end if
    end do
    if (iflag .gt. 0) then
      call rwarn ('The optimisation for the Gaussian approximation needs more&
         & iterations')
    end if
    call logcondyzhs_sp (hs, nu, y1, y2, mean, n, tsq)
  end select
  fcyz = -fcyz
  do i = 1, n
    prec(:i,i) = Ups(:i,i)/ssq ! Upper triangular elements
    prec(i,i) = prec(i,i) - hs(i)
  end do
end subroutine gaussaprx


subroutine poster (fssq, meang, prechg, ssq, ssqdfh, ssqdfsc, &
   y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  ! Gives the approximation to the log posterior of the nuissance
  ! parameters f(psi|y,theta) propto f(y,z,psi,theta)/f(z|y,psi,theta);
  ! psi = ssq here
  ! meang is the mean of the Gaussian approximation for that ssq
  ! prechg is the Chol of the precision matrix of the Gaussian
  ! approximation for that ssq.
  implicit none
  integer, intent(in) :: n, ifam
  logical, intent(in) :: lmxi
  double precision, intent(in) :: y1(n), y2(n), Ups(n, n), ldh_Ups, nu, &
     xi(n), ssq, tsq, tsqdf, ssqdfh, ssqdfsc
  double precision, intent(out) :: fssq, meang(n), prechg(n,n)
  double precision :: ldh, fpr
  logical, external :: disnan ! LAPACK auxiliary routine.
  if (disnan(ssq)) call rexit("poster - ssq entered is NaN.")
  fpr = -(1d0+ssqdfh)*log(ssq) - .5d0*ssqdfsc/ssq ! Prior
  call gaussaprx (meang, prechg, fssq, y1, y2, Ups, ldh_Ups, &
     nu, xi, lmxi, ssq, tsq, tsqdf, n, ifam)
  call factorpd(n,prechg,ldh)
  fssq = fssq - ldh + fpr
contains
  subroutine factorpd(n,x,ldh)
    integer, intent(in) :: n
    double precision, intent(inout) :: x(n,n)
    double precision, intent(out) :: ldh
    integer i
    call dpotrf ('u',n,x,n,i)
    if (i .ne. 0) then
      call rexit("poster - Non positive definite matrix")
    end if
    ldh = 0d0
    do i = 1,n
      ldh = ldh + log(x(i,i))
    end do
  end subroutine factorpd
end subroutine poster


subroutine posterlog (fval, meang, prechg, logssq, ssqdfh, ssqdfsc, &
   y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  ! Gives the approximation to the log posterior of the nuissance
  ! parameters f(psi|y,theta) propto f(y,z,psi,theta)/f(z|y,psi,theta);
  ! psi = log(ssq) here
  ! meang is the mean of the Gaussian approximation for that ssq
  ! prechg is the Chol of the precision matrix of the Gaussian
  ! approximation for that psi.
  implicit none
  integer, intent(in) :: n, ifam
  logical, intent(in) :: lmxi
  double precision, intent(in) :: y1(n), y2(n), Ups(n, n), ldh_Ups, nu, &
     xi(n), logssq, tsq, tsqdf, ssqdfh, ssqdfsc
  double precision, intent(out) :: fval, meang(n), prechg(n,n)
  double precision ssq
  logical, external :: disnan
  if (disnan(logssq)) call rexit("posterlog - logssq entered is NaN.")
  ssq = exp(logssq)
  call poster (fval, meang, prechg, ssq, ssqdfh, ssqdfsc, &
     y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  fval = fval + logssq
end subroutine posterlog

subroutine optlogssq (tval, tprc, pdfval, meang, prechg, ssqdfh, ssqdfsc, &
   y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  ! Maximise the posterior for t = log(ssq). Uses 2nd degree polynomial
  ! interpolation to the posterior of log(ssq).
  implicit none
  integer, intent(in) :: n, ifam
  logical, intent(in) :: lmxi
  double precision, intent(in) :: y1(n), y2(n), Ups(n, n), ldh_Ups, nu, &
     xi(n), tsq, tsqdf, ssqdfh, ssqdfsc
  double precision, intent(inout) :: tval
  double precision, intent(out) :: tprc, pdfval, meang(n), prechg(n,n)
  double precision pdfpnts(3), tnew, tpnts(3), mtmp(n), ptmp(n,n), pdfnew
  double precision, parameter :: eps = 1d-3
  integer, parameter :: maxiter = 100
  integer i, ii
  logical lconv, ok
  tpnts(1) = .8*tval - 15d-1
  tpnts(2) = tval
  tpnts(3) = 1.2*tval + 15d-1
  do i = 1, 3
    call posterlog (pdfpnts(i), mtmp, ptmp, tpnts(i), ssqdfh, ssqdfsc, &
       y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
    if (i .eq. 1 .or. pdfpnts(i) .gt. pdfval) then
      tval = tpnts(i)
      pdfval = pdfpnts(i)
      meang = mtmp
      prechg = ptmp
    end if
  end do
  !!print*, 0, real(tpnts), real(pdfpnts), real(pdfval)
  do i = 1, maxiter
    call polmax(tnew, ii, ok, tpnts, pdfpnts) ! Find new point
    call posterlog (pdfnew, mtmp, ptmp, tnew, ssqdfh, ssqdfsc, &
       y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
    lconv = .false.
    if (ok) lconv = abs(pdfnew - pdfval) .lt. eps !eps*(1d0+abs(pdfval))
    if (ok .and. .not. lconv) lconv = minval(abs(tnew-tpnts)) .lt. 15d-1
    if (lconv) then
      exit
    end if
    tpnts(ii) = tnew        ! Replace with new value
    tval = tnew             ! Current best guess
    pdfpnts(ii) = pdfnew
    if (pdfnew .gt. pdfval) then
      pdfval = pdfnew
      meang = mtmp
      prechg = ptmp
    end if
  !!print*, i, real(tpnts), real(pdfpnts), real(pdfval)
  end do
  if (.not. lconv) then
    call rwarn("optlogssq - Was not able to find the maximum posterior.")
  end if
  tprc = -.5d0*polcf2(tpnts, pdfpnts)
  if (tprc .le. 0d0) then
    call rexit("optlogssq - Computed non-positive precision.")
  end if
contains
  subroutine polmax(xmax, irem, ok, x, y)
    ! Find the argmax of a 2nd degree interpolated polynomial.
    implicit none
    double precision, intent(in), dimension(3) :: x, y
    double precision, intent(out) :: xmax
    integer, intent(out) :: irem
    logical, intent(out) :: ok
    double precision :: xm(3), xp(3), a, b
    xm(1) = x(2) - x(3)
    xm(2) = x(3) - x(1)
    xm(3) = x(1) - x(2)
    xp(1) = x(2) + x(3)
    xp(2) = x(3) + x(1)
    xp(3) = x(1) + x(2)
    b = xm(1)*xm(2)*xm(3)
    xm = y*xm
    a = xm(1) + xm(2) + xm(3)
    b = a*b
    ok = b .gt. 0d0 ! Check if the polynomial has a max
    if (ok) then
      xp = xm*xp
      b = xp(1) + xp(2) + xp(3)
      xmax = .5d0*b/a
      irem = minloc(y,1)
    else ! polynomial is the otherway around
      xmax = (x(1) + x(2) + x(3))/3d0
      xm = abs(x - xmax)
      irem = minloc(xm,1)
      if (xm(irem) .eq. 0d0) then
        irem = maxloc(xm,1)
        if (xm(irem) .eq. 0d0) then
          call rexit ("polmax - Cannot locate max")
        end if
      end if
    end if
  end subroutine polmax

  double precision function polcf2(x, y)
    ! Find the coefficient of x^2 on the interolated polynomial.
    implicit none
    double precision, intent(in), dimension(3) :: x, y
    double precision :: xm(3), xp(3)
    xm(1) = x(2) - x(3)
    xm(2) = x(3) - x(1)
    xm(3) = x(1) - x(2)
    xp(1) = xm(2)*xm(3)
    xp(2) = xm(3)*xm(1)
    xp(3) = xm(1)*xm(2)
    xm = y/xp
    polcf2 = -sum(xm)
  end function polcf2

  double precision function variance(x,n)
    integer, intent(in) :: n
    double precision, intent(in) :: x(1:n)
    integer i
    double precision mn, df
    mn = x(1)
    variance = 0d0
    do i = 2, n
      df = x(i) - mn
      mn = mn + df/i
      variance = variance + df*(x(i) - mn)
    end do
    variance = variance/(n-1)
  end function variance
end subroutine optlogssq

subroutine gridposter (np, tg, twght, meang, prechg, ssqdfh, ssqdfsc, &
   ssqin, y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  ! Create an appropriate grid of size np for the posterior of ssq. The
  ! grid size is 2*np+1 and the midpoint corresponds to the maximiser of
  ! the approximate posterior.
  implicit none
  integer, intent(in) :: np, n, ifam
  logical, intent(in) :: lmxi
  double precision, intent(in) :: y1(n), y2(n), Ups(n, n), ldh_Ups, nu, &
     xi(n), tsq, tsqdf, ssqdfh, ssqdfsc, ssqin
  double precision, intent(out) :: tg(np+np+1), twght(np+np+1), &
     meang(n,np+np+1), prechg(n,n,np+np+1)
  double precision tmx, tsd, tlo, thi, step, sfctr, tst
  integer i, m, e
  double precision, parameter :: zq = 4d0, eps = -6.5
  m = np + 1
  e = m + np
  sfctr = 1d0 - 1d0/np
  ! Find the middle point
  tmx = log(ssqin)
  call optlogssq (tmx, tsd, twght(m), meang(:,m), prechg(:,:,m), ssqdfh, &
     ssqdfsc, y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  tg(m) = tmx
  tsd = zq/sqrt(tsd)
  ! Lo
  tlo = tmx - tsd
  step = tsd/np
  tg(1) = tlo
  do i = 1, 20
    call posterlog (twght(1), meang(:,1), prechg(:,:,1), tg(1), ssqdfh, &
       ssqdfsc, y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
    tst = twght(1) - twght(m)
    if (tst .gt. eps) exit
    tlo = tlo + step
    tg(1) = tlo
    step = step*sfctr
  end do
  do i = 2, np
    tg(i) = tg(i-1) + step
    call posterlog (twght(i), meang(:,i), prechg(:,:,i), tg(i), ssqdfh, &
       ssqdfsc, y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
!     twght(i) = twght(i) - twght(m)
  end do
  ! Hi
  thi = tmx + tsd
  step = tsd/np
  tg(e) = thi
  do i = 1, 20
    call posterlog (twght(e), meang(:,e), prechg(:,:,e), tg(e), ssqdfh, &
       ssqdfsc, y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
    tst = twght(e) - twght(m)
    if (tst .gt. eps) exit
    thi = thi - step
    tg(e) = thi
    step = step*sfctr
  end do
  do i = e-1, np+2, -1
    tg(i) = tg(i+1) - step
    call posterlog (twght(i), meang(:,i), prechg(:,:,i), tg(i), ssqdfh, &
       ssqdfsc, y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
!     twght(i) = twght(i) - twght(m)
  end do
!   twght(m) = 0d0
  !! twght = twght - log(sum(exp(twght))) ! XXX Not normalised.
end subroutine gridposter


subroutine llikpars2 (fval, nu, phi, nsq, kappa, y1, y2, F, betm0, betQ0, &
   ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, ssqin, ifam, icf)
  use covfun, only: create_spcor, calc_cov
  use betaprior
  implicit none
  integer, intent(in) :: np, n, p, ifam, icf
  double precision, intent(in) :: ssqin, ssqdf, ssqsc, &
     y1(n), y2(n), phi, nu, nsq, kappa, dm(n,n), F(n,p), &
     betm0(p), betQ0(p,p), tsq, tsqdf
  double precision, intent(out) :: fval
  double precision Ups(n,n), ldh_Ups, xi(n), modeldfh, ssqdfh, ssqdfsc
  integer i
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p)
  double precision tg(np+np+1), logw(np+np+1), &
     meang(n,np+np+1), prechg(n,n,np+np+1)
  double precision ssqst
  call create_spcor(icf,n)
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)
  call calc_cov (phi,nsq,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
  ssqdfh = .5d0*ssqdf
  ssqdfsc = ssqdf*ssqsc
  ssqst = ssqstart(y1,y2,nu,Ups,n,ifam)
  ! TODO ssqin or ssqst?
  ! ssqst = ssqin
  call gridposter (np, tg, logw, meang, prechg, ssqdfh, ssqdfsc, &
     ssqst, y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  i = np+np+1
  fval = trapezoid(i, tg, logw)
contains
  function trapezoid (n, x, f)
    ! Compute the log integral of exp(f) along x
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), f(n)
    double precision trapezoid
    double precision fpf(n-1), xmx(n-1), fmx, fmm(n), ln2
    parameter (ln2 = 0.6931471805599453d0)
    fmx = maxval(f)
    fmm = f - fmx
    fmm = exp(fmm)
    fpf = fmm(2:n) + fmm(1:n-1)
    xmx = x(2:n) - x(1:n-1)
    trapezoid = dot_product(xmx,fpf)
    trapezoid = log(trapezoid) + fmx - ln2
  end function trapezoid

  function ssqstart (y1,y2,nu,Ups,n,ifam)
    use linkfcns, only: flink_ga
    use modelfcns, only: mustart, flink_sp => flink
    implicit none
    integer, intent(in) :: n, ifam
    double precision, intent(in) :: y1(n), y2(n), nu, Ups(n,n)
    double precision ssqstart
    double precision, dimension(n) :: z, tmp
    double precision zbar
    if (ifam .eq. 0) then ! Gaussian transformed
      z = flink_ga(y1,nu) ! Starting value
    else
      tmp = mustart(y1,y2)
      z = flink_sp(tmp,nu) ! Starting value
    end if
    zbar = sum(z)/n
    z = z - zbar
    call dsymv("u",n,1d0,Ups,n,z,1,0d0,tmp,1)
    ssqstart = dot_product(z,tmp)
    ssqstart = ssqstart/n
  end function ssqstart
end subroutine llikpars2


subroutine llikparsval (fval, nu, phi, nsq, kappa, y1, y2, F, betm0, betQ0, &
   ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, ssqin, ifam, icf)
  use modelfcns, only: create_model
  implicit none
  integer, intent(in) :: np, n, p, ifam, icf
  double precision, intent(in) :: ssqin, ssqdf, ssqsc, &
     y1(n), y2(n), phi, nu, nsq, kappa, dm(n,n), F(n,p), &
     betm0(p), betQ0(p,p), tsq, tsqdf
  double precision, intent(out) :: fval
  call create_model (ifam)
  call llikpars2 (fval, nu, phi, nsq, kappa, y1, y2, F, betm0, betQ0, &
   ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, ssqin, ifam, icf)
end subroutine llikparsval

subroutine llikparscalc (fval, nu, phi, nsq, kappa, npars, &
   y1, y2, F, betm0, betQ0, &
   ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, ssqin, ifam, icf)
  ! Calls subroutine llikpars2 repeatedly for each parameter input.
  use modelfcns
  implicit none
  integer, intent(in) :: np, n, p, ifam, icf, npars
  double precision, intent(in) :: ssqin, ssqdf, ssqsc, &
     y1(n), y2(n), phi(npars), nu(npars), nsq(npars), kappa(npars), &
     dm(n,n), F(n,p), betm0(p), betQ0(p,p), tsq, tsqdf
  double precision, intent(out) :: fval(npars)
  integer i
  call create_model (ifam)
  do i = 1, npars
    call llikpars2 (fval(i), nu(i), phi(i), nsq(i), kappa(i), &
       y1, y2, F, betm0, betQ0, &
       ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, ssqin, ifam, icf)
  end do
end subroutine llikparscalc


! Local Variables:
! compile-command: "gfortran -c -fpic -Wunused-parameter -Wall \
!   -pedantic -o skelpnts.o skelpnts.f90"
! End:

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
  use linkfcns, only: flink_ga
  use modelfcns, only: mustart, flink_sp => flink
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
  parameter ( maxiter = 1500, opfac = 1d7, oppgt = 0d0 )

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
  fcyz = -fcyz ! Because it was minus before
  do i = 1, n
    prec(:i,i) = Ups(:i,i)/ssq ! Upper triangular elements
    prec(i,i) = prec(i,i) - hs(i)
  end do

contains

  subroutine logpdfzdz (fc, gr, z, Ups, ldh_Ups, xi, lmxi, ssq, n)
!! log-pdf of z and its derivative after integrating out beta. The 2*pi
!! constant is removed.
    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: lmxi
    double precision, intent(in) :: z(n), Ups(n,n), ldh_Ups, xi(n), ssq
    double precision, intent(out) :: fc, gr(n)
    double precision zmxi(n)
    if (lmxi) then
      zmxi = z - xi
    else
      zmxi = z
    end if
    call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,gr,1) ! gr = Ups*(z-xi)
    gr = -gr/ssq ! gr = -Ups*(z-x)/ssq
    fc = -.5d0*n*log(ssq) + ldh_Ups + .5d0*dot_product(zmxi,gr)
  end subroutine logpdfzdz

  subroutine logcondyzdz_gt (fc, gr, nu, y1, y2, z, n, tsq)
!! TODO respdf is not needed here. Can be computed outside the fcn. XXX
!! log-pdf of y|z and its derivative w.r.t. z.
!! Gaussian transformed
!! tsq is tsqdf*tsqsc on input
!! respdf is n + tsqdf on input
    use linkfcns, only: invlink_ga
    use linkdz, only: invlinkdz_ga
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz
    fc = tsq
    do i = 1, n
      par = invlink_ga(z(i),nu)
      pardz = invlinkdz_ga(z(i),nu)
      par = y1(i) - par
      gr(i) = y2(i)*par
      fc = fc + gr(i)*par
      gr(i) = gr(i)*pardz
    end do
    gr = gr/fc ! XXX Missing respdf factor. Must be computed outside the fcn
    fc = -.5d0*log(fc) ! XXX Also missing respdf = n + tsqdf
  end subroutine logcondyzdz_gt

  subroutine logcondyzhs_gt (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
!! Gaussian transformed
!! tsq is tsqdf*tsqsc on input
    use linkfcns, only: invlink_ga
    use linkdz, only: invlinkdz_ga
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision par, pardz, fc
    fc = tsq
    do i = 1, n
      par = invlink_ga(z(i),nu)
      pardz = invlinkdz_ga(z(i),nu)
      hs(i) = y2(i)*par*pardz*pardz
      par = y1(i) - par
      fc = fc + y2(i)*par*par
    end do
    hs = -hs/fc ! XXX Missing respdf factor. Must be computed outside the fcn
  end subroutine logcondyzhs_gt

  subroutine logcondyzdz_sp (fc, gr, nu, y1, y2, z, n, tsq)
!! log-pdf of y|z and its derivative w.r.t. z.
    use modelfcns, only: invlink, invlinkdz, logpdfy, logpdfydlnk
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz
    fc = 0d0
    do i = 1, n
      par = invlink(z(i),nu)
      pardz = invlinkdz(z(i),nu)
      fc = fc + logpdfy(y1(i),y2(i),par)
      gr(i) = logpdfydlnk(y1(i),y2(i),par) * pardz
    end do
    fc = fc/tsq
    gr = gr/tsq
  end subroutine logcondyzdz_sp

  subroutine logcondyzhs_sp (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
    use modelfcns, only: invlink, invlinkdz, logpdfyhlnk, logpdfydlnk, invlinkhz
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision par, pardz
    do i = 1, n
      par = invlink(z(i),nu)
      pardz = invlinkdz(z(i),nu)
      hs(i) = logpdfyhlnk(y1(i),y2(i),par) * pardz*pardz
      hs(i) = hs(i) + logpdfydlnk(y1(i),y2(i),par) * invlinkhz(z(i),nu)
    end do
    hs = hs/tsq
  end subroutine logcondyzhs_sp
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
    ! The returned matrix x is upper triangular s.t. U'*U = X
    integer, intent(in) :: n
    double precision, intent(inout) :: x(n,n)
    double precision, intent(out) :: ldh
    integer i
!!     open(11,file='postermat.txt')
!!     write(11,*) x
!!     close(11)
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

subroutine aprxposterssq (fssq, meang, prechg, dz_dnu, dz_dphi, &
   ssq, nu, phi, omg, kappa, y1, y2, F, betm0, betQ0, &
   ssqdf, ssqsc, dm, tsq, tsqdf, n, p, ifam, icf)
  use covfun
  use betaprior
  use modelfcns
  use calcbd_fcns, only: qform, traceAB, cor_dcov
  implicit none
  integer, intent(in) :: n, p, ifam, icf
  double precision, intent(in) :: ssq, ssqdf, ssqsc, &
     y1(n), y2(n), phi, nu, omg, kappa, dm(n,n), F(n,p), &
     betm0(p), betQ0(p,p), tsq, tsqdf
  double precision, intent(out) :: fssq, meang(n), prechg(n,n), &
     dz_dnu(n), dz_dphi(n)
  double precision Ups(n,n), ldh_Ups, xi(n), modeldfh, ssqdfh, ssqdfsc
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), par(n), varh(n,n)
  integer j
  call create_model (ifam)
  call create_spcor (icf,n)
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)
  call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
  ssqdfh = .5d0*ssqdf
  ssqdfsc = ssqdf*ssqsc
  call poster (fssq, meang, prechg, ssq, ssqdfh, ssqdfsc, &
     y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  par = invlink(meang,nu)
  varh = prechg
  call dtrtri ('u','n',n,varh,n,j)
  dz_dnu = fdz_dnu (meang, par, ssq, varh, y1, y2, Ups, &
     nu, xi, lmxi, tsq, tsqdf, n)
  dz_dphi = fdz_dphi (meang, par, ssq, varh, y1, y2, dm, Ups, &
     nu, xi, lmxi, tsq, tsqdf, n, icf)

contains
!   function fdz_dnu (z, par, ssq, prechg, y1, y2, Ups, &
!      nu, xi, lmxi, tsq, tsqdf, n)
!     implicit none
!     integer, intent(in) :: n
!     logical, intent(in) :: lmxi
!     double precision, intent(in) :: y1(n), y2(n), Ups(n,n), nu, &
!        xi(n), tsq, tsqdf, ssq, z(n), prechg(n,n), par(n)
!     double precision, dimension(n) :: fdz_dnu
!     fdz_dnu = logpdfyhlnk(y1,y2,par) * invlinkdz(z,nu) * invlinkdn(z,nu) &
!        + logpdfydlnk(y1,y2,par) * invlinkdzdn(z,nu)
!     call dtrsv ('u','t','n',n,prechg,n,fdz_dnu,1)
!     call dtrsv ('u','n','n',n,prechg,n,fdz_dnu,1)
!   end function fdz_dnu

  function fdz_dnu (z, par, ssq, varh, y1, y2, Ups, &
     nu, xi, lmxi, tsq, tsqdf, n)
    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: lmxi
    double precision, intent(in) :: y1(n), y2(n), Ups(n,n), nu, &
       xi(n), tsq, tsqdf, ssq, z(n), varh(n,n), par(n)
    double precision, dimension(n) :: fdz_dnu
    fdz_dnu = logpdfyhlnk(y1,y2,par) * invlinkdz(z,nu) * invlinkdn(z,nu) &
       + logpdfydlnk(y1,y2,par) * invlinkdzdn(z,nu)
    fdz_dnu = fdz_dnu/tsq
    call dtrmv ('u','t','n',n,varh,n,fdz_dnu,1)
    call dtrmv ('u','n','n',n,varh,n,fdz_dnu,1)
  end function fdz_dnu

  function fdz_dphi (z, par, ssq, varh, y1, y2, dm, Ups, &
     nu, xi, lmxi, tsq, tsqdf, n, icf)
    implicit none
    integer, intent(in) :: n, icf
    logical, intent(in) :: lmxi
    double precision, intent(in) :: y1(n), y2(n), nu, dm(n,n), &
       xi(n), tsq, tsqdf, ssq, z(n), varh(n,n), par(n), Ups(n,n)
    double precision, dimension(n) :: fdz_dphi
    double precision zmxi(n), UpsDTUps(n,n), DT(n,n), &
       DTUps(n,n)
    DT = cor_dcov(n,dm,phi,omg,kappa,1)
    call fill_symmetric_matrix(DT,n)
    call dsymm ('r','u',n,n,1d0,Ups,n,DT,n,0d0,DTUps,n)
    call dsymm ('l','u',n,n,1d0,Ups,n,DTUps,n,0d0,UpsDTUps,n)
    UpsDTUps = UpsDTUps/ssq
    if (lmxi) then
      zmxi = z - xi
    else
      zmxi = z
    end if
    call dsymv ('u',n,1d0,UpsDTUps,n,zmxi,1,0d0,fdz_dphi,1) ! Ups*(z-xi)
    call dtrmv ('u','t','n',n,varh,n,fdz_dphi,1)
    call dtrmv ('u','n','n',n,varh,n,fdz_dphi,1)
  end function fdz_dphi
end subroutine aprxposterssq


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
  double precision, parameter :: eps = 1d-4
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
  ! print*, 0, real(tpnts), real(pdfpnts), real(pdfval)
  do i = 1, maxiter
    call polmax(tnew, ii, ok, tpnts, pdfpnts) ! Find new point
    call posterlog (pdfnew, mtmp, ptmp, tnew, ssqdfh, ssqdfsc, &
       y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
    lconv = .false.
    if (ok) lconv = abs(pdfnew - pdfval) .lt. eps*(1d0+sqrt(abs(pdfval)))
    if (ok .and. .not. lconv) lconv = minval(abs(tnew-tpnts)) .lt. &
       eps*1e2*(1d0+sqrt(abs(tnew)))
    ! print*, real(tnew), real(pdfnew), lconv
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
  ! print*, i, real(tpnts), real(pdfpnts), real(pdfval)
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
!!   open(11, file = "gridposter_data.txt")
!!   write(11,*) tg
!!   write(11,*) twght
!!   close(11)
!   twght(m) = 0d0
  !! twght = twght - log(sum(exp(twght))) ! XXX Not normalised.
end subroutine gridposter

subroutine aloglik (np, logssqg, out, meang, prechg, ssqdfh, ssqdfsc, &
   y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  implicit none
  integer, intent(in) :: np, n, ifam
  logical, intent(in) :: lmxi
  double precision, intent(in) :: y1(n), y2(n), Ups(n,n), nu, &
     xi(n), tsq, tsqdf, ssqdfh, ssqdfsc, ldh_Ups
  double precision, intent(in) :: logssqg(np+np+1)
  double precision, intent(out) :: out(np+np+1), meang(n,np+np+1), &
     prechg(n,n,np+np+1)
  integer i
  do i = 1, np+np+1
    call posterlog (out(i),meang(:,i),prechg(:,:,i),logssqg(i),ssqdfh,ssqdfsc,&
       y1,y2,Ups,ldh_Ups,nu,xi,lmxi,tsq,tsqdf,n,ifam)
  end do
end subroutine aloglik


subroutine aloglik_dnu (np, logssqg, dnu, meang, prechg, ssqdfh, ssqdfsc, &
   y1, y2, Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  use modelfcns
  implicit none
  integer, intent(in) :: np, n, ifam
  logical, intent(in) :: lmxi
  double precision, intent(in) :: y1(n), y2(n), Ups(n, n), nu, &
     xi(n), tsq, tsqdf, ssqdfh, ssqdfsc
  double precision, intent(in) :: logssqg(np+np+1), &
     meang(n,np+np+1), prechg(n,n,np+np+1)
  double precision, intent(out) :: dnu(np+np+1)
  double precision par(n), z(n), dz_dnu(n), dpar_dnu(n), &
     dpz_dz(n), ssq
  double precision dpym_dpar(n), dpym_hpar(n), dpym_3par(n)
  double precision invlink_dz(n), invlink_dn(n), invlink_hz(n), &
     invlink_dzdn(n), invlink_3z(n), invlink_hzdn(n)
  double precision dH(n)
  double precision varh(n,n)
  integer i, j
  do i = 1, np+np+1
    z = meang(:,i)
    ssq = exp(logssqg(i))
    varh = prechg(:,:,i)
    call dtrtri ('u','n',n,varh,n,j)
    if (j .ne. 0) call rexit ("aloglik_dnu - Non-invertible precision.")
    par = invlink(z,nu)
    invlink_dz = invlinkdz(z,nu)
    invlink_dn = invlinkdn(z,nu)
    invlink_hz = invlinkhz(z,nu)
    invlink_dzdn = invlinkdzdn(z,nu)
    invlink_3z = invlink3z(z,nu)
    invlink_hzdn = invlinkhzdn(z,nu)
    dpym_dpar = logpdfydlnk(y1,y2,par)
    dpym_hpar = logpdfyhlnk(y1,y2,par)
    dpym_3par = logpdfy3lnk(y1,y2,par)
    dpz_dz = logpdfzdz(z, Ups, xi, lmxi, ssq, n)
    dz_dnu = fdz_dnu(z, par, ssq, varh, &
       y1, y2, Ups, nu, xi, lmxi, tsq, tsqdf, n)
    dpar_dnu = invlink_dn + invlink_dz*dz_dnu
    dH = -dpym_3par*invlink_dz*invlink_dz*dpar_dnu &
       - 3d0*dpym_hpar*invlink_dz*invlink_hz*dz_dnu &
       - 2d0*dpym_hpar*invlink_dz*invlink_dzdn &
       - dpym_hpar*invlink_hz*invlink_dn &
       - dpym_dpar*invlink_3z*dz_dnu &
       - dpym_dpar*invlink_hzdn
    dH = dH/tsq
    dnu(i) = dot_product(dpym_dpar, dpar_dnu)/tsq &
       + dot_product(dpz_dz, dz_dnu) - .5d0*traceH(varh,dH,n)
  end do

contains
  function logpdfzdz (z, Ups, xi, lmxi, ssq, n) result (gr)
!! log-pdf of z and its derivative after integrating out beta. The 2*pi
!! constant is removed.
    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: lmxi
    double precision, intent(in) :: z(n), Ups(n,n), xi(n), ssq
    double precision :: gr(n)
    double precision zmxi(n)
    if (lmxi) then
      zmxi = z - xi
    else
      zmxi = z
    end if
    call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,gr,1) ! gr = Ups*(z-xi)
    gr = -gr/ssq ! gr = -Ups*(z-x)/ssq
  end function logpdfzdz

  function fdz_dnu (z, par, ssq, varh, y1, y2, Ups, &
     nu, xi, lmxi, tsq, tsqdf, n)
    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: lmxi
    double precision, intent(in) :: y1(n), y2(n), Ups(n,n), nu, &
       xi(n), tsq, tsqdf, ssq, z(n), varh(n,n), par(n)
    double precision, dimension(n) :: fdz_dnu
    fdz_dnu = logpdfyhlnk(y1,y2,par) * invlinkdz(z,nu) * invlinkdn(z,nu) &
       + logpdfydlnk(y1,y2,par) * invlinkdzdn(z,nu)
    fdz_dnu = fdz_dnu/tsq
    call dtrmv ('u','t','n',n,varh,n,fdz_dnu,1)
    call dtrmv ('u','n','n',n,varh,n,fdz_dnu,1)
  end function fdz_dnu

  function traceH (varh, dh, n)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: varh(n,n), dh(n)
    double precision traceH
    integer i
    double precision v(n)
    traceH = 0d0
    do i = 1, n
      v(1:i) = dh(1:i)*varh(1:i,i)
      traceH = traceH + dot_product(varh(1:i,i),v(1:i))
    end do
  end function traceH
end subroutine aloglik_dnu


subroutine aloglik_dcov (np, logssqg, dcov, ideriv, &
   meang, prechg, ssqdfh, ssqdfsc, &
   y1, y2, dm, phi, omg, kappa, Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  use modelfcns
  use covfun
  use calcbd_fcns, only: qform, traceAB, cor_dcov
  implicit none
  integer, intent(in) :: np, n, ifam, ideriv
  logical, intent(in) :: lmxi
  double precision, intent(in) :: y1(n), y2(n), Ups(n, n), nu, &
     xi(n), tsq, tsqdf, ssqdfh, ssqdfsc
  double precision, intent(in) :: dm(n,n), phi, omg, kappa
  double precision, intent(in) :: logssqg(np+np+1), &
     meang(n,np+np+1), prechg(n,n,np+np+1)
  double precision, intent(out) :: dcov(np+np+1)
  double precision par(n), z(n), dz_dcov(n), dpar_dcov(n), &
     dpz_dz(n), ssq
  double precision dpym_dpar(n), dpym_hpar(n), dpym_3par(n)
  double precision invlink_dz(n), invlink_hz(n), invlink_3z(n)
  double precision dH1(n), dH(n,n)
  double precision varh(n,n), DT(n,n), trUpsDTh, UpsDTUps(n,n), DTUps(n,n)
  double precision dpz_dcov
  integer i, j
  DT = cor_dcov(n,dm,phi,omg,kappa,ideriv)
  trUpsDTh = .5d0*traceAB(Ups,DT,n)
  call fill_symmetric_matrix(DT,n)
  call dsymm ('r','u',n,n,1d0,Ups,n,DT,n,0d0,DTUps,n)
  call dsymm ('l','u',n,n,1d0,Ups,n,DTUps,n,0d0,UpsDTUps,n)
  do i = 1, np+np+1
    z = meang(:,i)
    ssq = exp(logssqg(i))
    dH = UpsDTUps/ssq
    varh = prechg(:,:,i)
    call dtrtri ('u','n',n,varh,n,j)
    if (j .ne. 0) call rexit ("aloglik_dcov - Non-invertible precision.")
    par = invlink(z,nu)
    invlink_dz = invlinkdz(z,nu)
    invlink_hz = invlinkhz(z,nu)
    invlink_3z = invlink3z(z,nu)
    dpym_dpar = logpdfydlnk(y1,y2,par)
    dpym_hpar = logpdfyhlnk(y1,y2,par)
    dpym_3par = logpdfy3lnk(y1,y2,par)
    dpz_dz = logpdfzdz(z, Ups, xi, lmxi, ssq, n)
    dpz_dcov = logpdfzdcov(z, xi, lmxi, dH, trUpsDth, n)
    dz_dcov = fdz_dcov(z, par, ssq, varh, &
       y1, y2, dH, nu, xi, lmxi, tsq, tsqdf, n)
    dpar_dcov = invlink_dz*dz_dcov
    dH1 = -dpym_3par*invlink_dz*invlink_dz*dpar_dcov &
       - 3d0*dpym_hpar*invlink_dz*invlink_hz*dz_dcov &
       - dpym_dpar*invlink_3z*dz_dcov
    dH1 = dH1/tsq
    dH = -dH
    do j = 1, n
      dH(j,j) = dH(j,j) + dH1(j)
    end do
    dcov(i) = dpz_dcov &
       + dot_product(dpym_dpar, dpar_dcov)/tsq &
       + dot_product(dpz_dz, dz_dcov) &
       - .5d0*traceH(varh,dH,n)
  end do

contains
  function logpdfzdz (z, Ups, xi, lmxi, ssq, n) result (gr)
!! log-pdf of z and its derivative after integrating out beta. The 2*pi
!! constant is removed.
    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: lmxi
    double precision, intent(in) :: z(n), Ups(n,n), xi(n), ssq
    double precision :: gr(n)
    double precision zmxi(n)
    if (lmxi) then
      zmxi = z - xi
    else
      zmxi = z
    end if
    call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,gr,1) ! gr = Ups*(z-xi)
    gr = -gr/ssq ! gr = -Ups*(z-x)/ssq
  end function logpdfzdz

  function logpdfzdcov (z, xi, lmxi, UpsDTUps, trUpsDTh, n)
!! log-pdf of z and its derivative after integrating out beta. The 2*pi
!! constant is removed.
    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: lmxi
    double precision, intent(in) :: z(n), xi(n), UpsDTUps(n,n), trUpsDTh
    double precision :: logpdfzdcov
    double precision zmxi(n), zUDTUz
    if (lmxi) then
      zmxi = z - xi
    else
      zmxi = z
    end if
    zUDTUz = .5d0*qform(zmxi,UpsDTUps,n)
    logpdfzdcov = -trUpsDTh + zUDTUz
  end function logpdfzdcov

  function fdz_dcov (z, par, ssq, varh, y1, y2, UpsDTUps, &
     nu, xi, lmxi, tsq, tsqdf, n)
    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: lmxi
    double precision, intent(in) :: y1(n), y2(n), UpsDTUps(n,n), nu, &
       xi(n), tsq, tsqdf, ssq, z(n), varh(n,n), par(n)
    double precision, dimension(n) :: fdz_dcov
    double precision zmxi(n)
    if (lmxi) then
      zmxi = z - xi
    else
      zmxi = z
    end if
    call dsymv ('u',n,1d0,UpsDTUps,n,zmxi,1,0d0,fdz_dcov,1) ! Ups*(z-xi)
    call dtrmv ('u','t','n',n,varh,n,fdz_dcov,1)
    call dtrmv ('u','n','n',n,varh,n,fdz_dcov,1)
  end function fdz_dcov

  function traceH (varh, dH, n)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: varh(n,n), dH(n,n)
    double precision traceH
    integer i
    double precision v(n,n)
    v = dH
    call dtrmm ('r','u','n','n',n,n,1d0,varh,n,v,n)
    traceH = 0d0
    do i = 1, n
      traceH = traceH + dot_product(varh(1:i,i),v(1:i,i))
    end do
  end function traceH
end subroutine aloglik_dcov


subroutine llikpars2 (fval, gval, lderiv, &
   nu, phi, omg, kappa, y1, y2, F, betm0, betQ0, &
   ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, ssqin, ifam, icf)
  use covfun, only: create_spcor, calc_cov
  use betaprior
  implicit none
  integer, intent(in) :: np, n, p, ifam, icf
  logical, intent(in) :: lderiv(4) ! Which derivatives to calculate?
  double precision, intent(in) :: ssqin, ssqdf, ssqsc, &
     y1(n), y2(n), phi, nu, omg, kappa, dm(n,n), F(n,p), &
     betm0(p), betQ0(p,p), tsq, tsqdf
  double precision, intent(out) :: fval, gval(4)
  double precision Ups(n,n), ldh_Ups, xi(n), modeldfh, ssqdfh, ssqdfsc
  integer i
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p)
  double precision tg(np+np+1), logw(np+np+1), &
     meang(n,np+np+1), prechg(n,n,np+np+1)
  double precision, dimension(np+np+1) :: w_dnu, w_dphi, w_dnsq, w_dkap
  double precision ssqst
  call create_spcor(icf,n)
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)
  call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
  ssqdfh = .5d0*ssqdf
  ssqdfsc = ssqdf*ssqsc
  ssqst = ssqstart(y1,y2,nu,Ups,n,ifam)
  ! TODO ssqin or ssqst?
  ! ssqst = ssqin
  i = np+np+1
  call gridposter (np, tg, logw, meang, prechg, ssqdfh, ssqdfsc, &
     ssqst, y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  fval = trapezoid_lf(i, tg, logw)
  gval = 0d0
  if (lderiv(1)) then
    call aloglik_dnu (np, tg, w_dnu, meang, prechg, ssqdfh, ssqdfsc, &
       y1, y2, Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
    w_dnu = w_dnu*exp(logw - fval)
    gval(1) = trapezoid_f(i, tg, w_dnu)
  end if
  if (lderiv(2)) then
    call aloglik_dcov (np, tg, w_dphi, 1, &
       meang, prechg, ssqdfh, ssqdfsc, &
       y1, y2, dm, phi, omg, kappa, Ups, nu, &
       xi, lmxi, tsq, tsqdf, n, ifam)
    w_dphi = w_dphi*exp(logw - fval)
    gval(2) = trapezoid_f(i, tg, w_dphi)
  end if
  if (lderiv(3)) then
    call aloglik_dcov (np, tg, w_dnsq, 2, &
       meang, prechg, ssqdfh, ssqdfsc, &
       y1, y2, dm, phi, omg, kappa, Ups, nu, &
       xi, lmxi, tsq, tsqdf, n, ifam)
    w_dnsq = w_dnsq*exp(logw - fval)
    gval(3) = trapezoid_f(i, tg, w_dnsq)
  end if
  if (lderiv(4)) then
    call aloglik_dcov (np, tg, w_dkap, 3, &
       meang, prechg, ssqdfh, ssqdfsc, &
       y1, y2, dm, phi, omg, kappa, Ups, nu, &
       xi, lmxi, tsq, tsqdf, n, ifam)
    w_dkap = w_dkap*exp(logw - fval)
    gval(4) = trapezoid_f(i, tg, w_dkap)
  end if
  ! print*, real(fval), real(gval)

contains
  function trapezoid_lf (n, x, f)
    ! Compute the log integral of exp(f) along x
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(1:n), f(1:n)
    double precision trapezoid_lf
    double precision fpf(1:n-1), xmx(1:n-1), fmx, fmm(1:n), ln2
    parameter (ln2 = 0.6931471805599453d0)
    fmx = maxval(f)
    fmm = f - fmx
    fmm = exp(fmm)
    fpf = fmm(2:n) + fmm(1:n-1)
    xmx = x(2:n) - x(1:n-1)
    trapezoid_lf = dot_product(xmx,fpf)
    trapezoid_lf = log(trapezoid_lf) + fmx - ln2
  end function trapezoid_lf

  function trapezoid_f (n, x, f)
    ! Compute the integral of f along x
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(1:n), f(1:n)
    double precision trapezoid_f
    double precision fpf(1:n-1), xmx(1:n-1)
    fpf = f(2:n) + f(1:n-1)
    xmx = x(2:n) - x(1:n-1)
    trapezoid_f = dot_product(xmx,fpf)*.5d0
  end function trapezoid_f

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

subroutine llikparsval (fval, gval, ideriv, &
   nu, phi, omg, kappa, y1, y2, F, betm0, betQ0, &
   ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, ssqin, ifam, icf)
  use modelfcns, only: create_model
  implicit none
  integer, intent(in) :: np, n, p, ifam, icf, ideriv(4)
  double precision, intent(in) :: ssqin, ssqdf, ssqsc, &
     y1(n), y2(n), phi, nu, omg, kappa, dm(n,n), F(n,p), &
     betm0(p), betQ0(p,p), tsq, tsqdf
  double precision, intent(out) :: fval, gval(4)
  logical lderiv(4)
  call create_model (ifam)
  lderiv = ideriv .ne. 0
  call llikpars2 (fval, gval, lderiv, &
     nu, phi, omg, kappa, y1, y2, F, betm0, betQ0, &
   ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, ssqin, ifam, icf)
end subroutine llikparsval

subroutine aloglikval (fval, gval, &
   nu, phi, omg, kappa, y1, y2, F, betm0, betQ0, &
   ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, logssqg, ifam, icf)
  use modelfcns, only: create_model
  use covfun, only: create_spcor, calc_cov
  use betaprior
  implicit none
  integer, intent(in) :: np, n, p, ifam, icf
  double precision, intent(in) :: logssqg(np+np+1), ssqdf, ssqsc, &
     y1(n), y2(n), phi, nu, omg, kappa, dm(n,n), F(n,p), &
     betm0(p), betQ0(p,p), tsq, tsqdf
  double precision, intent(out) :: fval(np+np+1), gval(np+np+1)
  double precision Ups(n,n), ldh_Ups, xi(n), modeldfh, ssqdfh, ssqdfsc
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p)
  double precision meang(n,np+np+1), prechg(n,n,np+np+1)
  call create_model (ifam)
  call create_spcor(icf,n)
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)
  call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
  ssqdfh = .5d0*ssqdf
  ssqdfsc = ssqdf*ssqsc
  call aloglik (np, logssqg, fval, meang, prechg, ssqdfh, ssqdfsc, &
     y1, y2, Ups, ldh_Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
  call aloglik_dnu (np, logssqg, gval, meang, prechg, ssqdfh, ssqdfsc, &
     y1, y2, Ups, nu, xi, lmxi, tsq, tsqdf, n, ifam)
end subroutine aloglikval

subroutine llikparscalc (fval, nu, phi, omg, kappa, npars, &
   y1, y2, F, betm0, betQ0, &
   ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, ssqin, ifam, icf)
  ! Calls subroutine llikpars2 repeatedly for each parameter input.
  use modelfcns
  implicit none
  integer, intent(in) :: np, n, p, ifam, icf, npars
  double precision, intent(in) :: ssqin, ssqdf, ssqsc, &
     y1(n), y2(n), phi(npars), nu(npars), omg(npars), kappa(npars), &
     dm(n,n), F(n,p), betm0(p), betQ0(p,p), tsq, tsqdf
  double precision, intent(out) :: fval(npars)
  integer i
  logical lderiv(4)
  double precision gval(4)
  call create_model (ifam)
  lderiv = .false.
  do i = 1, npars
    call llikpars2 (fval(i), gval, lderiv, &
       nu(i), phi(i), omg(i), kappa(i), y1, y2, F, betm0, betQ0, &
       ssqdf, ssqsc, dm, tsq, tsqdf, n, p, np, ssqin, ifam, icf)
  end do
end subroutine llikparscalc


! Local Variables:
! compile-command: "gfortran -c -fpic -Wunused-parameter -Wall \
!   -pedantic -o skelpnts.o skelpnts.f90"
! End:

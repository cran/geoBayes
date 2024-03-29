!!! Compute the BF standard errors

!!!   Compute BF and SE
!!!
!!! @param logbf        The log(BF) (output)
!!! @param Sig  The covariance matrix of exp(logbf)
!!! @param SE   The SE value (output)
!!! @param VT1  The first term of variance
!!! @param VT2  The second term of variance
!!! @param iref Which model is the reference model
!!! @param phi  A vector of values to compute the SE
!!! @param omg  A vector of values to compute the SE
!!! @param nu   A vector of values to compute the SE
!!! @param kappa        A vector of values to compute the SE
!!! @param philist      A vector of values to compute the BF
!!! @param omglist      A vector of values to compute the BF
!!! @param nulist       A vector of values to compute the BF
!!! @param kappalist    A vector of values to compute the BF
!!! @param sample1      The sample to be used to compute the BF
!!! @param Nout1        The sample size from each model
!!! @param Ntot1        sum(Nout1)
!!! @param sample2      The sample to be used to compute the SE
!!! @param Nout2        The sample size from each model
!!! @param Ntot2        sum(Nout2)
!!! @param y    Vector of observations
!!! @param l    Vector of replications
!!! @param F    Design matrix
!!! @param dm   Distance matrix
!!! @param betm0        Prior mean for beta
!!! @param betQ0        Prior precision for beta
!!! @param ssqdf        Prior df for sigma^2
!!! @param ssqsc        Prior scale for sigma^2
!!! @param tsqdf        Prior df for tau^2 if transformed Gaussian model
!!! @param tsq  Prior scale for tau^2 for transf Ga model or dispersion ow
!!! @param icf	Correlation function
!!! @param n	Spatial dimension
!!! @param p	length of beta
!!! @param nnew	length of phi
!!! @param kg	number of models
!!! @param ifam	family
!!! @param imeth	Method to estimate the BF
!!! @param nb1	Number of batches to create when computing the 1st SE term
!!! @param nb2	Number of batches to create when computing the 2nd SE term
!!! @param ibvmeth	Which method to use for batch variance calculation
!!! @return SE, logbf, Sig, VT1, VT2

subroutine bfse_no (bf, logbfnew, Sig, SE, VT1, VT2, iref, &
   phi, omg, nu, kappa, &
   philist, omglist, nulist, kappalist, &
   sample1, Nout1, Ntot1, sample2, Nout2, Ntot2, &
   y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, nnew, kg, ifam, imeth, nb1, nb2, ibvmeth, itr)
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf, iref, nb1(kg), nnew, nb2(kg), ibvmeth, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     sample1(n,Ntot1), sample2(n,Ntot2), y(n), l(n), F(n,p), offset(n), &
     dm(n,n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, &
     kappalist(kg), phi(nnew), omg(nnew), nu(nnew), kappa(nnew)
  double precision, intent(out) :: bf(kg), SE(nnew), Sig(1:kg-1,1:kg-1), &
     VT1(nnew), VT2(nnew), logbfnew(nnew)
  double precision llik1(Ntot1,kg), llik2(Ntot2,kg), &
     llikn(Ntot2,nnew)
  double precision Bet(kg,kg), OOmg(kg,kg)

  ! Compute log-likelihood values
  call llikfcn_no (llik1, philist, omglist, nulist, kappalist, sample1, &
     Ntot1, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, kg, ifam, itr)

  call llikfcn_no (llik2, philist, omglist, nulist, kappalist, sample2, &
     Ntot2, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, kg, ifam, itr)

  call llikfcn_no (llikn, phi, omg, nu, kappa, sample2, &
     Ntot2, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, nnew, ifam, itr)

  call bfsecalc (bf, logbfnew, Sig, SE, VT1, VT2, iref, &
   llik1, llik2, llikn, Nout1, Ntot1, Nout2, Ntot2, &
   nnew, kg, imeth, nb1, nb2, ibvmeth, Bet, OOmg)
end subroutine bfse_no

subroutine bfse_mu (bf, logbfnew, Sig, SE, VT1, VT2, iref, &
   phi, omg, nu, kappa, &
   philist, omglist, nulist, kappalist, &
   sample1, Nout1, Ntot1, sample2, Nout2, Ntot2, &
   y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, nnew, kg, ifam, imeth, nb1, nb2, ibvmeth, itr)
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf, iref, nb1(kg), nnew, nb2(kg), ibvmeth, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     sample1(n,Ntot1), sample2(n,Ntot2), y(n), l(n), F(n,p), offset(n), &
     dm(n,n), betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, &
     kappalist(kg), phi(nnew), omg(nnew), nu(nnew), kappa(nnew)
  double precision, intent(out) :: bf(kg), SE(nnew), Sig(1:kg-1,1:kg-1), &
     VT1(nnew), VT2(nnew), logbfnew(nnew)
  double precision llik1(Ntot1,kg), llik2(Ntot2,kg), &
     llikn(Ntot2,nnew)
  double precision Bet(kg,kg), OOmg(kg,kg)

  ! Compute log-likelihood values
  call llikfcn_mu (llik1, philist, omglist, nulist, kappalist, sample1, &
     Ntot1, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, kg, ifam, itr)

  call llikfcn_mu (llik2, philist, omglist, nulist, kappalist, sample2, &
     Ntot2, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, kg, ifam, itr)

  call llikfcn_mu (llikn, phi, omg, nu, kappa, sample2, &
     Ntot2, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, nnew, ifam, itr)

  call bfsecalc (bf, logbfnew, Sig, SE, VT1, VT2, iref, &
     llik1, llik2, llikn, Nout1, Ntot1, Nout2, Ntot2, &
     nnew, kg, imeth, nb1, nb2, ibvmeth, Bet, OOmg)
end subroutine bfse_mu

subroutine bfse_wo (bf, logbfnew, Sig, SE, VT1, VT2, iref, &
   phi, omg, nu, kappa, &
   philist, omglist, nulist, kappalist, &
   sample1, Nout1, Ntot1, sample2, Nout2, Ntot2, &
   y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, nnew, kg, ifam, imeth, nb1, nb2, ibvmeth, itr)
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf, iref, nb1(kg), nnew, nb2(kg), ibvmeth, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     sample1(n,Ntot1), sample2(n,Ntot2), y(n), l(n), F(n,p), offset(n), &
     dm(n,n), betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, &
     kappalist(kg), phi(nnew), omg(nnew), nu(nnew), kappa(nnew)
  double precision, intent(out) :: bf(kg), SE(nnew), Sig(1:kg-1,1:kg-1), &
     VT1(nnew), VT2(nnew), logbfnew(nnew)
  double precision llik1(Ntot1,kg), llik2(Ntot2,kg), &
     llikn(Ntot2,nnew)
  double precision Bet(kg,kg), OOmg(kg,kg)

  ! Compute log-likelihood values
  call llikfcn_wo (llik1, philist, omglist, nulist, kappalist, sample1, &
     Ntot1, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, kg, ifam, itr)

  call llikfcn_wo (llik2, philist, omglist, nulist, kappalist, sample2, &
     Ntot2, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, kg, ifam, itr)

  call llikfcn_wo (llikn, phi, omg, nu, kappa, sample2, &
     Ntot2, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, nnew, ifam, itr)

  call bfsecalc (bf, logbfnew, Sig, SE, VT1, VT2, iref, &
   llik1, llik2, llikn, Nout1, Ntot1, Nout2, Ntot2, &
   nnew, kg, imeth, nb1, nb2, ibvmeth, Bet, OOmg)
end subroutine bfse_wo

subroutine bfse_tr (bf, logbfnew, Sig, SE, VT1, VT2, iref, &
   phi, omg, nu, kappa, &
   philist, omglist, nulist, kappalist, &
   sample1, Nout1, Ntot1, sample2, Nout2, Ntot2, &
   y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, nnew, kg, ifam, imeth, nb1, nb2, ibvmeth, itr)
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf, iref, nb1(kg), nnew, nb2(kg), ibvmeth, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     sample1(n,Ntot1), sample2(n,Ntot2), y(n), l(n), F(n,p), offset(n), &
     dm(n,n), betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, &
     kappalist(kg), phi(nnew), omg(nnew), nu(nnew), kappa(nnew)
  double precision, intent(out) :: bf(kg), SE(nnew), Sig(1:kg-1,1:kg-1), &
     VT1(nnew), VT2(nnew), logbfnew(nnew)
  double precision llik1(Ntot1,kg), llik2(Ntot2,kg), &
     llikn(Ntot2,nnew)
  double precision Bet(kg,kg), OOmg(kg,kg)

  ! Compute log-likelihood values
  call llikfcn_tr (llik1, philist, omglist, nulist, kappalist, sample1, &
     Ntot1, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, kg, ifam, itr)

  call llikfcn_tr (llik2, philist, omglist, nulist, kappalist, sample2, &
     Ntot2, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, kg, ifam, itr)

  call llikfcn_tr (llikn, phi, omg, nu, kappa, sample2, &
     Ntot2, y, l, F, offset, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, icf, &
     n, p, nnew, ifam, itr)

  call bfsecalc (bf, logbfnew, Sig, SE, VT1, VT2, iref, &
   llik1, llik2, llikn, Nout1, Ntot1, Nout2, Ntot2, &
   nnew, kg, imeth, nb1, nb2, ibvmeth, Bet, OOmg)
end subroutine bfse_tr



subroutine bfsecalc (bf, logbfnew, Sig, SE, VT1, VT2, iref, &
   llik1, llik2, llikn, Nout1, Ntot1, Nout2, Ntot2, &
   nnew, kg, imeth, nb1, nb2, ibvmeth, Bet, OOmg)
  use bmargin
  implicit none
  integer, intent(in) :: kg, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, iref, nb1(kg), nnew, nb2(kg), ibvmeth
  double precision, intent(in) :: &
     llik1(Ntot1,kg), llik2(Ntot2,kg), llikn(Ntot2,nnew)
  double precision, intent(out) :: bf(kg), SE(nnew), Sig(1:kg-1,1:kg-1), &
     VT1(nnew), VT2(nnew), OOmg(kg,kg), logbfnew(nnew)
  double precision eta(kg), logbf(kg)
  double precision logYY(Ntot1,kg), YY(Ntot1,kg), logVV(Ntot2,nnew), &
     VV(Ntot2,nnew), VVcol(Ntot2)
  double precision Bet(kg,kg), BMP(kg,kg), BD(kg,1:kg-1), &
     OOmgBD(kg,1:kg-1), cvec(1:kg-1,nnew), logcall(kg,nnew), &
     SigC(1:kg-1)
  double precision NNratio(kg), logN(kg)
  integer i, j, kgm1, ia, ie
  double precision tmp

  if (kg .lt. 2) then
    logbf = 0d0
    bf = 1d0
    go to 9
  end if

  ! Estimate BF
  logN = log(dble(Nout1))
  eta = logN
  select case (imeth)
  case (1)
    call revlogistic (eta,llik1,kg,Ntot1,Nout1)
  case (2)
    call revlogistic (eta,llik1,kg,Ntot1,Nout1)
    call mengwong (eta,llik1,kg,Ntot1,Nout1)
  end select
  ! eta_i = log(N_i) - log(bf_i)
  logbf = logN - eta ! = log(bf_i)
  logbf = logbf - logbf(iref)
  bf = exp(logbf)
  eta = eta - sum(eta)/dble(kg) ! zeta

  logYY = logp(llik1,eta,Ntot1,kg)
  YY = exp(logYY)

  OOmg = 0d0
  ie = 0
  do j = 1, kg
    ia = ie + 1
    ie = ie + Nout1(j)
    Bet = bmmcvrmat(ibvmeth,YY(ia:ie,:),Nout1(j),kg,nb1(j))
    OOmg = OOmg + Bet*dble(Nout1(j))/dble(Ntot1)
  end do

  do j = 1, kg
    do i = 1, j-1
      Bet(i,j) = -dot_product(YY(:,i),YY(:,j))/dble(Ntot1)
    end do
    Bet(j,j) = dot_product(YY(:,j),1d0-YY(:,j))/dble(Ntot1)
  end do

  BMP = Bet + 1d0/dble(kg)
  call dpotrf("u",kg,BMP,kg,i)
  call dpotri("u",kg,BMP,kg,i)
  BMP = BMP - 1d0/dble(kg)

  j = 0
  eta = symcol(BMP,kg,iref)
  do i = 1, kg
    if (i .eq. iref) cycle
    j = j + 1
    BD(:,j) = (eta - symcol(BMP,kg,i))*bf(i) ! B * D
  end do

9 kgm1 = kg - 1
  logN = log(dble(Nout2))
  eta = logN - logbf
  logbfnew = logbfnew_calc(llik2,llikn,eta,Ntot2,kg,nnew)
  OOmgBD = 0d0
  Sig = 0d0
  SigC = 0d0
  VT1 = 0d0

  if (kgm1 .gt. 0) then
    call dsymm("l","u",kg,kgm1,1d0,OOmg,kg,BD,kg,0d0,OOmgBD,kg)
    call dgemm("t","n",kgm1,kgm1,kg,1d0,BD,kg,OOmgBD,kg,0d0,Sig,kgm1)

    logcall = logc(llik2,llikn,eta,Ntot2,kg,nnew) + spread(eta+eta-logN,2,nnew)
    cvec(:iref-1,:) = exp(logcall(:iref-1,:))
    cvec(iref:,:) = exp(logcall(iref+1:,:))

    do i = 1, nnew
      call dsymv ("u",kgm1,1d0,Sig,kgm1,cvec(:,i),1,0d0,SigC,1)
      VT1(i) = dot_product(cvec(:,i),SigC)
    end do
    VT1 = dble(Ntot2)*VT1/dble(Ntot1)
  end if

  logVV = logb(llik2,llikn,eta,Ntot2,kg,nnew) + log(dble(Ntot2))
  VV = exp(logVV)

  NNratio = dble(Nout2)/dble(Ntot2)
  VT2 = 0d0
  do i = 1, nnew
    VVcol = VV(:,i)
    ie = 0
    do j = 1, kg
      ia = ie + 1
      ie = ie + Nout2(j)
      tmp = bmmcvr(ibvmeth,VVcol(ia:ie),Nout2(j),nb2(j))
      VT2(i) = VT2(i) + tmp*NNratio(j)
    end do
  end do

  SE = sqrt(VT1 + VT2)
contains
  function symcol (A,n,i)
    ! Gives the ith column of a symmetric matrix when only its upper triangle
    ! is filled
    implicit none
    integer, intent(in) :: n, i
    double precision, intent(in) :: A(n,n)
    double precision symcol(n)
    symcol(:i) = A(:i,i)
    symcol(i+1:) = A(i,i+1:)
  end function symcol

  function logp (llik,eta,Ntot,kg)
    use flogsumexp, only: logrsumexp
    implicit none
    integer, intent(in) :: Ntot, kg
    double precision, intent(in) :: llik(Ntot,kg), eta(kg)
    double precision logp(Ntot,kg)
    double precision lliketa(Ntot,kg), lgdenom(Ntot)
    integer i, j
    do j = 1, kg
      do i = 1, Ntot
        lliketa(i,j) = llik(i,j) + eta(j)
      end do
    end do
    lgdenom = logrsumexp(lliketa,Ntot,kg)
    logp = spread(lgdenom,2,kg)
    logp = lliketa - logp
  end function logp

  function logb (llik1,llik2,eta,Ntot,kg,nnew)
    use flogsumexp, only: logrsumexp
    implicit none
    integer, intent(in) :: Ntot, kg, nnew
    double precision, intent(in) :: llik1(Ntot,kg), llik2(Ntot,nnew), eta(kg)
    double precision logb(Ntot,nnew)
    double precision lliketa(Ntot,kg), lgdenom(Ntot)
    integer i, j
    do j = 1, kg
      do i = 1, Ntot
        lliketa(i,j) = llik1(i,j) + eta(j)
      end do
    end do
    lgdenom = logrsumexp(lliketa,Ntot,kg)
    logb = spread(lgdenom,2,nnew)
    logb = llik2 - logb
  end function logb

  function logc (llik1,llik2,eta,Ntot,kg,nnew)
    ! Computes log c(h)_j without the N_j/d_j^2 term
    ! eta(j) = log(Nout2(j)) - log(BF(j))
    use flogsumexp
    implicit none
    integer, intent(in) :: Ntot, kg, nnew
    double precision, intent(in) :: llik1(Ntot,kg), llik2(Ntot,nnew), eta(kg)
    double precision logc(kg,nnew)
    double precision lliketa(Ntot,kg), lgdenom(Ntot), llikrat(Ntot,kg)
    integer i
    lliketa = spread(eta,1,Ntot)
    lliketa = llik1 + lliketa
    lgdenom = logrsumexp(lliketa,Ntot,kg)
    lgdenom = lgdenom + lgdenom
    llikrat = spread(lgdenom,2,kg)
    llikrat = llik1 - llikrat
    do i = 1, nnew
      lliketa = spread(llik2(:,i),2,kg)
      lliketa = lliketa + llikrat
      logc(:,i) = logcsumexp(lliketa,Ntot,kg)
    end do
  end function logc

  function bmmcvrmat_st (x,n,k,b)
    ! Split matrix x(n,k) to b batches of roughly equal size and compute the
    ! MC variance. The output is a k*k matrix.
    implicit none
    integer, intent(in) :: n, k, b
    double precision, intent(in) :: x(n,k)
    double precision bmmcvrmat_st(k,k)
    integer bsz(b), m, r, i
    double precision bmean(b,k), xmean(k), dbsz(b)
    m = n/b ! Minimum size for each batch
    r = mod(n,b)
    bsz(:r) = m + 1
    bsz(r+1:) = m ! Size of each batch
    dbsz = dble(bsz)
    do i = 1, k
      bmean(:,i) = batchmeans(x(:,i),n,b)
    end do
    xmean = matmul(dbsz,bmean)/dble(n)
    do i = 1, k
      bmean(:,i) = bmean(:,i) - xmean(i)
    end do
    bmmcvrmat_st = 0d0
    do i = 1, b
      xmean = bmean(i,:)
      call dsyr ("u",k,dbsz(i),xmean,1,bmmcvrmat_st,k)
    end do
    bmmcvrmat_st = bmmcvrmat_st/dble(b-1)
  end function bmmcvrmat_st

  function bmmcvr_st(x,n,b)
    ! Split vector x(n) to b batches of roughly equal size and compute the
    ! MC variance. The output is a scalar.
    implicit none
    integer, intent(in) :: n, b
    double precision, intent(in) :: x(n)
    double precision bmmcvr_st
    integer bsz(b), m, r
    double precision bmean(b), xmean, dbsz(b)
    m = n/b
    r = mod(n,b)
    bsz(:r) = m + 1
    bsz(r+1:) = m
    dbsz = dble(bsz)
    bmean = batchmeans(x,n,b)
    xmean = dot_product(bmean,dbsz)/dble(n)
    bmean = bmean - xmean
    bmean = bmean*bmean
    bmmcvr_st = dot_product(bmean,dbsz)/dble(b-1)
    ! bmmcvr_st = bmmcvr_st/dble(n)
  end function bmmcvr_st

  function bmmcvrmat (i,x,n,k,b)
    ! Split matrix x(n,k) to b batches of roughly equal size and compute the
    ! MC variance. The output is a k*k matrix.
    implicit none
    integer, intent(in) :: i, n, k, b
    double precision, intent(in) :: x(n,k)
    double precision bmmcvrmat(k,k)
    select case (i)
    case (1) ! Standard
      bmmcvrmat = bmmcvrmat_st(x,n,k,b)
    case (2) ! Spectral Tukey-Hanning
      bmmcvrmat = bmmcvrmat_th(x,n,k,b)
    case (3) ! Spectral modified Bartlett
      bmmcvrmat = bmmcvrmat_mb(x,n,k,b)
    end select
  end function bmmcvrmat

  function batchmeans (x,n,b)
    ! Split vector x(n) to b batches of roughly equal size and compute the
    ! mean of each batch. The output is a vector of size b.
    implicit none
    integer, intent(in) :: n, b
    double precision, intent(in) :: x(n)
    double precision batchmeans(b)
    integer m, i, ia, ie, r, d

    d = n/b ! Minimum batch size
    if (d .eq. 0) then
      batchmeans(1:n) = x
      batchmeans(n+1:b) = (d*b)/d
      return
    end if
    r = n - d*b
    m = d + 1
    ie = 0
    do i = 1, r
      ia = ie + 1
      ie = ie + m
      batchmeans(i) = mean(x(ia:ie),m)
    end do
    m = d
    do i = r+1, b
      ia = ie + 1
      ie = ie + m
      batchmeans(i) = mean(x(ia:ie),m)
    end do
  end function batchmeans

  function bmmcvrmat_th (x,n,k,b)
    implicit none
    integer, intent(in) :: n, k, b
    double precision, intent(in) :: x(n,k)
    double precision bmmcvrmat_th(k,k)
    integer i, j
    double precision xmm(n,k), x1(k), x2(k), wj
    double precision, parameter :: pi = 3.1415926535897932384626433832795d0
    x1 = sum(x,1)/dble(n)
    do j = 1, k
      do i = 1, n
        xmm(i,j) = x(i,j) - x1(j) ! Y_i - W
      end do
    end do
    bmmcvrmat_th = 0d0
    do i = 1, n ! j = 0
      x1 = xmm(i,:)
      call dsyr ('u',k,1d0,x1,1,bmmcvrmat_th,k)
    end do
    do j = 1, b - 1
      wj = .5d0*(1d0 + cos(pi*dble(j)/dble(b)))
      do i = 1, n - j
        x1 = xmm(i,:)
        x2 = xmm(i+j,:)
        call dsyr2 ('u',k,wj,x1,1,x2,1,bmmcvrmat_th,k)
      end do
    end do
    bmmcvrmat_th = bmmcvrmat_th/dble(n)
  end function bmmcvrmat_th

  function bmmcvr_th(x,n,b)
    ! Compute the univariate spectral variance estimator using the
    ! Tukey-Hanning lag window.
    implicit none
    integer, intent(in) :: n, b
    double precision, intent(in) :: x(n)
    double precision bmmcvr_th
    integer i, j
    double precision xmm(n), wj
    double precision, parameter :: pi = 3.1415926535897932384626433832795d0
    wj = sum(x)/dble(n) ! mean(x,n)
    xmm = x - wj
    bmmcvr_th = dot_product(xmm,xmm) ! j = 0
    do j = 1, b - 1
      wj = 1d0 + cos(pi*dble(j)/dble(b)) ! Without the 0.5 factor because
      ! it cancels with the factor 2 below.
      do i = 1, n - j
        bmmcvr_th = bmmcvr_th + wj*xmm(i)*xmm(i+j)
      end do
    end do
    bmmcvr_th = bmmcvr_th/dble(n)
  end function bmmcvr_th

  function bmmcvr (i,x,n,b)
    ! Compute the univariate spectral variance estimator using the
    ! Tukey-Hanning lag window.
    implicit none
    integer, intent(in) :: i, n, b
    double precision, intent(in) :: x(n)
    double precision bmmcvr
    select case (i)
    case (1) ! Standard
      bmmcvr = bmmcvr_st(x,n,b)
    case (2) ! Spectral Tukey-Hanning
      bmmcvr = bmmcvr_th(x,n,b)
    case (3) ! Spectral modified Bartlett
      bmmcvr = bmmcvr_mb(x,n,b)
    case default
      bmmcvr = -huge(bmmcvr)
    end select
  end function bmmcvr

  function bmmcvrmat_mb (x,n,k,b)
    implicit none
    integer, intent(in) :: n, k, b
    double precision, intent(in) :: x(n,k)
    double precision bmmcvrmat_mb(k,k)
    integer i, j
    double precision xmm(n,k), x1(k), x2(k), wj
    x1 = sum(x,1)/dble(n)
    do j = 1, k
      do i = 1, n
        xmm(i,j) = x(i,j) - x1(j) ! Y_i - W
      end do
    end do
    bmmcvrmat_mb = 0d0
    do i = 1, n ! j = 0
      x1 = xmm(i,:)
      call dsyr ('u',k,1d0,x1,1,bmmcvrmat_mb,k)
    end do
    do j = 1, b - 1
      wj = 1d0 - dble(j)/dble(b)
      do i = 1, n - j
        x1 = xmm(i,:)
        x2 = xmm(i+j,:)
        call dsyr2 ('u',k,wj,x1,1,x2,1,bmmcvrmat_mb,k)
      end do
    end do
    bmmcvrmat_mb = bmmcvrmat_mb/dble(n)
  end function bmmcvrmat_mb

  function bmmcvr_mb(x,n,b)
    ! Compute the univariate spectral variance estimator using the
    ! modified Bartlett lag window.
    implicit none
    integer, intent(in) :: n, b
    double precision, intent(in) :: x(n)
    double precision bmmcvr_mb
    integer i, j
    double precision xmm(n), wj
    wj = mean(x,n)
    xmm = x - wj
    bmmcvr_mb = dot_product(xmm,xmm) ! j = 0
    do j = 1, b - 1
      wj = 2d0*(1d0 - dble(j)/dble(b)) ! With an extra factor 2 beacause
      ! of the sum below
      do i = 1, n - j
        bmmcvr_mb = bmmcvr_mb + wj*xmm(i)*xmm(i+j)
      end do
    end do
    bmmcvr_mb = bmmcvr_mb/dble(n)
  end function bmmcvr_mb

  double precision function mean(x, m)
    use interfaces, only: isfinite
    implicit none
    integer, intent(in) :: m
    double precision, intent(in) :: x(m)
    integer i
    double precision tot, tot0, mm
    mean = 0d0
    tot = 0d0
    mm = dble(m)
    do i = 1, m
      tot0 = tot
      tot = tot + x(i)
      if (isfinite(tot) .eq. 0) then
        mean = mean + tot0/mm
        tot = x(i)
      end if
    end do
    mean = mean + tot/mm
  end function mean

  function logbfnew_calc (llik2,llikn,eta,Ntot,kg,nnew)
    use flogsumexp, only: logrsumexp, logcsumexp
    implicit none
    integer, intent(in) :: Ntot, kg, nnew
    double precision, intent(in) :: llik2(Ntot,kg), llikn(Ntot,nnew), &
       eta(kg)
    double precision logbfnew_calc(nnew)
    double precision logb(Ntot,nnew)
    double precision lliketa(Ntot,kg), lgdenom(Ntot)
    integer i, j
    if (nnew .gt. 0) then
      do j = 1, kg
        do i = 1, Ntot
          lliketa(i,j) = llik2(i,j) + eta(j)
        end do
      end do
      lgdenom = logrsumexp(lliketa,Ntot,kg)
      logb = spread(lgdenom,2,nnew)
      logb = llikn - logb
      logbfnew_calc = logcsumexp(logb,Ntot,nnew)
    end if
  end function logbfnew_calc
end subroutine bfsecalc



subroutine logbfcalc (logbf, logbfnew, iref, &
   llik1, llik2, llikn, Nout1, Ntot1, Nout2, Ntot2, &
   nnew, kg, imeth)
  use bmargin
  implicit none
  integer, intent(in) :: kg, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, nnew, iref
  double precision, intent(in) :: &
     llik1(Ntot1,kg), llik2(Ntot2,kg), llikn(Ntot2,nnew)
  double precision, intent(out) :: logbf(kg), logbfnew(nnew)
  double precision eta(kg), logN(kg)

  if (kg .lt. 2) then
    logbf = 0d0
    go to 9
  end if

  ! Estimate BF
  logN = log(dble(Nout1))
  eta = logN
  select case (imeth)
  case (1)
    call revlogistic (eta,llik1,kg,Ntot1,Nout1)
  case (2)
    call revlogistic (eta,llik1,kg,Ntot1,Nout1)
    call mengwong (eta,llik1,kg,Ntot1,Nout1)
  end select
  ! eta_i = log(N_i) - log(bf_i)
  logbf = logN - eta ! = log(bf_i)
  logbf = logbf - logbf(iref)

9 if (nnew > 0) then
    logN = log(dble(Nout2))
    eta = logN - logbf
    logbfnew = logbfnew_calc(llik2,llikn,eta,Ntot2,kg,nnew)
  end if

contains
  function logbfnew_calc (llik2,llikn,eta,Ntot,kg,nnew)
    use flogsumexp, only: logrsumexp, logcsumexp
    implicit none
    integer, intent(in) :: Ntot, kg, nnew
    double precision, intent(in) :: llik2(Ntot,kg), llikn(Ntot,nnew), &
       eta(kg)
    double precision logbfnew_calc(nnew)
    double precision logb(Ntot,nnew)
    double precision lliketa(Ntot,kg), lgdenom(Ntot)
    integer i, j
    do j = 1, kg
      do i = 1, Ntot
        lliketa(i,j) = llik2(i,j) + eta(j)
      end do
    end do
    lgdenom = logrsumexp(lliketa,Ntot,kg)
    logb = spread(lgdenom,2,nnew)
    logb = llikn - logb
    logbfnew_calc = logcsumexp(logb,Ntot,nnew)
  end function logbfnew_calc
end subroutine logbfcalc


subroutine bfsecalcef (Ef, VT2, logbfnew, evec, iref, &
   logbf,llik2,llikn,fval,kg,nnew,nf,Nout2,Ntot2,nb,ibvmeth)
  use flogsumexp, only: logrsumexp, logcsumexp, logsumexpv
  implicit none
  integer, intent(in) :: kg, nnew, Nout2(kg), Ntot2, nf, ibvmeth, nb(kg), &
     iref 
  double precision, intent(in) :: logbf(kg), llik2(Ntot2,kg), &
     llikn(Ntot2,nnew), fval(Ntot2,nf)
  double precision, intent(out) :: Ef(nnew,nf), VT2(nnew,nf), &
     logbfnew(nnew), evec(kg-1,nnew,nf)
  integer i, ia, ie, j
  double precision zeta(kg), logVV(Ntot2,nnew), lliketa(Ntot2,kg), &
     logdenom(Ntot2), VV(Ntot2,nnew), WW(Ntot2,nnew,nf), NNratio(kg), &
     logweight(Ntot2,nnew), weight(Ntot2,nnew), WWwght(Ntot2,nnew,nf), &
     logbfnorm(kg), GamWW(nnew,nf), GamVW(nnew,nf), GamVV(nnew), &
     invbfnew(nnew), logV1(Ntot2,kg-1), V1(Ntot2,kg-1), cvec(kg-1,nnew), &
     bf(kg), bbn(kg-1), invbfnewsq(nnew)
  NNratio = dble(Nout2)/dble(Ntot2)
  logbfnorm = logbf - logbf(iref)
  bf = exp(logbfnorm)
  do i = 1, iref - 1
    bbn(i) = NNratio(i)/(bf(i)*bf(i)*Ntot2)
  end do
  do i = iref+1, kg
    bbn(i-1) = NNratio(i)/(bf(i)*bf(i)*Ntot2)
  end do
  zeta = -logbfnorm + log(NNratio)
  do i = 1, kg
    lliketa(:,i) = llik2(:,i) + zeta(i)
  end do
  logdenom = logrsumexp(lliketa,Ntot2,kg)
  do i = 1, iref-1
    logV1(:,i) = llik2(:,i) - logdenom
  end do
  do i = iref+1, kg
    logV1(:,i-1) = llik2(:,i) - logdenom
  end do
  V1 = exp(logV1)
  do i = 1, nnew
    logVV(:,i) = llikn(:,i) - logdenom
  end do
  VV = exp(logVV)
  logbfnew = logcsumexp(logVV,Ntot2,nnew)
  j = kg-1
  call dgemm ("t","n",j,nnew,Ntot2,1d0,V1,Ntot2,VV,Ntot2,0d0,cvec,j)
  do i = 1, nnew
    cvec(:,i) = cvec(:,i)*bbn
  end do
  do i = 1, nf
    do j = 1, nnew
      WW(:,j,i) = fval(:,i)*VV(:,j)
    end do
  end do
  do i = 1, nnew
    logweight(:,i) = logVV(:,i) - logbfnew(i)
  end do
  weight = exp(logweight) ! Importance weights
  do i = 1, nf
    do j = 1, nnew
      WWwght(:,j,i) = fval(:,i)*weight(:,j)
    end do
  end do
  Ef = sum(WWwght,1) ! Estimation of expectation
  logbfnew = logbfnew - log(dble(Ntot2))
  invbfnew = exp(-logbfnew)
  j = kg-1
  i = nnew*nf
  call dgemm ("t","n",j,i,Ntot2,1d0,V1,Ntot2,WW,Ntot2,0d0,evec,j)
  do j = 1, nf
    do i = 1, nnew
      evec(:,i,j) = evec(:,i,j)*bbn - cvec(:,i)*Ef(i,j)
      evec(:,i,j) = evec(:,i,j)*invbfnew(i)
    end do
  end do
  GamWW = 0d0
  GamVW = 0d0
  GamVV = 0d0
  ie = 0
  do i = 1, kg
    ia = ie + 1
    ie = ie + Nout2(i)
    select case (ibvmeth)
    case (1)
      call covmat_accum_bm (GamWW, GamVW, GamVV,WW(ia:ie,:,:),VV(ia:ie,:),&
         NNratio(i), &
         Nout2(i),nnew,nf,nb(i))
    case (2)
      call covmat_accum_th (GamWW, GamVW, GamVV,WW(ia:ie,:,:),VV(ia:ie,:),&
         NNratio(i), &
         Nout2(i),nnew,nf,nb(i))
    case (3)
      call covmat_accum_mb (GamWW, GamVW, GamVV,WW(ia:ie,:,:),VV(ia:ie,:),&
         NNratio(i), &
         Nout2(i),nnew,nf,nb(i))
    end select
  end do
  invbfnewsq = invbfnew*invbfnew
  do i = 1, nf
    VT2(:,i) = GamWW(:,i) - (2*Ef(:,i))*GamVW(:,i) + (Ef(:,i)*Ef(:,i))*GamVV
    VT2(:,i) = VT2(:,i)*invbfnewsq
  end do
  
contains

  subroutine covmat_accum_bm (GamWW, GamVW, GamVV, WW, VV, NNratio, &
     Nout,nnew,nf,nb)
    implicit none
    integer, intent(in) :: Nout, nnew, nf, nb
    double precision, intent(in) :: WW(Nout,nnew,nf), VV(Nout,nnew), &
       NNratio
    double precision, intent(inout) :: GamWW(nnew,nf), GamVW(nnew,nf), &
       GamVV(nnew)
    integer i, j, k, l, batch_size(nb), ia, ie
    double precision Wb(nb,nnew,nf), Vb(nb,nnew), Wm(nnew,nf), Vm(nnew), &
       Wbm(nb,nnew,nf), Vbm(nb,nnew), s, wgt
    batch_size = Nout/nb
    i = mod(Nout,nb)
    batch_size(:i) = batch_size(:i) + 1
    do l = 1, nf
      do k = 1, nnew
        Wm(k,l) = 0d0
        ie = 0
        do i = 1, nb
          ia = ie + 1
          ie = ie + batch_size(i)
          s = 0d0
          do j = ia, ie
            s = s + WW(j,k,l)
          end do
          Wm(k,l) = Wm(k,l) + s
          Wb(i,k,l) = s/dble(batch_size(i))
        end do
        Wm(k,l) = Wm(k,l)/dble(Nout)
        Wbm(:,k,l) = Wb(:,k,l) - Wm(k,l)
      end do
    end do
    do k = 1, nnew
      Vm(k) = 0d0
      ie = 0
      do i = 1, nb
        ia = ie + 1
        ie = ie + batch_size(i)
        s = 0d0
        do j = ia, ie
          s = s + VV(j,k)
        end do
        Vm(k) = Vm(k) + s
        Vb(i,k) = s/dble(batch_size(i))
      end do
      Vm(k) = Vm(k)/dble(Nout)
      Vbm(:,k) = Vb(:,k) - Vm(k)
    end do
    do l = 1, nf
      do k = 1, nnew
        do i = 1, nb
          wgt = NNratio*dble(batch_size(i))/dble(nb-1)
          GamWW(k,l) = GamWW(k,l) + wgt*Wbm(i,k,l)*Wbm(i,k,l)
        end do
      end do
    end do
    do k = 1, nnew
      do i = 1, nb
        wgt = NNratio*dble(batch_size(i))/dble(nb-1)
        GamVV(k) = GamVV(k) + wgt*Vbm(i,k)*Vbm(i,k)
      end do
    end do
    do l = 1, nf
      do k = 1, nnew
        do i = 1, nb
          wgt = NNratio*dble(batch_size(i))/dble(nb-1)
          GamVW(k,l) = GamVW(k,l) + wgt*Vbm(i,k)*Wbm(i,k,l)
        end do
      end do
    end do
  end subroutine covmat_accum_bm
  
  subroutine covmat_accum_th (GamWW, GamVW, GamVV, WW, VV, NNratio, &
     Nout,nnew,nf,nb)
    implicit none
    integer, intent(in) :: Nout, nnew, nf, nb
    double precision, intent(in) :: WW(Nout,nnew,nf), VV(Nout,nnew), &
       NNratio
    double precision, intent(inout) :: GamWW(nnew,nf), GamVW(nnew,nf), &
       GamVV(nnew)
    integer i, j, k, l
    double precision wgt, Wm(nnew,nf), Vm(nnew), Wmm(Nout,nnew,nf), &
       Vmm(Nout,nnew)
    double precision, parameter :: pi = 3.1415926535897932384626433832795d0 
    Vm = sum(VV,1)/dble(Nout)
    Wm = sum(WW,1)/dble(Nout)
    do k = 1, nf
      do j = 1, nnew
        do i = 1, Nout
          Wmm(i,j,k) = WW(i,j,k) - Wm(j,k)
        end do
      end do
    end do
    do j = 1, nnew
      do i = 1, Nout
        Vmm(i,j) = VV(i,j) - Vm(j)
      end do
    end do
    ! j = 0
    wgt = NNratio/dble(Nout)
    do l = 1, nf
      do k = 1, nnew
        do i = 1, Nout
          GamWW(k,l) = GamWW(k,l) + wgt*Wmm(i,k,l)*Wmm(i,k,l)
        end do
      end do
    end do
    do k = 1, nnew
      do i = 1, Nout
        GamVV(k) = GamVV(k) + wgt*Vmm(i,k)*Vmm(i,k)
      end do
    end do
    do l = 1, nf
      do k = 1, nnew
        do i = 1, Nout
          GamVW(k,l) = GamVW(k,l) + wgt*Vmm(i,k)*Wmm(i,k,l)
        end do
      end do
    end do
    do j = 1, nb - 1
      wgt = NNratio*(1d0 + cos((pi*j)/nb))/Nout ! w/o the .5 factor because it
      ! cancels
      do l = 1, nf
        do k = 1, nnew
          do i = 1, Nout-j
            GamWW(k,l) = GamWW(k,l) + wgt*Wmm(i,k,l)*Wmm(i+j,k,l)
          end do
        end do
      end do
      do k = 1, nnew
        do i = 1, Nout-j
          GamVV(k) = GamVV(k) + wgt*Vmm(i,k)*Vmm(i+j,k)
        end do
      end do
      wgt = .5d0*wgt
      do l = 1, nf
        do k = 1, nnew
          do i = 1, Nout-j
            GamVW(k,l) = GamVW(k,l) + wgt*(Vmm(i,k)*Wmm(i+j,k,l) &
               + Vmm(i+j,k)*Wmm(i,k,l))
          end do
        end do
      end do
    end do
  end subroutine covmat_accum_th

  subroutine covmat_accum_mb (GamWW, GamVW, GamVV, WW, VV, NNratio, &
     Nout,nnew,nf,nb)
    implicit none
    integer, intent(in) :: Nout, nnew, nf, nb
    double precision, intent(in) :: WW(Nout,nnew,nf), VV(Nout,nnew), &
       NNratio
    double precision, intent(inout) :: GamWW(nnew,nf), GamVW(nnew,nf), &
       GamVV(nnew)
    integer i, j, k, l
    double precision wgt, Wm(nnew,nf), Vm(nnew), Wmm(Nout,nnew,nf), &
       Vmm(Nout,nnew)
    Vm = sum(VV,1)/dble(Nout)
    Wm = sum(WW,1)/dble(Nout)
    do k = 1, nf
      do j = 1, nnew
        do i = 1, Nout
          Wmm(i,j,k) = WW(i,j,k) - Wm(j,k)
        end do
      end do
    end do
    do j = 1, nnew
      do i = 1, Nout
        Vmm(i,j) = VV(i,j) - Vm(j)
      end do
    end do
    ! j = 0
    wgt = NNratio/dble(Nout)
    do l = 1, nf
      do k = 1, nnew
        do i = 1, Nout
          GamWW(k,l) = GamWW(k,l) + wgt*Wmm(i,k,l)*Wmm(i,k,l)
        end do
      end do
    end do
    do k = 1, nnew
      do i = 1, Nout
        GamVV(k) = GamVV(k) + wgt*Vmm(i,k)*Vmm(i,k)
      end do
    end do
    do l = 1, nf
      do k = 1, nnew
        do i = 1, Nout
          GamVW(k,l) = GamVW(k,l) + wgt*Vmm(i,k)*Wmm(i,k,l)
        end do
      end do
    end do
    do j = 1, nb - 1
      wgt = 2*NNratio*(1d0 - dble(j)/nb)/Nout ! w/ a factor because only
      ! one term is computed
      do l = 1, nf
        do k = 1, nnew
          do i = 1, Nout-j
            GamWW(k,l) = GamWW(k,l) + wgt*Wmm(i,k,l)*Wmm(i+j,k,l)
          end do
        end do
      end do
      do k = 1, nnew
        do i = 1, Nout-j
          GamVV(k) = GamVV(k) + wgt*Vmm(i,k)*Vmm(i+j,k)
        end do
      end do
      wgt = .5d0*wgt
      do l = 1, nf
        do k = 1, nnew
          do i = 1, Nout-j
            GamVW(k,l) = GamVW(k,l) + wgt*(Vmm(i,k)*Wmm(i+j,k,l) &
               + Vmm(i+j,k)*Wmm(i,k,l))
          end do
        end do
      end do
    end do
  end subroutine covmat_accum_mb
  
  
  double precision function mean(x, m)
    use interfaces, only: isfinite
    implicit none
    integer, intent(in) :: m
    double precision, intent(in) :: x(m)
    integer i
    double precision tot, tot0, mm
    mean = 0d0
    tot = 0d0
    mm = dble(m)
    do i = 1, m
      tot0 = tot
      tot = tot + x(i)
      if (isfinite(tot) .eq. 0) then
        mean = mean + tot0/mm
        tot = x(i)
      end if
    end do
    mean = mean + tot/mm
  end function mean
end subroutine bfsecalcef

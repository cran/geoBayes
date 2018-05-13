!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Commentary:
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcb_no_st (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, zsample, weights, n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr)

  use modelfcns, condyz_sp => condyz
  use interfaces
  use flogsumexp
  use covfun
  use condyz, only: condyz_gt
  use betaprior
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf, itr(n)
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), zsample(n, Ntot), weights(Ntot), &
     betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n, p), dm(n, n)
  double precision, intent(out) :: bfact(n_nu, n_cov)
  logical lmxi
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh
  integer i, j, k
  double precision logfy(n_nu, Ntot), lfz
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, llikw(n_nu, Ntot), xi(n)

  call create_model(ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  ! Calculate log f(y|z,nu)
  select case (ifam)
  case (0)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_gt(n, y, l, zsample(:, j), nu(k), tsqdfsc, &
           respdfh)
      end do
      call rchkusr
    end do
  case default
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_sp(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  end select

  do k = 1, n_cov
    call calc_cov (phi(k),nsq(k),dm,F,betQ0,&
       kappa(k),n,p,T,TiF,FTF,Ups,ldh_Ups)
    do j = 1, Ntot
      ! Calculate unnormalised log-likelihood at sampled points
      lfz = logpdfz(n, zsample(:, j), Ups, ldh_Ups, xi, lmxi, &
         ssqdfsc, modeldfh)
      do i = 1, n_nu
        llikw(i, j) = logfy(i, j) + lfz - weights(j)
      end do
    end do
    bfact(:, k) = logrsumexp(llikw, n_nu, Ntot)
  end do
end subroutine calcb_no_st


subroutine calcb_wo_st (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, wsample, weights, n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr)

  use modelfcns
  use flogsumexp
  use covfun
  use betaprior
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf, itr(n)
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), wsample(n,Ntot), weights(Ntot), &
     betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n,p), dm(n,n)
  double precision, intent(out) :: bfact(n_nu,n_cov)
  logical lmxi
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh, zsample(n)
  integer i, j, k, m
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, llikw(n_nu,Ntot), xi(n), lfw

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  select case (ifam)
  case (0)
    call rexit ("This method has not been implemented.")
  case default
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         kappa(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          zsample = transfw(wsample(:,j), nu(k))
          lfw = jointyz(n, zsample, y, l, Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, tsq, modeldfh)
          do m = 1, n
            lfw = lfw - logitrwdz(zsample(m),nu(k))
          end do
          llikw(k,j) = lfw - weights(j)
        end do
      end do
      bfact(:, i) = logrsumexp(llikw, n_nu, Ntot)
    end do
  end select
end subroutine calcb_wo_st


subroutine calcb_mu_st (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, musample, weights, n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr)

  use modelfcns
  use interfaces
  use flogsumexp
  use covfun
  use pdfmu, only: logpdfmu_ga
  use betaprior
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf, itr(n)
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), musample(n, Ntot), weights(Ntot), &
     betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n, p), dm(n, n)
  double precision, intent(out) :: bfact(n_nu, n_cov)
  logical lmxi
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh
  integer i, j, k
  double precision lfmu
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, llikw(n_nu, Ntot), xi(n)

  call create_model(ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  select case (ifam)
  case (0)
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         kappa(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_ga(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw(k, j) = lfmu - weights(j)
        end do
      end do
      bfact(:, i) = logrsumexp(llikw, n_nu, Ntot)
    end do
  case default
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         kappa(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu(n, musample(:,j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw(k,j) = lfmu - weights(j)
        end do
      end do
      bfact(:,i) = logrsumexp(llikw, n_nu, Ntot)
    end do
  end select
end subroutine calcb_mu_st


subroutine calcb_tr_st (bfact, philist, nulist, nsqlist, kappalist, &
   icf, n_cov, n_nu, Ntot, sample, weights, n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr)

  use modelfcns, logpdfzf => logpdfz, condymu_mf => condymu
  use flogsumexp
  use covfun
  use betaprior
  use condymu, only: condymu_gt
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf, itr(n)
  double precision, intent(in) :: philist(n_cov), nsqlist(n_cov), &
     kappalist(n_cov), nulist(n_nu), sample(n,Ntot), weights(Ntot), &
     betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n,p), dm(n,n)
  double precision, intent(out) :: bfact(n_nu,n_cov)
  logical lmxi
  double precision ssqdfsc, respdfh, modeldfh, tsqval
  integer i, j, k
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, lglk(Ntot), xi(n)
  double precision nu, phi, nsq, kappa
  double precision zsam(n), msam(n), jsam(n), sam(n)

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  do i = 1, n_cov
    phi = philist(i)
    nsq = nsqlist(i)
    kappa = kappalist(i)
    call calc_cov (phi,nsq,dm,F,betQ0,&
       kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
    do j = 1, n_nu
      nu = nulist(j)
      do k = 1, Ntot
        call rchkusr
        sam = sample(:,k)
        where (itr == 0)
          zsam = sam
          msam = invlink(zsam,nu)
          jsam = 0d0
        elsewhere (itr == 1)
          msam = sam
          zsam = flink(msam,nu)
          jsam = logilinkdz(zsam,nu)
        elsewhere (itr == 2)
          zsam = transfw(sam,nu)
          msam = invlink(zsam,nu)
          jsam = logitrwdz(zsam,nu)
        end where
        lglk(k) = logpdfzf(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh) &
           + condymuf(ifam,n,y,l,msam,tsqval,respdfh) - sum(jsam) - weights(k)
      end do
      bfact(j,i) = logsumexpv(lglk,Ntot)
    end do
  end do

contains

  pure function condymu_sp (n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymu_sp
    condymu_sp = condymu_mf(n,y1,y2,mu,tsqval)
  end function condymu_sp

  pure function condymuf (ifam, n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: n, ifam
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymuf
    select case (ifam)
    case (0)
      condymuf = condymu_gt(n, y1, y2, mu, tsqval, respdfh)
    case default
      condymuf = condymu_sp(n, y1, y2, mu, tsqval, respdfh)
    end select
  end function condymuf
end subroutine calcb_tr_st



!!!!!!!!!!!!!!!!!!!!!!!!! Control variates method !!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcb_no_cv (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, zsample, weights, QRin, n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr)

  use modelfcns, condyz_sp => condyz
  use interfaces
  use flogsumexp
  use covfun
  use condyz, only: condyz_gt
  use betaprior
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf, itr(n)
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), zsample(n, Ntot), weights(Ntot), QRin(Ntot), &
     betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n, p), dm(n, n)
  double precision, intent(out) :: bfact(n_nu, n_cov)
  logical lmxi
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh
  integer i, j, k
  double precision logfy(n_nu, Ntot), lfz, dNtot
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ycv(n_nu, Ntot), llikw, xi(n), betareg(n_nu)

  call create_model(ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)
  dNtot = log(dble(Ntot))

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr
  ! Calculate log f(y|z,nu)
  select case (ifam)
  case (0)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_gt(n, y, l, zsample(:, j), nu(k), tsqdfsc, &
           respdfh)
      end do
      call rchkusr
    end do
  case default
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_sp(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  end select

  do k = 1, n_cov
    call rchkusr
    call calc_cov (phi(k),nsq(k),dm,F,betQ0,&
       kappa(k),n,p,T,TiF,FTF,Ups,ldh_Ups)
    do j = 1, Ntot
      ! Calculate unnormalised log-likelihood at sampled points
      lfz = logpdfz(n, zsample(:,j), Ups, ldh_Ups, xi, lmxi, &
         ssqdfsc, modeldfh)
      do i = 1, n_nu
        llikw = logfy(i,j) + lfz - weights(j)
        ycv(i,j) = exp(llikw + dNtot)
      end do
    end do
    betareg = matmul(ycv,QRin)
    where (betareg .gt. 0d0)
      bfact(:,k) = log(betareg)
    elsewhere
      bfact(:,k) = -huge(1d0)
    end where
  end do
end subroutine calcb_no_cv




subroutine calcb_wo_cv (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, wsample, weights, QRin, n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr)

  use modelfcns
  use flogsumexp
  use covfun
  use betaprior
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf, itr(n)
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), wsample(n,Ntot), weights(Ntot), QRin(Ntot), &
     betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n,p), dm(n,n)
  double precision, intent(out) :: bfact(n_nu,n_cov)
  logical lmxi
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh
  integer i, j, k, m
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, llikw, ycv(n_nu,Ntot), xi(n), lfw, zsample(n), betareg(n_nu)
  double precision dNtot

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)
  dNtot = log(dble(Ntot))

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  select case (ifam)
  case (0)
    call rexit ("This method has not been implemented.")
  case default
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         kappa(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          zsample = transfw(wsample(:,j),nu(k))
          lfw = jointyz(n, zsample, y, l, Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, tsq, modeldfh)
          do m = 1, n
            lfw = lfw - logitrwdz(zsample(m),nu(k))
          end do
          llikw = lfw - weights(j)
          ycv(k,j) = exp(llikw + dNtot)
        end do
      end do
      betareg = matmul(ycv,QRin)
      where (betareg .gt. 0d0)
        bfact(:,i) = log(betareg)
      elsewhere
        bfact(:,i) = -huge(1d0)
      end where
    end do
  end select
end subroutine calcb_wo_cv



subroutine calcb_mu_cv (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, musample, weights, QRin, n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr)

  use modelfcns
  use interfaces
  use flogsumexp
  use covfun
  use pdfmu, only: logpdfmu_ga
  use betaprior
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf, itr(n)
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), musample(n,Ntot), weights(Ntot), QRin(Ntot), &
     betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n, p), dm(n, n)
  double precision, intent(out) :: bfact(n_nu, n_cov)
  logical lmxi
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh
  integer i, j, k
  double precision lfmu, dNtot
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ycv(n_nu, Ntot), llikw, xi(n), betareg(n_nu)

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)
  dNtot = log(dble(Ntot))

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  select case (ifam)
  case (0)
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         kappa(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_ga(n, musample(:,j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw = lfmu - weights(j)
          ycv(k, j) = exp(llikw + dNtot)
        end do
      end do
      betareg = matmul(ycv,QRin)
      where (betareg .gt. 0d0)
        bfact(:,i) = log(betareg)
      elsewhere
        bfact(:,i) = -huge(1d0)
      end where
    end do
  case default
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         kappa(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu(n, musample(:,j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw = lfmu - weights(j)
          ycv(k, j) = exp(llikw + dNtot)
        end do
      end do
      betareg = matmul(ycv,QRin)
      where (betareg .gt. 0d0)
        bfact(:,i) = log(betareg)
      elsewhere
        bfact(:,i) = -huge(1d0)
      end where
    end do
  end select
end subroutine calcb_mu_cv



subroutine calcb_tr_cv (bfact, philist, nulist, nsqlist, kappalist, &
   icf, n_cov, n_nu, Ntot, sample, weights, QRin, &
   n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr)

  use modelfcns, logpdfzf => logpdfz, condymu_mf => condymu
  use flogsumexp
  use covfun
  use betaprior
  use condymu, only: condymu_gt
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf, itr(n)
  double precision, intent(in) :: philist(n_cov), nsqlist(n_cov), &
     kappalist(n_cov), nulist(n_nu), sample(n,Ntot), weights(Ntot), &
     betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n,p), dm(n,n)
  double precision, intent(in) :: QRin(Ntot)
  double precision, intent(out) :: bfact(n_nu,n_cov)
  logical lmxi
  double precision ssqdfsc, respdfh, modeldfh, tsqval
  integer i, j, k
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, lglk(Ntot), xi(n)
  double precision nu, phi, nsq, kappa
  double precision zsam(n), msam(n), jsam(n), sam(n)
  double precision dNtot, ycv(Ntot), betareg

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select
  dNtot = log(dble(Ntot))

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  do i = 1, n_cov
    phi = philist(i)
    nsq = nsqlist(i)
    kappa = kappalist(i)
    call calc_cov (phi,nsq,dm,F,betQ0,&
       kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
    do j = 1, n_nu
      nu = nulist(j)
      do k = 1, Ntot
        call rchkusr
        sam = sample(:,k)
        where (itr == 0)
          zsam = sam
          msam = invlink(zsam,nu)
          jsam = 0d0
        elsewhere (itr == 1)
          msam = sam
          zsam = flink(msam,nu)
          jsam = logilinkdz(zsam,nu)
        elsewhere (itr == 2)
          zsam = transfw(sam,nu)
          msam = invlink(zsam,nu)
          jsam = logitrwdz(zsam,nu)
        end where
        lglk(k) = logpdfzf(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh) &
           + condymuf(ifam,n,y,l,msam,tsqval,respdfh) - sum(jsam) - weights(k)
        ycv(k) = exp(lglk(k) + dNtot)
      end do
      betareg = dot_product(ycv,QRin)
      if (betareg .gt. 0d0) then
        bfact(j,i) = log(betareg)
      else
        bfact(j,i) = -huge(1d0)
      end if
    end do
  end do

contains

  pure function condymu_sp (n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymu_sp
    condymu_sp = condymu_mf(n,y1,y2,mu,tsqval)
  end function condymu_sp

  pure function condymuf (ifam, n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: n, ifam
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymuf
    select case (ifam)
    case (0)
      condymuf = condymu_gt(n, y1, y2, mu, tsqval, respdfh)
    case default
      condymuf = condymu_sp(n, y1, y2, mu, tsqval, respdfh)
    end select
  end function condymuf
end subroutine calcb_tr_cv

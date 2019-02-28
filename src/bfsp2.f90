!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Commentary: Compute the Bayes factors using the z sample
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bfsp_no (weights, zcv, logbf, lglk1, lglk2, &
   philist, omglist, nulist, &
   zsample1, Nout1, Ntot1, zsample2, Nout2, Ntot2, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   kappalist, icf, n, p, kg, ifam, imeth, itr)
  use modelfcns, jointyz_sp => jointyz
  use interfaces
  use linkfcns
  use flogsumexp
  use bmargin
  use covfun
  use jointyz, only: jointyz_gt
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     zsample1(n,Ntot1), zsample2(n,Ntot2), y(n), l(n), F(n, p), &
     dm(n,n), betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: logbf(kg), lglk1(Ntot1,kg), &
     lglk2(Ntot2,kg), weights(Ntot2), zcv(Ntot2,kg)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, ssqdfsc, modeldfh, &
     tsqdfsc, respdfh, xi(n), eta(kg), mxlglk, lglketa(Ntot2, kg), nu
  integer i, j, m
  double precision zsam(n)

  call create_model (ifam)
  call create_spcor (icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  select case (ifam)
  case (0)
    do i = 1, kg
      nu = nulist(i)
      call calc_cov (philist(i),omglist(i),dm,F,betQ0,&
         kappalist(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,i) = jointyz_gt(n, zsample1(:,j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,i) = jointyz_gt(n, zsample2(:,j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
      end do
    end do
  case default
    do i = 1, kg
      nu = nulist(i)
      call calc_cov (philist(i),omglist(i),dm,F,betQ0,&
         kappalist(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,i) = jointyz_sp(n, zsample1(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,i) = jointyz_sp(n, zsample2(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  end select
  mxlglk = maxval(lglk1)
  lglk1 = lglk1 - mxlglk
  ! lglk2 = lglk2 - mxlglk XXX Don't do this; may err SE calc different fcn

  ! Use the sample in reverse logistic regression
  eta = log(dble(Nout1))
  select case (imeth)
  case (1)
    call revlogistic (eta,lglk1,kg,Ntot1,Nout1)
  case (2)
    call revlogistic (eta,lglk1,kg,Ntot1,Nout1)
    call mengwong (eta,lglk1,kg,Ntot1,Nout1)
  end select
  ! eta_i = log(N_i) - log(bf_i)
  logbf = log(dble(Nout1)) - eta

  if (Ntot2 .eq. 0) return

  ! Compute weights
  eta = log(dble(Nout2)) - logbf
  lglketa = spread(eta,1,Ntot2) + lglk2
  weights = logrsumexp(lglketa,Ntot2,kg)

  ! Compute control variates
  lglketa = lglketa - spread(weights,2,kg) &
     + spread(log(dble(Ntot2)/dble(Nout2)),1,Ntot2)

  zcv(:,2:kg) = spread(lglketa(:,1),2,kg-1)
  zcv(:,2:kg) = zcv(:,2:kg) - lglketa(:,2:kg)
  where (zcv(:,2:kg) > 0d0)
    zcv(:,2:kg) = -exp(lglketa(:,2:kg) + flogexpm1(zcv(:,2:kg)))
  elsewhere (zcv(:,2:kg) < 0d0)
    zcv(:,2:kg) = exp(lglketa(:,2:kg) + flog1mexp(zcv(:,2:kg)))
  end where
  zcv(:,1) = 1d0
end subroutine bfsp_no






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Commentary: Compute the Bayes factors using the mu sample
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bfsp_mu (weights, zcv, logbf, lglk1, lglk2, &
   philist, omglist, nulist, &
   musample1, Nout1, Ntot1, musample2, Nout2, Ntot2, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   kappalist, icf, n, p, kg, ifam, imeth, itr)

  use modelfcns
  use interfaces
  use flogsumexp
  use bmargin
  use covfun
  use pdfmu, only: logpdfmu_ga
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     musample1(n, Ntot1), musample2(n, Ntot2), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: logbf(kg), lglk1(Ntot1, kg), &
     lglk2(Ntot2, kg), weights(Ntot2), zcv(Ntot2, kg)
  logical lmxi
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, tsqdfsc, respdfh, xi(n), eta(kg), mxlglk, &
     lglketa(Ntot2, kg), nu
  integer i, j

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  select case (ifam)
  case (0)
    do i = 1, kg
      nu = nulist(i)
      call calc_cov (philist(i),omglist(i),dm,F,betQ0,&
         kappalist(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,i) = logpdfmu_ga(n, musample1(:, j), Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,i) = logpdfmu_ga(n, musample2(:, j), Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, modeldfh)
      end do
    end do
  case default
    do i = 1, kg
      nu = nulist(i)
      call calc_cov (philist(i),omglist(i),dm,F,betQ0,&
         kappalist(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,i) = logpdfmu(n, musample1(:, j), Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,i) = logpdfmu(n, musample2(:, j), Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, modeldfh)
      end do
    end do
  end select
  mxlglk = maxval(lglk1)
  lglk1 = lglk1 - mxlglk
  ! lglk2 = lglk2 - mxlglk XXX Don't do this; may err SE calc different fcn

  ! Use the sample in reverse logistic regression
  eta = log(dble(Nout1))
  select case (imeth)
  case (1)
!     open(11,file='eta.txt'); write(11,*) eta; close(11)
!     open(11,file='lglk1.txt'); write(11,*) lglk1; close(11)
!     open(11,file='kg.txt'); write(11,*) kg; close(11)
!     open(11,file='Ntot1.txt'); write(11,*) Ntot1; close(11)
!     open(11,file='Nout1.txt'); write(11,*) Nout1; close(11)
    call revlogistic (eta,lglk1,kg,Ntot1,Nout1)
  case (2)
    call revlogistic (eta,lglk1,kg,Ntot1,Nout1)
    call mengwong (eta,lglk1,kg,Ntot1,Nout1)
  end select
  ! eta_i = log(N_i) - log(bf_i)
  logbf = log(dble(Nout1)) - eta

  if (Ntot2 .eq. 0) return

  ! Compute weights
  eta = log(dble(Nout2)) - logbf
  lglketa = spread(eta,1,Ntot2) + lglk2
  weights = logrsumexp(lglketa,Ntot2,kg)

  ! Compute control variates
  lglketa = lglketa - spread(weights,2,kg) &
     + spread(log(dble(Ntot2)/dble(Nout2)),1,Ntot2)

  zcv(:,2:kg) = spread(lglketa(:,1),2,kg-1)
  zcv(:,2:kg) = zcv(:,2:kg) - lglketa(:,2:kg)
  where (zcv(:,2:kg) > 0d0)
    zcv(:,2:kg) = -exp(lglketa(:,2:kg) + flogexpm1(zcv(:,2:kg)))
  elsewhere (zcv(:,2:kg) < 0d0)
    zcv(:,2:kg) = exp(lglketa(:,2:kg) + flog1mexp(zcv(:,2:kg)))
  end where
  zcv(:,1) = 1d0
end subroutine bfsp_mu



subroutine bfsp_wo (weights, zcv, logbf, lglk1, lglk2, &
   philist, omglist, nulist, &
   sample1, Nout1, Ntot1, sample2, Nout2, Ntot2, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   kappalist, icf, n, p, kg, ifam, imeth, itr)
  use modelfcns, logpdfzf => logpdfz, condymuf => condymu
  use interfaces
  use linkfcns
  use flogsumexp
  use bmargin
  use covfun
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     sample1(n,Ntot1), sample2(n,Ntot2), y(n), l(n), F(n, p), &
     dm(n,n), betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: logbf(kg), lglk1(Ntot1,kg), &
     lglk2(Ntot2,kg), weights(Ntot2), zcv(Ntot2,kg)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, ssqdfsc, modeldfh, &
     xi(n), eta(kg), lglketa(Ntot2, kg)
  integer i, j
  double precision zsam(n), msam(n), jsam(n), sam(n)
  double precision nu, phi, omg, kappa

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  do i = 1, kg
    nu = nulist(i)
    phi = philist(i)
    omg = omglist(i)
    kappa = kappalist(i)
    call calc_cov (phi,omg,dm,F,betQ0,&
       kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
    do j = 1, Ntot1
      call rchkusr
      sam = sample1(:,j)
      zsam = transfw(sam,nu)
      msam = invlink(zsam,nu)
      jsam = loginvtrwdz(zsam,nu)
      lglk1(j,i) = logpdfzf(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh) &
         + condymuf(n,y,l,msam,tsq) - sum(jsam)
    end do
    do j = 1, Ntot2
      call rchkusr
      sam = sample2(:,j)
      zsam = transfw(sam,nu)
      msam = invlink(zsam,nu)
      jsam = loginvtrwdz(zsam,nu)
      lglk2(j,i) = logpdfzf(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh) &
         + condymuf(n,y,l,msam,tsq) - sum(jsam)
    end do
  end do

  ! Use the sample in reverse logistic regression
  eta = log(dble(Nout1))
  select case (imeth)
  case (1)
    call revlogistic (eta,lglk1,kg,Ntot1,Nout1)
  case (2)
    call revlogistic (eta,lglk1,kg,Ntot1,Nout1)
    call mengwong (eta,lglk1,kg,Ntot1,Nout1)
  end select
  ! eta_i = log(N_i) - log(bf_i)
  logbf = log(dble(Nout1)) - eta

  if (Ntot2 .eq. 0) return

  ! Compute weights
  eta = log(dble(Nout2)) - logbf
  lglketa = spread(eta,1,Ntot2) + lglk2
  weights = logrsumexp(lglketa,Ntot2,kg)

  ! Compute control variates
  lglketa = lglketa - spread(weights,2,kg) &
     + spread(log(dble(Ntot2)/dble(Nout2)),1,Ntot2)

  zcv(:,2:kg) = spread(lglketa(:,1),2,kg-1)
  zcv(:,2:kg) = zcv(:,2:kg) - lglketa(:,2:kg)
  where (zcv(:,2:kg) > 0d0)
    zcv(:,2:kg) = -exp(lglketa(:,2:kg) + flogexpm1(zcv(:,2:kg)))
  elsewhere (zcv(:,2:kg) < 0d0)
    zcv(:,2:kg) = exp(lglketa(:,2:kg) + flog1mexp(zcv(:,2:kg)))
  end where
  zcv(:,1) = 1d0
end subroutine bfsp_wo


subroutine bfsp_tr (weights, zcv, logbf, lglk1, lglk2, &
   philist, omglist, nulist, &
   sample1, Nout1, Ntot1, sample2, Nout2, Ntot2, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   kappalist, icf, n, p, kg, ifam, imeth, itr)
  use modelfcns, logpdfzf => logpdfz, condymu_mf => condymu
  use interfaces
  use linkfcns
  use flogsumexp
  use bmargin
  use covfun
  use condymu, only: condymu_gt
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     sample1(n,Ntot1), sample2(n,Ntot2), y(n), l(n), F(n, p), &
     dm(n,n), betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: logbf(kg), lglk1(Ntot1,kg), &
     lglk2(Ntot2,kg), weights(Ntot2), zcv(Ntot2,kg)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, ssqdfsc, modeldfh, &
     tsqval, respdfh, xi(n), eta(kg), lglketa(Ntot2, kg)
  integer i, j
  double precision zsam(n), msam(n), jsam(n), sam(n)
  double precision nu, phi, omg, kappa

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
    !!condymuf => condymu_gt
  case default
    tsqval = tsq
    !!condymuf => condymu_sp
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  do i = 1, kg
    nu = nulist(i)
    phi = philist(i)
    omg = omglist(i)
    kappa = kappalist(i)
    call calc_cov (phi,omg,dm,F,betQ0,&
       kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
    do j = 1, Ntot1
      call rchkusr
      sam = sample1(:,j)
      where (itr == 0)
        zsam = sam
        msam = invlink(zsam,nu)
        jsam = 0d0
      elsewhere (itr == 1)
        msam = sam
        zsam = flink(msam,nu)
        jsam = loginvlinkdz(zsam,nu)
      elsewhere (itr == 2)
        zsam = transfw(sam,nu)
        msam = invlink(zsam,nu)
        jsam = loginvtrwdz(zsam,nu)
      end where
      lglk1(j,i) = logpdfzf(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh) &
         + condymuf(ifam, n,y,l,msam,tsqval,respdfh) - sum(jsam)
    end do
    do j = 1, Ntot2
      call rchkusr
      sam = sample2(:,j)
      where (itr == 0)
        zsam = sam
        msam = invlink(zsam,nu)
        jsam = 0d0
      elsewhere (itr == 1)
        msam = sam
        zsam = flink(msam,nu)
        jsam = loginvlinkdz(zsam,nu)
      elsewhere (itr == 2)
        zsam = transfw(sam,nu)
        msam = invlink(zsam,nu)
        jsam = loginvtrwdz(zsam,nu)
      end where
      lglk2(j,i) = logpdfzf(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh) &
         + condymuf(ifam, n,y,l,msam,tsqval,respdfh) - sum(jsam)
    end do
  end do

  ! Use the sample in reverse logistic regression
  eta = log(dble(Nout1))
  select case (imeth)
  case (1)
    call revlogistic (eta,lglk1,kg,Ntot1,Nout1)
  case (2)
    call revlogistic (eta,lglk1,kg,Ntot1,Nout1)
    call mengwong (eta,lglk1,kg,Ntot1,Nout1)
  end select
  ! eta_i = log(N_i) - log(bf_i)
  logbf = log(dble(Nout1)) - eta

  if (Ntot2 .eq. 0) return

  ! Compute weights
  eta = log(dble(Nout2)) - logbf
  lglketa = spread(eta,1,Ntot2) + lglk2
  weights = logrsumexp(lglketa,Ntot2,kg)

  ! Compute control variates
  lglketa = lglketa - spread(weights,2,kg) &
     + spread(log(dble(Ntot2)/dble(Nout2)),1,Ntot2)

  zcv(:,2:kg) = spread(lglketa(:,1),2,kg-1)
  zcv(:,2:kg) = zcv(:,2:kg) - lglketa(:,2:kg)
  where (zcv(:,2:kg) > 0d0)
    zcv(:,2:kg) = -exp(lglketa(:,2:kg) + flogexpm1(zcv(:,2:kg)))
  elsewhere (zcv(:,2:kg) < 0d0)
    zcv(:,2:kg) = exp(lglketa(:,2:kg) + flog1mexp(zcv(:,2:kg)))
  end where
  zcv(:,1) = 1d0

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
    integer, intent(in) :: ifam
    integer, intent(in) :: n
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
end subroutine bfsp_tr

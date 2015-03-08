!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Compute the Bayes factors using the z sample
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bfspz (weights, zcv, logbf, lglk1, lglk2, &
   philist, nsqlist, nulist, &
   zsample1, Nout1, Ntot1, zsample2, Nout2, Ntot2, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   kappalist, icf, n, p, kg, ifam, imeth)
  use interfaces
  use linkfcn
  use flogsumexp
  use bmargin
  use covfun
  use pdfy
  use jointyz
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf
  double precision, intent(in) :: philist(kg), nsqlist(kg), nulist(kg), &
     zsample1(n, Ntot1), zsample2(n, Ntot2), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: logbf(kg), lglk1(Ntot1, kg), &
     lglk2(Ntot2, kg), weights(Ntot2), zcv(Ntot2, kg)
  logical lup(n, n), lmxi
  double precision T(n, n), TiF(n, p), FTF(p, p), TFFT(n, p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, &
     tsqdfsc, respdfh, xi(n), eta(kg), mxlglk, lglketa(Ntot2, kg), nu
  integer i, ii, j

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  j=0
  do i = 1, p
    if (betQ0(i,i) /= 0d0) j = j + 1
  end do
  if (j == 0) then ! Flat prior
    modeldfh = .5d0*(n - p + ssqdf)
    xi = 0d0
    lmxi = .false. 
  else ! Normal prior
    modeldfh = .5d0*(n + ssqdf)
    xi = matmul(F,betm0)
    lmxi = any(xi .ne. 0d0)
  end if

  do i = 1, n
    lup(:i-1,i) = .true.
    lup(i:,i) = .false. 
  end do

  select case (ifam)
  case (0) ! Transformed Gaussian
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointyz_gt(n, zsample1(:,j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointyz_gt(n, zsample2(:,j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
      end do
    end do
  case (1) ! Gaussian
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointyz_ga(n, zsample1(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointyz_ga(n, zsample2(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (2) ! Binomial
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointyz_bi(n, zsample1(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointyz_bi(n, zsample2(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (3) ! Poisson
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointyz_po(n, zsample1(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointyz_po(n, zsample2(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (4) ! Gamma
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointyz_gm(n, zsample1(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointyz_gm(n, zsample2(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (5) ! Binomial Asymmetric
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointyz_ba(n, zsample1(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointyz_ba(n, zsample2(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (6) ! Binomial Asymmetric Decreasing
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointyz_bd(n, zsample1(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointyz_bd(n, zsample2(:,j), &
           y, l, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case default
    call rexit ("Unrecognised family")
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
end subroutine bfspz






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Commentary: Compute the Bayes factors using the mu sample
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bfspmu (weights, zcv, logbf, lglk1, lglk2, &
   philist, nsqlist, nulist, &
   musample1, Nout1, Ntot1, musample2, Nout2, Ntot2, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   kappalist, icf, n, p, kg, ifam, imeth)
  use interfaces
  use flogsumexp
  use bmargin
  use covfun
  use jointymu
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf
  double precision, intent(in) :: philist(kg), nsqlist(kg), nulist(kg), &
     musample1(n, Ntot1), musample2(n, Ntot2), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: logbf(kg), lglk1(Ntot1, kg), &
     lglk2(Ntot2, kg), weights(Ntot2), zcv(Ntot2, kg)
  logical lup(n, n), lmxi
  double precision T(n, n), TiF(n, p), FTF(p, p), TFFT(n, p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, tsqdfsc, respdfh, xi(n), eta(kg), mxlglk, &
     lglketa(Ntot2, kg), nu
     
  integer i, ii, j

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  j=0
  do i = 1, p
    if (betQ0(i,i) /= 0d0) j = j + 1
  end do
  if (j == 0) then ! Flat prior
    modeldfh = .5d0*(n - p + ssqdf)
    lmxi = .false.
  else ! Normal prior
    modeldfh = .5d0*(n + ssqdf)
    xi = matmul(F,betm0)
    lmxi = any(xi .ne. 0d0)
  end if

  do i = 1, n
    lup(:i-1,i) = .true.
    lup(i:,i) = .false. 
  end do

  select case (ifam)
  case (0) ! Transformed Gaussian
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointymu_gt(n, musample1(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointymu_gt(n, musample2(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
      end do
    end do
  case (1) ! Gaussian
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointymu_ga(n, musample1(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointymu_ga(n, musample2(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (2) ! Binomial
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointymu_bi(n, musample1(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointymu_bi(n, musample2(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (3) ! Poisson
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointymu_po(n, musample1(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointymu_po(n, musample2(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (4) ! Gamma
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointymu_gm(n, musample1(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointymu_gm(n, musample2(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (5) ! Binomial Asymmetric
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointymu_ba(n, musample1(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointymu_ba(n, musample2(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (6) ! Binomial Asymmetric Decreasing
    do ii = 1, kg
      nu = nulist(ii)
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        lglk1(j,ii) = jointymu_bd(n, musample1(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
      do j = 1, Ntot2
        call rchkusr
        lglk2(j,ii) = jointymu_bd(n, musample2(:, j), y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case default
    call rexit ("Unrecognised family")
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
end subroutine bfspmu

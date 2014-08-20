!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Compute Bayes factors
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bfsp (weights, zcv, logbf, lglk1, lglk2, &
   philist, nsqlist, nulist, &
   zsample1, Nout1, Ntot1, zsample2, Nout2, Ntot2, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   kappalist, icf, n, p, kg, ifam, imeth)
  use interfaces
  use linkfcn
  use flogsumexp
  use bmargin
  use covfun
  implicit none
  integer, intent(in) :: n, p, kg, ifam, imeth, Nout1(kg), Ntot1, &
     Nout2(kg), Ntot2, icf
  double precision, intent(in) :: philist(kg), nsqlist(kg), nulist(kg), &
     zsample1(n, Ntot1), zsample2(n, Ntot2), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: logbf(kg), lglk1(Ntot1, kg), &
     lglk2(Ntot2, kg), weights(Ntot2), zcv(Ntot2, kg)
  logical lup(n, n)
  double precision T(n, n), TiF(n, p), FTF(p, p), TFFT(n, p), Ups(n, n), &
     ldh_Ups, zmxi1(n, Ntot1), zmxi2(n, Ntot2), ssqdfsc, modeldfh, &
     tsqdfsc, respdfh, xi(n), eta(kg), lfy, mu, mxlglk, Upsz(n), zUz, &
     lglketa(Ntot2, kg)
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
    zmxi1 = zsample1
    zmxi2 = zsample2
  else ! Normal prior
    modeldfh = .5d0*(n + ssqdf)
    xi = matmul(F,betm0)
    zmxi1 = spread(xi, 2, Ntot1)
    zmxi1 = zsample1 - zmxi1
    zmxi2 = spread(xi, 2, Ntot2)
    zmxi2 = zsample2 - zmxi2
  end if

  lup = .false.
  do i = 2, n
    lup(:i-1,i) = .true.
  end do

  select case (ifam)
  case (0) ! Transformed Gaussian
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        call dsymv ('u',n,1d0,Ups,n,zmxi1(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
        zUz = dot_product(zmxi1(:,j),Upsz) + ssqdfsc
        lfy = tsqdfsc
        do i = 1, n
          mu = invlink_ga(zsample1(i,j),nulist(ii))
          mu = y(i) - mu
          lfy = lfy + l(i)*mu*mu
        end do
        lglk1(j,ii) = ldh_Ups - modeldfh*log(zUz) - respdfh*log(lfy)
      end do
      do j = 1, Ntot2
        call rchkusr
        call dsymv ('u',n,1d0,Ups,n,zmxi2(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
        zUz = dot_product(zmxi2(:,j),Upsz) + ssqdfsc
        lfy = tsqdfsc
        do i = 1, n
          mu = invlink_ga(zsample2(i,j),nulist(ii))
          mu = y(i) - mu
          lfy = lfy + l(i)*mu*mu
        end do
        lglk2(j,ii) = ldh_Ups - modeldfh*log(zUz) - respdfh*log(lfy)
      end do
    end do
  case (1) ! Gaussian
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        call dsymv ('u',n,1d0,Ups,n,zmxi1(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
        zUz = dot_product(zmxi1(:,j),Upsz) + ssqdfsc
        lfy = 0d0
        do i = 1, n
          mu = invlink_ga(zsample1(i,j),nulist(ii))
          lfy = lfy + y(i)*mu - .5d0*l(i)*mu*mu
        end do
        lglk1(j,ii) = ldh_Ups - modeldfh*log(zUz) + lfy/tsq
      end do
      do j = 1, Ntot2
        call rchkusr
        call dsymv ('u',n,1d0,Ups,n,zmxi2(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
        zUz = dot_product(zmxi2(:,j),Upsz) + ssqdfsc
        lfy = 0d0
        do i = 1, n
          mu = invlink_ga(zsample2(i,j),nulist(ii))
          lfy = lfy + y(i)*mu - .5d0*l(i)*mu*mu
        end do
        lglk2(j,ii) = ldh_Ups - modeldfh*log(zUz) + lfy/tsq
      end do
    end do
  case (2) ! Binomial
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        call dsymv ('u',n,1d0,Ups,n,zmxi1(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
        zUz = dot_product(zmxi1(:,j),Upsz) + ssqdfsc
        lfy = 0d0
        do i = 1, n
          mu = invlink_bi(zsample1(i,j),nulist(ii))
          lfy = lfy + y(i)*mu + l(i)*flog1mexp(mu)
        end do
        lglk1(j,ii) = ldh_Ups - modeldfh*log(zUz) + lfy/tsq
      end do
      do j = 1, Ntot2
        call rchkusr
        call dsymv ('u',n,1d0,Ups,n,zmxi2(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
        zUz = dot_product(zmxi2(:,j),Upsz) + ssqdfsc
        lfy = 0d0
        do i = 1, n
          mu = invlink_bi(zsample2(i,j),nulist(ii))
          lfy = lfy + y(i)*mu + l(i)*flog1mexp(mu)
        end do
        lglk2(j,ii) = ldh_Ups - modeldfh*log(zUz) + lfy/tsq
      end do
    end do
  case (3) ! Poisson
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        call dsymv ('u',n,1d0,Ups,n,zmxi1(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
        zUz = dot_product(zmxi1(:,j),Upsz) + ssqdfsc
        lfy = 0d0
        do i = 1, n
          mu = invlink_po(zsample1(i,j),nulist(ii))
          lfy = lfy + y(i)*mu - l(i)*exp(mu)
        end do
        lglk1(j,ii) = ldh_Ups - modeldfh*log(zUz) + lfy/tsq
      end do
      do j = 1, Ntot2
        call rchkusr
        call dsymv ('u',n,1d0,Ups,n,zmxi2(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
        zUz = dot_product(zmxi2(:,j),Upsz) + ssqdfsc
        lfy = 0d0
        do i = 1, n
          mu = invlink_po(zsample2(i,j),nulist(ii))
          lfy = lfy + y(i)*mu - l(i)*exp(mu)
        end do
        lglk2(j,ii) = ldh_Ups - modeldfh*log(zUz) + lfy/tsq
      end do
    end do
  case (4) ! Gamma
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
      do j = 1, Ntot1
        call rchkusr
        call dsymv ('u',n,1d0,Ups,n,zmxi1(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
        zUz = dot_product(zmxi1(:,j),Upsz) + ssqdfsc
        lfy = 0d0
        do i = 1, n
          mu = invlink_gm(zsample1(i,j),nulist(ii))
          lfy = lfy + y(i)*mu - l(i)*log(-mu)
        end do
        lglk1(j,ii) = ldh_Ups - modeldfh*log(zUz) + lfy/tsq
      end do
      do j = 1, Ntot2
        call rchkusr
        call dsymv ('u',n,1d0,Ups,n,zmxi2(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
        zUz = dot_product(zmxi2(:,j),Upsz) + ssqdfsc
        lfy = 0d0
        do i = 1, n
          mu = invlink_gm(zsample2(i,j),nulist(ii))
          lfy = lfy + y(i)*mu - l(i)*log(-mu)
        end do
        lglk2(j,ii) = ldh_Ups - modeldfh*log(zUz) + lfy/tsq
      end do
    end do
  end select
  mxlglk = maxval(lglk1)
  lglk1 = lglk1 - mxlglk
  ! lglk2 = lglk2 - mxlglk
  
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
end subroutine bfsp

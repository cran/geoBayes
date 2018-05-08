module mcmcfcns
  implicit none
  private
  double precision, parameter :: bigpos = huge(1d0), bigneg = -bigpos
  public ini_mcmc, end_mcmc, sample_ssq, sample_tsq, sample_beta, &
     sample_cov, sample_z, sample_z0, samplez_gt

contains
  subroutine ini_mcmc (lglk, z, p0, phi, nsq, y1, y2, F, kappa, icf, dm, &
     betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, dft, n, p, ifam, &
     betQm0, zmxi, T, TiF, FTF, Ups, Upsz, zUz, ldh_Ups, &
     modeldfh, ssqdfsc, respdf, tsqdfsc, tsqyy, lnewcov)
    ! tsq is either tsq or tsqsc
    use interfaces
    use covfun
    use betaprior
    use modelfcns
    use linkfcns, only: invlink_ga
    use pdfy, only: logpdfy_gt
    implicit none
    integer, intent(in) :: n, p, ifam, icf
    double precision, intent(in) :: y1(n), y2(n), dm(n,n), &
       F(n,p), kappa, betm0(p), betQ0(p,p), ssqdf, ssqsc, &
       tsqdf, tsq, phi, nsq, dft
    double precision, intent(out) :: zUz, ldh_Ups, modeldfh, ssqdfsc, &
       respdf, tsqdfsc, tsqyy
    double precision, intent(inout) :: z(n)
    logical, intent(out) :: lnewcov
    double precision, intent(out) :: lglk, p0(n), betQm0(p), zmxi(n), &
       T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), Upsz(n)
    logical lmxi
    integer i

    call create_spcor(icf,n)

    ! Starting value for z
    z = mustart(y1,y2)
    z = flink(z,dft)

    ssqdfsc = ssqdf*ssqsc

    ! Determine flat or normal prior
    call betapriorz (modeldfh, zmxi, lmxi, betm0, betQ0, F, n, p, ssqdf)
    if (lmxi) then
      zmxi = z - zmxi
    else
      zmxi = z
    end if

    if (betQ0(1,1) .gt. 0d0) then
      call dsymv ('u',p,1d0,betQ0,p,betm0,1,0d0,betQm0,1) ! betQm0 = Q0*m0
    else
      betQm0 = 0d0
    end if

    call calc_cov (phi,nsq,dm,F,betQ0, &
       kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)

    call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
    zUz = dot_product(zmxi,Upsz) + ssqdfsc

    lnewcov = .true. ! Indicates whether to compute prediction matrix

    lglk = ldh_Ups - modeldfh*log(zUz)

    select case (ifam)
    case (0)
      tsqdfsc = tsqdf*tsq
      respdf = n + tsqdf
      tsqyy = tsqdfsc
      do i = 1, n
        p0(i) = invlink_ga(z(i),dft)
        tsqyy = tsqyy + logpdfy_gt(y1(i), y2(i), p0(i))
      end do
      lglk = lglk - .5d0*respdf*log(tsqyy)
    case default
      tsqyy = 0d0
      do i = 1, n
        p0(i) = invlink(z(i),dft)
        tsqyy = tsqyy + logpdfy(y1(i), y2(i), p0(i))
      end do
      lglk = lglk + tsqyy/tsq
    end select

    call rngini
  end subroutine ini_mcmc

  subroutine end_mcmc
    use interfaces, only: rngend
    call rngend
  end subroutine end_mcmc


  subroutine sample_ssq (ssq, modeldfh, zUz)
    use interfaces, only: randgamma
    implicit none
    double precision, intent(in) :: zUz, modeldfh
    double precision, intent(out) :: ssq
    ssq = randgamma(modeldfh)
    ssq = zUz/(ssq+ssq)
  end subroutine sample_ssq

  subroutine sample_tsq (tsq, respdf, tsqyy)
    use interfaces, only: randchisq
    implicit none
    double precision, intent(in) :: respdf, tsqyy
    double precision, intent(out) :: tsq
    tsq = randchisq(respdf)
    tsq = tsqyy/tsq
  end subroutine sample_tsq

  subroutine sample_beta (beta, z, ssq, n, p, betQm0, TiF, FTF)
    use interfaces, only: randnorm
    implicit none
    integer, intent(in) :: n, p
    double precision, intent(in) :: z(n), betQm0(p), TiF(n,p), FTF(p,p)
    double precision, intent(inout) :: beta(p), ssq
    integer j
    beta = betQm0
    call dgemv ('t',n,p,1d0,TiF,n,z,1,1d0,beta,1)
    call dtrmv ('u','t','n',p,FTF,p,beta,1)
    do j = 1, p
      beta(j) = sqrt(ssq)*randnorm() + beta(j)
    end do
    call dtrmv ('u','n','n',p,FTF,p,beta,1)
  end subroutine sample_beta

  subroutine sample_cov (lglk, phi, nsq, phipars, nsqpars, phisc, nsqsc, dm, &
     F, betQ0, n, p, kappa, icf, acc, zmxi, T, TiF, FTF, Ups, &
     Upsz, lnewcov, zUz, ldh_Ups, modeldfh, ssqdfsc)
    use interfaces
    use covfun
    implicit none
    integer, intent(in) :: n, p, icf
    logical, intent(inout) :: lnewcov
    integer, intent(inout):: acc
    double precision, intent(in) :: modeldfh, dm(n,n), F(n,p), &
       kappa, betQ0(p,p), zmxi(n), ssqdfsc
    double precision, intent(in) :: phipars(4), nsqpars(4), phisc, nsqsc
    double precision, intent(inout) :: phi, nsq, lglk, T(n,n), TiF(n,p), &
       FTF(p,p), Ups(n,n), Upsz(n), zUz, ldh_Ups
    double precision tr, ll, u, phi2, nsq2, up3
    double precision FTF2(p,p), Ups2(n,n), Upsz2(n), &
       zUz2, ldh_Ups2, T2(n,n), TiF2(n,p)

    call create_spcor(icf,n)

    if (phisc .le. 0d0 .and. nsqsc .le. 0d0) return

    tr = 0d0

!! Sample phi
    if (phisc .gt. 0d0) then
      ll = phi - phipars(4)
      u = randnorm()*phisc ! u ~ N(0,phisc^2)
      phi2 = ll*exp(u) + phipars(4) ! log(phi2-p4) ~ N(log(phi-p4), phisc^2)
      ! Add logratio of priors
      up3 = u*phipars(3)
      if (up3 .lt. 0d0) then
        tr = tr + phipars(2)*u + exp(phipars(3)*log(ll/phipars(1)) &
           + flog1mexp(up3))
      else
        tr = tr + phipars(2)*u - exp(phipars(3)*log(ll/phipars(1)) &
           + flogexpm1(up3))
      end if
    else
      phi2 = phi
    end if

!! Sample nsq
    if (nsqsc .gt. 0d0) then
      ll = nsq - nsqpars(4)
      u = randnorm()*nsqsc ! u ~ N(0,nsqsc^2)
      nsq2 = ll*exp(u) + nsqpars(4) ! log(nsq2-p4) ~ N(log(nsq-p4), nsqsc^2)
      ! Add logratio of priors
      up3 = u*nsqpars(3)
      if (up3 .lt. 0d0) then
        tr = tr + nsqpars(2)*u + exp(nsqpars(3)*log(ll/nsqpars(1)) &
           + flog1mexp(up3))
      else
        tr = tr + nsqpars(2)*u - exp(nsqpars(3)*log(ll/nsqpars(1)) &
           + flogexpm1(up3))
      end if
    else
      nsq2 = nsq
    end if

    if (tr .le. bigneg .or. tr .ne. tr) return

!! Compute covariance matrices for the new phi, nsq
    call calc_cov (phi2,nsq2,dm,F,betQ0, &
       kappa,n,p,T2,TiF2,FTF2,Ups2,ldh_Ups2)

!! Compute log likelihood ratio
    call dsymv ('u',n,1d0,Ups2,n,zmxi,1,0d0,Upsz2,1) ! Upsz = Ups*(z-xi)
    zUz2 = dot_product(zmxi,Upsz2) + ssqdfsc
    ll = ldh_Ups2 - ldh_Ups - modeldfh*(log(zUz2) - log(zUz))
!! Add the logratio and transition probabilities
    tr = ll + tr

    u = randunif()
    u = log(u)
    if (u .lt. tr) then ! Update phi, nsq
      lnewcov = .true.
      lglk = lglk + ll
      phi = phi2
      nsq = nsq2
      FTF = FTF2
      T = T2
      TiF = TiF2
      Ups = Ups2
      Upsz = Upsz2
      zUz = zUz2
      ldh_Ups = ldh_Ups2
      acc = acc + 1
    end if
  end subroutine sample_cov

  subroutine sample_z0 (z0,z,beta,ssq,phi,nsq,n0,n,p,dmdm0,F,F0,kappa,icf,&
     T, z0_ups, TC, FCTF, lnewcov)
    use interfaces
    use covfun
    implicit none
    integer, intent(in) :: n, n0, p, icf
    double precision, intent(in) :: z(n), beta(p), ssq, phi, nsq, &
       dmdm0(n,n0), F(n,p), F0(n0,p), kappa, T(n,n)
    double precision, intent(out) :: z0(n0)
    double precision, intent(inout) :: TC(n,n0), FCTF(n0,p), &
       z0_ups(n0)
    logical, intent(inout) :: lnewcov
    integer j
    double precision z0_sd(n0), z0_mean(n0)
    call create_spcor(icf,0)
    if (lnewcov) then
      call calc_cov_pred(z0_ups, TC, FCTF, phi, nsq, dmdm0, F, &
         F0, kappa, T, n, n0, p)
      lnewcov = .false.
    end if
    call dgemv ('t',n,n0,1d0,TC,n,z,1,0d0,z0_mean,1) ! z0_mean = C'*T^{-1}*z
    call dgemv ('n',n0,p,1d0,FCTF,n0,beta,1,1d0,z0_mean,1) ! + (F0-C'T^-1F)*beta
    z0_sd = sqrt(ssq)*z0_ups
    do j = 1, n0
      z0(j) = z0_sd(j)*randnorm() + z0_mean(j)
    end do
  end subroutine sample_z0



  subroutine sample_z (lglk,z,p0,y1,y2,dft,ssq,tsq,zmxi,Ups,Upsz,zUz, &
     modeldfh,n)
    use interfaces
    use modelfcns
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), dft, ssq, tsq, Ups(n,n), &
       modeldfh
    double precision, intent(inout) :: z(n), lglk, p0(n), zmxi(n), Upsz(n),&
       zUz
    integer j
    double precision uj(n), z_mean, u, zz, pp, ll

    do j = 1, n
      uj = (/Ups(1:j,j),Ups(j,j+1:n)/) ! jth column of Ups
      z_mean = z(j) - Upsz(j)/uj(j)
      u = randnorm()
      zz = u*sqrt(ssq/uj(j)) + z_mean
      pp = invlink(zz,dft)
      ll = logdffy(y1(j), y2(j), pp, p0(j))
      ll = ll/tsq
      if (ll .le. bigneg .or. ll .ne. ll) return
      u = randunif()
      u = log(u)
      if (u .lt. ll) then ! Update z(j)
        u = zz - z(j)
        z(j) = zz
        p0(j) = pp
        zmxi(j) = zmxi(j) + u
        Upsz = Upsz + uj*u
        zz = zUz
        zUz = zUz + 2d0*u*Upsz(j) - uj(j)*u*u
        lglk = lglk + ll - modeldfh*(log(zUz) - log(zz))
      end if
    end do
  end subroutine sample_z

!!! Sample GRF with Gaussian response from the transformed model, i.e.
!!! also estimate tsq
!!!
!!! @param z    The GRF sample
!!! @param ym   Average at each location
!!! @param l    Sample size at each location
!!! @param dft  Link parameter
!!! @param ssq  Variance parameter
!!! @param tsq  Dispersion parameter
!!! @param n	Number of locations
  subroutine samplez_gt (lglk,z,p0,ym,l,dft,ssq,zmxi,Ups,Upsz,zUz, &
     modeldfh,respdf,tsqyy,n)
    use interfaces
    use linkfcns, only: invlink_ga
    use pdfy, only: logdffy_gt
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: ym(n), l(n), dft, ssq, Ups(n,n), &
       modeldfh, respdf
    double precision, intent(inout) :: z(n), lglk, p0(n), zmxi(n), Upsz(n),&
       zUz, tsqyy
    integer j
    double precision uj(n), z_mean, u, zz, pp, ll, tsqyy2

    do j = 1, n
      uj = (/Ups(1:j,j),Ups(j,j+1:n)/) ! jth column of Ups
      z_mean = z(j) - Upsz(j)/uj(j)
      u = randnorm()
      zz = u*sqrt(ssq/uj(j)) + z_mean
      pp = invlink_ga(zz,dft) ! Mean of y
      tsqyy2 = tsqyy + logdffy_gt(ym(j), l(j), pp, p0(j))
      ll = -.5d0*respdf*(log(tsqyy2) - log(tsqyy))
      if (ll .le. bigneg .or. ll .ne. ll) return
      u = randunif()
      u = log(u)
      if (u .lt. ll) then ! Update z(j)
        tsqyy = tsqyy2
        u = zz - z(j)
        z(j) = zz
        p0(j) = pp
        zmxi(j) = zmxi(j) + u
        Upsz = Upsz + uj*u
        zz = zUz
        zUz = zUz + 2d0*u*Upsz(j) - uj(j)*u*u
        lglk = lglk + ll - modeldfh*(log(zUz) - log(zz))
      end if
    end do
  end subroutine samplez_gt
end module mcmcfcns

!!! mcsp2.f90 --- 
!! 
!! Author: Evangelos Evangelou
!! Created: Tue, 15 Jul, 2014 16:59 (BST)
!! Last-Updated: Wed, 20 Aug, 2014 19:06 (BST)
!!     Update #: 186
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mcmcfcns
  use linkfcn
  implicit none
  private
  logical lnewcov
  logical, pointer :: lup(:,:)
  double precision, pointer :: p0(:), xi(:), betQm0(:), zmxi(:)
  double precision, pointer :: FTF(:,:), TFFT(:,:), Ups(:,:), Upsz(:)
  double precision, pointer :: T(:,:), TiF(:,:)
  double precision :: zUz, ldh_Ups, modeldf, ssqdfsc, respdf, &
     tsqdfsc, tsqyy 
  double precision, pointer :: TC(:,:), FCTF(:,:), z0_ups(:)
  double precision, parameter :: bigpos = huge(1d0), bigneg = -bigpos
  public ini_mcmc, end_mcmc, sample_ssq, sample_tsq, sample_beta, &
     sample_cov, sample_z0, samplez_bi, samplez_po, samplez_gm, &
     samplez_ga, samplez_gt

contains
  subroutine ini_mcmc (lglk, z, phi, nsq, y1, y2, F, kappa, icf, &
     dm, &
     betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, dft, n, n0, p, ifam)
    ! tsq is either tsq or tsqsc
    use interfaces
    use covfun
    implicit none
    integer, intent(in) :: n, n0, p, ifam, icf
    double precision, intent(in) :: y1(n), y2(n), dm(n,n), &
       F(n,p), kappa, betm0(p), betQ0(p,p), ssqdf, ssqsc, &
       tsqdf, tsq, z(n), phi, nsq, dft
    double precision, intent(out) :: lglk
    integer i, j

    allocate ( lup(n,n), p0(n), xi(n), betQm0(p), zmxi(n) )
    allocate ( T(n,n), TiF(n,p), FTF(p,p), TFFT(n,p), Ups(n,n), Upsz(n) )
    allocate ( TC(n,n0), FCTF(n0,p), Z0_ups(n0) )

    ssqdfsc = ssqdf*ssqsc

    ! Determine flat or normal prior
    j=0
    do i = 1, p
      if (betQ0(i,i) /= 0d0) j = j + 1
    end do
    if (j == 0) then ! Flat prior
      modeldf = n - p + ssqdf
      xi = 0d0
      betQm0 = 0d0
      zmxi = z
    else ! Normal prior
      modeldf = n + ssqdf
      xi = matmul(F,betm0)
      call dsymv ('u',p,1d0,betQ0,p,betm0,1,0d0,betQm0,1) ! betQm0 = Q0*m0
      zmxi = z - xi
    end if

    lup = .false.
    do i = 2, n
      lup(:i-1,i) = .true.
    end do

    call calc_cov (phi,nsq,dm,F,betQ0, &
       lup,kappa,icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)

    call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
    zUz = dot_product(zmxi,Upsz) + ssqdfsc

    lnewcov = .true. ! Indicates whether to compute prediction matrix

    lglk = ldh_Ups - .5d0*modeldf*log(zUz)

    select case (ifam)
    case (0) ! Transformed Gaussian
      tsqdfsc = tsqdf*tsq
      respdf = n + tsqdf
      tsqyy = tsqdfsc
      do i = 1, n
        p0(i) = invlink_ga(z(i),dft)
        p0(i) = y1(i) - p0(i)
        p0(i) = p0(i)*p0(i)
        tsqyy = tsqyy + y2(i)*p0(i)
      end do
      lglk = lglk - .5d0*respdf*log(tsqyy)
    case (1) ! Gaussian
      do i = 1, n
        p0(i) = invlink_ga(z(i),dft)
        lglk = lglk + (y1(i)*p0(i) - y2(i)*.5d0*p0(i)*p0(i))/tsq
      end do
    case (2) ! Binomial
      do i = 1, n
        p0(i) = invlink_bi(z(i),dft)
        lglk = lglk + (y1(i)*p0(i) + y2(i)*flog1mexp(p0(i)))/tsq
      end do
    case (3) ! Poisson
      do i = 1, n
        p0(i) = invlink_po(z(i),dft)
        lglk = lglk + (y1(i)*p0(i) - y2(i)*exp(p0(i)))/tsq
      end do
    case (4) ! Gamma
      do i = 1, n
        p0(i) = invlink_gm(z(i),dft)
        lglk = lglk + (y1(i)*p0(i) - y2(i)*log(-p0(i)))/tsq
      end do
    end select

    call rngini
  end subroutine ini_mcmc

  subroutine end_mcmc
    use interfaces, only: rngend
    deallocate ( lup, p0, xi, betQm0, zmxi )
    deallocate ( T, TiF, FTF, TFFT, Ups, Upsz )
    deallocate ( TC, FCTF, Z0_ups )
    call rngend
  end subroutine end_mcmc


  subroutine sample_ssq (ssq)
    use interfaces, only: randchisq
    implicit none
    double precision, intent(out) :: ssq
    ssq = randchisq(modeldf)
    ssq = zUz/ssq
  end subroutine sample_ssq

  subroutine sample_tsq (tsq)
    use interfaces, only: randchisq
    implicit none
    double precision, intent(out) :: tsq
    tsq = randchisq(respdf)
    tsq = tsqyy/tsq
  end subroutine sample_tsq

  subroutine sample_beta (beta, z, ssq, n, p)
    use interfaces, only: randnorm
    implicit none 
    integer, intent(in) :: n, p
    double precision, intent(in) :: z(n)
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
     F, betQ0, n, p, kappa, icf, acc)
    use interfaces
    use covfun
    implicit none
    integer, intent(in) :: n, p, icf
    integer, intent(inout):: acc
    double precision, intent(in) :: dm(n,n), F(n,p), &
       kappa, betQ0(p,p)
    double precision, intent(in) :: phipars(4), nsqpars(4), phisc, nsqsc
    double precision, intent(inout) :: phi, nsq, lglk
    double precision tr, ll, u, phi2, nsq2
    double precision FTF2(p,p), TFFT2(n,p), Ups2(n,n), Upsz2(n), &
       zUz2, ldh_Ups2, T2(n,n), TiF2(n,p)

    if (phisc .le. 0d0 .and. nsqsc .le. 0d0) return

    tr = 0d0

!! Sample phi
    if (phisc .gt. 0d0) then
      ll = phi - phipars(4)
      u = randnorm()*phisc ! u ~ N(0,phisc^2)
      phi2 = ll*exp(u) + phipars(4) ! log(phi2-p4) ~ N(log(phi-p4), phisc^2)
      ! Add logratio of priors
      if (u .lt. 0d0) then
        tr = tr + phipars(2)*u + exp(phipars(3)*log(ll/phipars(1)) &
           + flog1mexp(u*phipars(3)))
      else
        tr = tr + phipars(2)*u - exp(phipars(3)*log(ll/phipars(1)) &
           + flogexpm1(u*phipars(3)))
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
      if (u .lt. 0d0) then
        tr = tr + nsqpars(2)*u + exp(nsqpars(3)*log(ll/nsqpars(1)) &
           + flog1mexp(u*nsqpars(3)))
      else
        tr = tr + nsqpars(2)*u - exp(nsqpars(3)*log(ll/nsqpars(1)) &
           + flogexpm1(u*nsqpars(3)))
      end if
    else
      nsq2 = nsq
    end if

    if (tr .le. bigneg .or. tr .ne. tr) return

!! Compute covariance matrices for the new phi, nsq
    call calc_cov (phi2,nsq2,dm,F,betQ0, &
       lup,kappa,icf,n,p,T2,TiF2,FTF2,TFFT2,Ups2,ldh_Ups2)

!! Compute log likelihood ratio
    call dsymv ('u',n,1d0,Ups2,n,zmxi,1,0d0,Upsz2,1) ! Upsz = Ups*(z-xi)
    zUz2 = dot_product(zmxi,Upsz2) + ssqdfsc
    ll = ldh_Ups2 - ldh_Ups -.5d0*modeldf*(log(zUz2) - log(zUz))
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
      TFFT = TFFT2
      Ups = Ups2
      Upsz = Upsz2
      zUz = zUz2
      ldh_Ups = ldh_Ups2
      acc = acc + 1
    end if
  end subroutine sample_cov

  subroutine sample_z0 (z0,z,beta,ssq,phi,nsq,n0,n,p,dmdm0,F,F0,kappa,icf)
    use interfaces
    use covfun
    implicit none
    integer, intent(in) :: n, n0, p, icf
    double precision, intent(in) :: z(n), beta(p), ssq, phi, nsq, &
       dmdm0(n,n0), F(n,p), F0(n0,p), kappa
    double precision, intent(out) :: z0(n0)
    integer j
    double precision z0_sd(n0), z0_mean(n0)
    if (n0 .eq. 0) return
    if (lnewcov) then
      call calc_cov_pred(z0_ups, TC, FCTF, phi, nsq, dmdm0, F, &
         F0, kappa, icf, T, n, n0, p)
      lnewcov = .false.
    end if
    call dgemv ('t',n,n0,1d0,TC,n,z,1,0d0,z0_mean,1) ! z0_mean = C'*T^{-1}*z
    call dgemv ('n',n0,p,1d0,FCTF,n0,beta,1,1d0,z0_mean,1) ! + (F0-C'T^-1F)*beta
    z0_sd = sqrt(ssq)*z0_ups
    do j = 1, n0
      z0(j) = z0_sd(j)*randnorm() + z0_mean(j)
    end do
  end subroutine sample_z0


!!! Sample GRF with Binomial response
!!! 
!!! @param z    The GRF sample
!!! @param ys   Number of successes at each location
!!! @param yf   Number of failures
!!! @param dft  Link parameter
!!! @param ssq  Variance parameter
!!! @param tsq  Dispersion parameter
!!! @param n	Number of locations
  subroutine samplez_bi (lglk,z,ys,yf,dft,ssq,tsq,n)
    use interfaces
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: ys(n), yf(n), dft, ssq, tsq
    double precision, intent(inout) :: z(n), lglk
    integer j
    double precision uj(n), z_mean, u, zz, pp, ll

    do j = 1, n
      uj = (/Ups(1:j,j),Ups(j,j+1:n)/) ! jth column of Ups
      z_mean = z(j) - Upsz(j)/uj(j)
      u = randnorm()
      zz = u*sqrt(ssq/uj(j)) + z_mean
      pp = invlink_bi(zz,dft)
      ll = ys(j)*(pp-p0(j)) + yf(j)*(flog1mexp(pp)-flog1mexp(p0(j)))
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
        lglk = lglk + ll - (.5d0*modeldf)*(log(zUz) - log(zz))
      end if
    end do
  end subroutine samplez_bi

!!! Sample GRF with Poisson response
!!! 
!!! @param z    The GRF sample
!!! @param yo   Oserved count at each location
!!! @param yt   Total time of observation (weight)
!!! @param dft  Link parameter
!!! @param ssq  Variance parameter
!!! @param tsq  Dispersion parameter
!!! @param n	Number of locations
  subroutine samplez_po (lglk,z,yo,yt,dft,ssq,tsq,n)
    use interfaces
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: yo(n), yt(n), dft, ssq, tsq
    double precision, intent(inout) :: z(n), lglk
    integer j
    double precision uj(n), z_mean, u, zz, pp, ll

    do j = 1, n
      uj = (/Ups(1:j,j),Ups(j,j+1:n)/) ! jth column of Ups
      z_mean = z(j) - Upsz(j)/uj(j)
      u = randnorm()
      zz = u*sqrt(ssq/uj(j)) + z_mean
      pp = invlink_po(zz,dft)
      ll = yo(j)*(pp-p0(j)) - yt(j)*(exp(pp)-exp(p0(j)))
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
        lglk = lglk + ll - (.5d0*modeldf)*(log(zUz) - log(zz))
      end if
    end do
  end subroutine samplez_po

!!! Sample GRF with Gamma response
!!! 
!!! @param z    The GRF sample
!!! @param ys   Sum of oservation at each location
!!! @param yn   Number of samples at each location (weight)
!!! @param dft  Link parameter
!!! @param ssq  Variance parameter
!!! @param tsq  Dispersion parameter
!!! @param n	Number of locations
  subroutine samplez_gm (lglk,z,ys,yn,dft,ssq,tsq,n)
    use interfaces
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: ys(n), yn(n), dft, ssq, tsq
    double precision, intent(inout) :: z(n), lglk
    integer j
    double precision uj(n), z_mean, u, zz, pp, ll

    do j = 1, n
      uj = (/Ups(1:j,j),Ups(j,j+1:n)/) ! jth column of Ups
      z_mean = z(j) - Upsz(j)/uj(j)
      u = randnorm()
      zz = u*sqrt(ssq/uj(j)) + z_mean
      pp = invlink_gm(zz,dft)
      ll = ys(j)*(pp - p0(j)) + yn(j)*(log(-pp) - log(-p0(j)))
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
        lglk = lglk + ll - (.5d0*modeldf)*(log(zUz) - log(zz))
      end if
    end do
  end subroutine samplez_gm

!!! Sample GRF with Gaussian response
!!! 
!!! @param z    The GRF sample
!!! @param ys   Total at each location
!!! @param yl   Number of samples
!!! @param dft  Link parameter
!!! @param ssq  Variance parameter
!!! @param tsq  Dispersion parameter
!!! @param n	Number of locations
  subroutine samplez_ga (lglk,z,ys,yl,dft,ssq,tsq,n)
    use interfaces
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: ys(n), yl(n), dft, ssq, tsq
    double precision, intent(inout) :: z(n), lglk
    integer j
    double precision uj(n), z_mean, u, zz, pp, ll, dfpp

    do j = 1, n
      uj = (/Ups(1:j,j),Ups(j,j+1:n)/) ! jth column of Ups
      z_mean = z(j) - Upsz(j)/uj(j)
      u = randnorm()
      zz = u*sqrt(ssq/uj(j)) + z_mean
      pp = invlink_ga(zz,dft) ! Mean of y
      dfpp = pp*pp - p0(j)*p0(j)
      ll = ys(j)*(pp-p0(j)) - yl(j)*.5d0*dfpp
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
        lglk = lglk + ll - (.5d0*modeldf)*(log(zUz) - log(zz))
      end if
    end do
  end subroutine samplez_ga

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
  subroutine samplez_gt (lglk,z,ym,l,dft,ssq,tsq,n)
    use interfaces
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: ym(n), l(n), dft, ssq, tsq
    double precision, intent(inout) :: z(n), lglk
    integer j
    double precision uj(n), z_mean, u, zz, pp, ll, lli, tsqyy2

    do j = 1, n
      uj = (/Ups(1:j,j),Ups(j,j+1:n)/) ! jth column of Ups
      z_mean = z(j) - Upsz(j)/uj(j)
      u = randnorm()
      zz = u*sqrt(ssq/uj(j)) + z_mean
      pp = invlink_ga(zz,dft) ! Mean of y
      pp = ym(j) - pp
      pp = pp*pp
      tsqyy2 = tsqyy + l(j)*(pp - p0(j))
      lli = -.5d0*respdf*(log(tsqyy2) - log(tsqyy))
      ll = -.5d0/tsq*(pp - p0(j))
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
        lglk = lglk + lli - (.5d0*modeldf)*(log(zUz) - log(zz))
      end if
    end do
  end subroutine samplez_gt
end module mcmcfcns



subroutine mcspsample (lglk, z, z0, beta, ssq, phi, nsq, acc, y, l, F, F0, &
   betm0, betQ0, ssqdf, ssqsc, phipars, phisc, nsqpars, nsqsc, kappa, icf, &
   dft, tsq, dm, dmdm0, Nout, Nbi, Nthin, n, n0, p, ifam)

  use mcmcfcns
  implicit none 
  integer, intent(in) :: Nout, Nbi, Nthin, n, n0, p, ifam, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), F0(n0,p), &
     dft, tsq, dm(n,n), dmdm0(n,n0), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, phipars(4), nsqpars(4), kappa, phisc, nsqsc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), nsq(Nout)
  double precision, intent(out) :: z0(n0,Nout), beta(p,Nout), ssq(Nout), &
     lglk(Nout)
  integer, intent(out) :: acc
  integer i, j
  double precision, parameter :: tsqdf = 0d0

  acc = 0
  i = 1
  call ini_mcmc(lglk(i),z(:,i),phi(i),nsq(i),y,l,F,kappa,icf,dm,&
     betm0,betQ0,ssqdf,ssqsc,tsqdf,tsq,dft,n,n0,p,ifam)
  call rchkusr
  select case (ifam) ! Which family?
  case (1) ! Gaussian
    do j = 1, max(Nbi,Nthin)
      call sample_cov(lglk(i),phi(i),nsq(i),phipars,nsqpars,phisc,nsqsc,&
         dm,F,betQ0,n,p,kappa,icf,acc)
      call sample_ssq(ssq(i))
      call samplez_ga(lglk(i),z(:,i),y,l,dft,ssq(i),tsq,n)
    end do
    call sample_beta(beta(:,i),z(:,i),ssq(i),n,p)
    call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),nsq(i),&
       n0,n,p,dmdm0,F,F0,kappa,icf)
    call rchkusr
    do i = 2, Nout
      lglk(i) = lglk(i-1)
      z(:,i) = z(:,i-1)
      phi(i) = phi(i-1)
      nsq(i) = nsq(i-1)
      do j = 1, Nthin
        call sample_cov(lglk(i),phi(i),nsq(i),phipars,nsqpars,phisc,nsqsc,&
           dm,F,betQ0,n,p,kappa,icf,acc)
        call sample_ssq(ssq(i))
        call samplez_ga(lglk(i),z(:,i),y,l,dft,ssq(i),tsq,n)
      end do
      call sample_beta(beta(:,i),z(:,i),ssq(i),n,p)
      call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),nsq(i),&
         n0,n,p,dmdm0,F,F0,kappa,icf)
      call rchkusr
    end do
  case (2) ! Binomial
    do j = 1, max(Nbi,Nthin)
      call sample_cov(lglk(i),phi(i),nsq(i),phipars,nsqpars,phisc,nsqsc,&
         dm,F,betQ0,n,p,kappa,icf,acc)
      call sample_ssq(ssq(i))
      call samplez_bi(lglk(i),z(:,i),y,l,dft,ssq(i),tsq,n)
    end do
    call sample_beta(beta(:,i),z(:,i),ssq(i),n,p)
    call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),nsq(i),&
       n0,n,p,dmdm0,F,F0,kappa,icf)
    call rchkusr
    do i = 2, Nout
      lglk(i) = lglk(i-1)
      z(:,i) = z(:,i-1)
      phi(i) = phi(i-1)
      nsq(i) = nsq(i-1)
      do j = 1, Nthin
        call sample_cov(lglk(i),phi(i),nsq(i),phipars,nsqpars,phisc,nsqsc,&
           dm,F,betQ0,n,p,kappa,icf,acc)
        call sample_ssq(ssq(i))
        call samplez_bi(lglk(i),z(:,i),y,l,dft,ssq(i),tsq,n)
      end do
      call sample_beta(beta(:,i),z(:,i),ssq(i),n,p)
      call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),nsq(i),&
         n0,n,p,dmdm0,F,F0,kappa,icf)
      call rchkusr
    end do
  case (3) ! Poisson
    do j = 1, max(Nbi,Nthin)
      call sample_cov(lglk(i),phi(i),nsq(i),phipars,nsqpars,phisc,nsqsc,&
         dm,F,betQ0,n,p,kappa,icf,acc)
      call sample_ssq(ssq(i))
      call samplez_po(lglk(i),z(:,i),y,l,dft,ssq(i),tsq,n)
    end do
    call sample_beta(beta(:,i),z(:,i),ssq(i),n,p)
    call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),nsq(i),&
       n0,n,p,dmdm0,F,F0,kappa,icf)
    call rchkusr
    do i = 2, Nout
      lglk(i) = lglk(i-1)
      z(:,i) = z(:,i-1)
      phi(i) = phi(i-1)
      nsq(i) = nsq(i-1)
      do j = 1, Nthin
        call sample_cov(lglk(i),phi(i),nsq(i),phipars,nsqpars,phisc,nsqsc,&
           dm,F,betQ0,n,p,kappa,icf,acc)
        call sample_ssq(ssq(i))
        call samplez_po(lglk(i),z(:,i),y,l,dft,ssq(i),tsq,n)
      end do
      call sample_beta(beta(:,i),z(:,i),ssq(i),n,p)
      call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),nsq(i),&
         n0,n,p,dmdm0,F,F0,kappa,icf)
      call rchkusr
    end do
  case (4) ! Gamma
    do j = 1, max(Nbi,Nthin)
      call sample_cov(lglk(i),phi(i),nsq(i),phipars,nsqpars,phisc,nsqsc,&
         dm,F,betQ0,n,p,kappa,icf,acc)
      call sample_ssq(ssq(i))
      call samplez_gm(lglk(i),z(:,i),y,l,dft,ssq(i),tsq,n)
    end do
    call sample_beta(beta(:,i),z(:,i),ssq(i),n,p)
    call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),nsq(i),&
       n0,n,p,dmdm0,F,F0,kappa,icf)
    call rchkusr
    do i = 2, Nout
      lglk(i) = lglk(i-1)
      z(:,i) = z(:,i-1)
      phi(i) = phi(i-1)
      nsq(i) = nsq(i-1)
      do j = 1, Nthin
        call sample_cov(lglk(i),phi(i),nsq(i),phipars,nsqpars,phisc,nsqsc,&
           dm,F,betQ0,n,p,kappa,icf,acc)
        call sample_ssq(ssq(i))
        call samplez_gm(lglk(i),z(:,i),y,l,dft,ssq(i),tsq,n)
      end do
      call sample_beta(beta(:,i),z(:,i),ssq(i),n,p)
      call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),nsq(i),&
         n0,n,p,dmdm0,F,F0,kappa,icf)
      call rchkusr
    end do
  end select
  call end_mcmc
end subroutine mcspsample


subroutine samplemulti (lglk, z, phi, nsq, y, l, F, &
   betm0, betQ0, ssqdf, ssqsc, kappa, icf, &
   nu, tsqdf, tsqin, dm, Ntot, Nout, Nbi, Nthin, n, p, kg, ifam)

  use mcmcfcns
  implicit none 
  integer, intent(in) :: n, p, ifam, kg, icf, &
     Nbi(kg), Nthin(kg), Nout(kg), Ntot
  double precision, intent(in) :: y(n), l(n), F(n,p), &
     nu(kg), tsqdf, tsqin, dm(n,n), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, phi(kg), nsq(kg), kappa(kg)
  double precision, intent(inout) :: z(n, Ntot)
  double precision, intent(out) :: lglk(Ntot)
  double precision :: ssq, tsq
  integer i, ii, j, k
  integer, parameter :: n0 = 0

  i = 0
  select case (ifam) ! Which family?
  case (0) ! Gaussian
    do k = 1, kg
      i = i + 1
      call ini_mcmc(lglk(i),z(:,i),phi(k),nsq(k),y,l,F,kappa(k),icf,&
         dm,&
         betm0,betQ0,ssqdf,ssqsc,tsqdf,tsqin,nu(k),n,n0,p,ifam)
      call rchkusr
      do j = 1, max(Nbi(k),Nthin(k))
        call sample_ssq(ssq)
        call sample_tsq(tsq)
        call samplez_gt(lglk(i),z(:,i),y,l,nu(k),ssq,tsq,n)
      end do
      call rchkusr
      do ii = 2, Nout(k)
        i = i + 1
        lglk(i) = lglk(i-1)
        z(:,i) = z(:,i-1)
        do j = 1, Nthin(k)
          call sample_ssq(ssq)
          call sample_tsq(tsq)
          call samplez_gt(lglk(i),z(:,i),y,l,nu(k),ssq,tsq,n)
        end do
        call rchkusr
      end do
      call end_mcmc
    end do
  case (1) ! Gaussian
    tsq = tsqin
    do k = 1, kg
      i = i + 1
      call ini_mcmc(lglk(i),z(:,i),phi(k),nsq(k),y,l,F,kappa(k),icf,&
         dm,&
         betm0,betQ0,ssqdf,ssqsc,tsqdf,tsq,nu(k),n,n0,p,ifam)
      call rchkusr
      do j = 1, max(Nbi(k),Nthin(k))
        call sample_ssq(ssq)
        call samplez_ga(lglk(i),z(:,i),y,l,nu(k),ssq,tsq,n)
      end do
      call rchkusr
      do ii = 2, Nout(k)
        i = i + 1
        lglk(i) = lglk(i-1)
        z(:,i) = z(:,i-1)
        do j = 1, Nthin(k)
          call sample_ssq(ssq)
          call samplez_ga(lglk(i),z(:,i),y,l,nu(k),ssq,tsq,n)
        end do
        call rchkusr
      end do
      call end_mcmc
    end do
  case (2) ! Binomial
    tsq = tsqin
    do k = 1, kg
      i = i + 1
      call ini_mcmc(lglk(i),z(:,i),phi(k),nsq(k),y,l,F,kappa(k),icf,&
         dm,&
         betm0,betQ0,ssqdf,ssqsc,tsqdf,tsq,nu(k),n,n0,p,ifam)
      call rchkusr
      do j = 1, max(Nbi(k),Nthin(k))
        call sample_ssq(ssq)
        call samplez_bi(lglk(i),z(:,i),y,l,nu(k),ssq,tsq,n)
      end do
      call rchkusr
      do ii = 2, Nout(k)
        i = i + 1
        lglk(i) = lglk(i-1)
        z(:,i) = z(:,i-1)
        do j = 1, Nthin(k)
          call sample_ssq(ssq)
          call samplez_bi(lglk(i),z(:,i),y,l,nu(k),ssq,tsq,n)
        end do
        call rchkusr
      end do
      call end_mcmc
    end do
  case (3) ! Poisson
    tsq = tsqin
    do k = 1, kg
      i = i + 1
      call ini_mcmc(lglk(i),z(:,i),phi(k),nsq(k),y,l,F,kappa(k),icf,&
         dm,&
         betm0,betQ0,ssqdf,ssqsc,tsqdf,tsq,nu(k),n,n0,p,ifam)
      call rchkusr
      do j = 1, max(Nbi(k),Nthin(k))
        call sample_ssq(ssq)
        call samplez_po(lglk(i),z(:,i),y,l,nu(k),ssq,tsq,n)
      end do
      call rchkusr
      do ii = 2, Nout(k)
        i = i + 1
        lglk(i) = lglk(i-1)
        z(:,i) = z(:,i-1)
        do j = 1, Nthin(k)
          call sample_ssq(ssq)
          call samplez_po(lglk(i),z(:,i),y,l,nu(k),ssq,tsq,n)
        end do
        call rchkusr
      end do
      call end_mcmc
    end do
  case (4) ! Gamma
    tsq = tsqin
    do k = 1, kg
      i = i + 1
      call ini_mcmc(lglk(i),z(:,i),phi(k),nsq(k),y,l,F,kappa(k),icf,&
         dm,&
         betm0,betQ0,ssqdf,ssqsc,tsqdf,tsq,nu(k),n,n0,p,ifam)
      call rchkusr
      do j = 1, max(Nbi(k),Nthin(k))
        call sample_ssq(ssq)
        call samplez_gm(lglk(i),z(:,i),y,l,nu(k),ssq,tsq,n)
      end do
      call rchkusr
      do ii = 2, Nout(k)
        i = i + 1
        lglk(i) = lglk(i-1)
        z(:,i) = z(:,i-1)
        do j = 1, Nthin(k)
          call sample_ssq(ssq)
          call samplez_gm(lglk(i),z(:,i),y,l,nu(k),ssq,tsq,n)
        end do
        call rchkusr
      end do
      call end_mcmc
    end do
  end select
end subroutine samplemulti


subroutine mcspsamtry (lglk, z, phi, nsq, acc, y, l, F, &
   betm0, betQ0, ssqdf, ssqsc, phipars, phisc, nsqpars, nsqsc, kappa, &
   icf, dft, tsq, dm, Nout, Npr, n, p, ifam)
  
  use mcmcfcns
  use msgmc
  implicit none 
  integer, intent(in) :: Nout, n, p, ifam, Npr, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), &
     dft, tsq, dm(n,n), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, phipars(4), nsqpars(4), kappa, phisc, nsqsc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), nsq(Nout)
  double precision, intent(out) :: lglk(Nout)
  integer, intent(out) :: acc
  integer i, ia
  double precision lglk1, ssq, phi1, nsq1, z1(n)
!   integer, parameter :: nl = 19
!   character(len=nl) :: label
  integer iap
  integer, parameter :: n0 = 0
  double precision :: tsqdf

  tsqdf = 0d0

!   write (label,'("Iteration",2x,"Acc% cov")')
!   call intpr (label,nl,0,0)
  call msgmca
!   call intpr ("-----------------------------------------------",nl,0,0)
  call msgmcl
  acc = 0
  ia = 0
  i = 0
  z1 = z(:,1)
  phi1 = phi(1)
  nsq1 = nsq(1)
  call ini_mcmc(lglk1,z1,phi1,nsq1,y,l,F,kappa,icf,dm,&
     betm0,betQ0,ssqdf,ssqsc,tsqdf,tsq,dft,n,n0,p,ifam)
  call rchkusr
  select case (ifam)
  case (1)
    do i = 1, Nout
      call sample_cov(lglk1,phi1,nsq1,phipars,nsqpars,phisc,nsqsc,dm,&
         F,betQ0,n,p,kappa,icf,ia)
      call sample_ssq(ssq)
      call samplez_ga(lglk1,z1,y,l,dft,ssq,tsq,n)
      lglk(i) = lglk1
      z(:,i) = z1
      phi(i) = phi1
      nsq(i) = nsq1
      if (Npr .gt. 0 .and. mod(i,Npr) .eq. 0) then
        iap = 100*ia/Npr
!         write (label,'(i9,2x,i8)') i, iap
!         call intpr (label,nl,0,0)
        call msgmci (i, iap)
        acc = acc + ia
        ia = 0
        call rchkusr
      end if
    end do
  case (2)
    do i = 1, Nout
      call sample_cov(lglk1,phi1,nsq1,phipars,nsqpars,phisc,nsqsc,dm,&
         F,betQ0,n,p,kappa,icf,ia)
      call sample_ssq(ssq)
      call samplez_bi(lglk1,z1,y,l,dft,ssq,tsq,n)
      lglk(i) = lglk1
      z(:,i) = z1
      phi(i) = phi1
      nsq(i) = nsq1
      if (Npr .gt. 0 .and. mod(i,Npr) .eq. 0) then
        iap = 100*ia/Npr
!         write (label,'(i9,2x,i8)') i, iap
!         call intpr (label,nl,0,0)
        call msgmci (i, iap)
        acc = acc + ia
        ia = 0
        call rchkusr
      end if
    end do
  case (3)
    do i = 1, Nout
      call sample_cov(lglk1,phi1,nsq1,phipars,nsqpars,phisc,nsqsc,dm,&
         F,betQ0,n,p,kappa,icf,ia)
      call sample_ssq(ssq)
      call samplez_po(lglk1,z1,y,l,dft,ssq,tsq,n)
      lglk(i) = lglk1
      z(:,i) = z1
      phi(i) = phi1
      nsq(i) = nsq1
      if (Npr .gt. 0 .and. mod(i,Npr) .eq. 0) then
        iap = 100*ia/Npr
!         write (label,'(i9,2x,i8)') i, iap
!         call intpr (label,nl,0,0)
        call msgmci (i, iap)
        acc = acc + ia
        ia = 0
        call rchkusr
      end if
    end do
  case (4)
    do i = 1, Nout
      call sample_cov(lglk1,phi1,nsq1,phipars,nsqpars,phisc,nsqsc,dm,&
         F,betQ0,n,p,kappa,icf,ia)
      call sample_ssq(ssq)
      call samplez_gm(lglk1,z1,y,l,dft,ssq,tsq,n)
      lglk(i) = lglk1
      z(:,i) = z1
      phi(i) = phi1
      nsq(i) = nsq1
      if (Npr .gt. 0 .and. mod(i,Npr) .eq. 0) then
        iap = 100*ia/Npr
!         write (label,'(i9,2x,i8)') i, iap
!         call intpr (label,nl,0,0)
        call msgmci (i, iap)
        acc = acc + ia
        ia = 0
        call rchkusr
      end if
    end do
  end select
  acc = acc + ia
  call end_mcmc
!   call intpr ("-----------------------------------------------",nl,0,0)
  call msgmcl
  iap = 100*acc/Nout
!   write (label,'("Avg acc %",2x,i8)') iap
!   call intpr (label,nl,0,0)
  call msgmce (iap)
!   call intpr ("-----------------------------------------------",nl,0,0)
  call msgmcl
end subroutine mcspsamtry


subroutine trgasample (lglk, z, z0, beta, ssq, tsq, phi, nsq, acc, y, l, F,&
   F0, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsqsc, phipars, phisc, nsqpars, &
   nsqsc, kappa, icf, dft, dm, dmdm0, Nout, Nbi, Nthin, n, n0, p)
  
  use mcmcfcns
  implicit none 
  integer, intent(in) :: Nout, Nbi, Nthin, n, n0, p, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), F0(n0,p), &
     dft, dm(n,n), dmdm0(n,n0), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, tsqdf, tsqsc, phipars(4), nsqpars(4), kappa, phisc, nsqsc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), nsq(Nout)
  double precision, intent(out) :: z0(n0,Nout), beta(p,Nout), ssq(Nout), &
     tsq(Nout), lglk(Nout)
  integer, intent(out) :: acc
  integer i, j
  integer, parameter :: ifam = 0

  acc = 0
  i = 1
  call ini_mcmc(lglk(i),z(:,i),phi(i),nsq(i),y,y,F,kappa,icf,dm,&
     betm0,betQ0,ssqdf,ssqsc,tsqdf,tsqsc,dft,n,n0,p,ifam)
  call rchkusr

  do j = 1, max(Nbi,Nthin)
    call sample_cov(lglk(i),phi(i),nsq(i),phipars,nsqpars,phisc,nsqsc,&
       dm,F,betQ0,n,p,kappa,icf,acc)
    call sample_ssq(ssq(i))
    call sample_tsq(tsq(i))
    call samplez_gt(lglk(i),z(:,i),y,l,dft,ssq(i),tsq(i),n)
  end do
  call sample_beta(beta(:,i),z(:,i),ssq(i),n,p)
  call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),nsq(i),&
       n0,n,p,dmdm0,F,F0,kappa,icf)
  call rchkusr
  do i = 2, Nout
    lglk(i) = lglk(i-1)
    z(:,i) = z(:,i-1)
    phi(i) = phi(i-1)
    nsq(i) = nsq(i-1)
    do j = 1, Nthin
      call sample_cov(lglk(i),phi(i),nsq(i),phipars,nsqpars,phisc,nsqsc,&
         dm,F,betQ0,n,p,kappa,icf,acc)
      call sample_ssq(ssq(i))
      call sample_tsq(tsq(i))
      call samplez_gt(lglk(i),z(:,i),y,l,dft,ssq(i),tsq(i),n)
    end do
    call sample_beta(beta(:,i),z(:,i),ssq(i),n,p)
    call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),nsq(i),&
       n0,n,p,dmdm0,F,F0,kappa,icf)
    call rchkusr
  end do

  call end_mcmc
end subroutine trgasample


subroutine trgasamtry (lglk, z, phi, nsq, acc, y, l, F, &
   betm0, betQ0, ssqdf, ssqsc, tsqdf, tsqsc, phipars, phisc, nsqpars, &
   nsqsc, kappa, icf, dft, dm, Nout, Npr, n, p)
  
  use mcmcfcns
  use msgmc
  implicit none 
  integer, intent(in) :: Nout, n, p, Npr, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), &
     dft, dm(n,n), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, tsqdf, tsqsc, phipars(4), nsqpars(4), kappa, phisc, nsqsc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), nsq(Nout)
  double precision, intent(out) :: lglk(Nout)
  integer, intent(out) :: acc
  integer i, ia
  double precision ssq, tsq, phi1, nsq1, z1(n), lglk1
!   integer, parameter :: nl = 19
!   character(len=nl) :: label
  integer iap
  integer, parameter :: n0 = 0, ifam = 0
  double precision :: dmdm0(n,n0), F0(n0,p)

  dmdm0 = 0d0
  F0 = 0d0
  
!   write (label,'("Iteration",2x,"Acc% cov")')
!   call intpr (label,nl,0,0)
  call msgmca
!   call intpr ("-----------------------------------------------",nl,0,0)
  call msgmcl
  acc = 0
  ia = 0
  i = 0
  z1 = z(:,1)
  phi1 = phi(1)
  nsq1 = nsq(1)
  call ini_mcmc(lglk1,z1,phi1,nsq1,y,y,F,kappa,icf,dm,&
     betm0,betQ0,ssqdf,ssqsc,tsqdf,tsqsc,dft,n,n0,p,ifam)
  call rchkusr

  do i = 1, Nout
    call sample_cov(lglk1,phi1,nsq1,phipars,nsqpars,phisc,nsqsc,dm,&
       F,betQ0,n,p,kappa,icf,ia)
    call sample_ssq(ssq)
    call sample_tsq(tsq)
    call samplez_gt(lglk1,z1,y,l,dft,ssq,tsq,n)
    lglk(i) = lglk1
    z(:,i) = z1
    phi(i) = phi1
    nsq(i) = nsq1
    if (Npr .gt. 0 .and. mod(i,Npr) .eq. 0) then
      iap = 100*ia/Npr
!       write (label,'(i9,2x,i8)') i, iap
!       call intpr (label,nl,0,0)
      call msgmci (i, iap)
      acc = acc + ia
      ia = 0
      call rchkusr
    end if
  end do
  acc = acc + ia
  call end_mcmc
!   call intpr ("-----------------------------------------------",nl,0,0)
  call msgmcl
  iap = 100*acc/Nout
!   write (label,'("Avg acc %",2x,i8)') iap
!   call intpr (label,nl,0,0)
  call msgmce (iap)
!   call intpr ("-----------------------------------------------",nl,0,0)
  call msgmcl
end subroutine trgasamtry

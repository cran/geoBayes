
!! Sample z0 many times when the covariance parameters change.
subroutine mcspz0mc (z0, Ns, z, beta, ssq, phi, nsq, F, F0, &
   betQ0, kappa, icf, dm, dmdm0, Nout, n, n0, p)
  
  use interfaces, only: rngini, rngend, randnorm
  use covfun, only: calc_cov, calc_cov_pred
  implicit none 
  integer, intent(in) :: Nout, n, n0, p, icf, Ns
  double precision, intent(in) :: F(n,p), F0(n0,p), &
     dm(n,n), dmdm0(n,n0), betQ0(p,p), kappa, ssq(Nout), &
     z(n,Nout), phi(Nout), nsq(Nout), beta(p, Nout)
  double precision, intent(out) :: z0(n0, Ns, Nout)
  integer i, j, k
  logical lup(n,n)
  double precision T(n,n), TiF(n,p), FTF(p,p), TFFT(n,p), Ups(n,n), &
     TC(n,n0), FCTF(n0,p), z0_ups(n0), ldh_Ups, z0_mean(n0), z0_sd(n0)

  do i = 1, n
    lup(:i-1,i) = .true.
    lup(i:,i) = .false. 
  end do

  call rngini
  do i = 1, Nout
    call rchkusr
    ! Compute covariance
    call calc_cov (phi(i),nsq(i),dm,F,betQ0, &
       lup,kappa,icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
    call calc_cov_pred(z0_ups, TC, FCTF, phi(i), nsq(i), dmdm0, F, &
       F0, kappa, icf, T, n, n0, p)
    call dgemv ('t',n,n0,1d0,TC,n,z,1,0d0,z0_mean,1) ! z0_mean = C'*T^{-1}*z
    call dgemv ('n',n0,p,1d0,FCTF,n0,beta(:,i),1,1d0,z0_mean,1) ! + (F0-C'T^-1F)*b
    z0_sd = sqrt(ssq(i))*z0_ups
    do j = 1, Ns
      do k = 1, n0
        z0(k,j,i) = z0_sd(k)*randnorm() + z0_mean(k)
      end do
    end do
  end do
  call rngend
end subroutine mcspz0mc

!! Sample z0 many times when the covariance parameters are the same.
subroutine mcspz0eb (z0, Ns, z, beta, ssq, phi, nsq, F, F0, &
   betQ0, kappa, icf, dm, dmdm0, Nout, n, n0, p)
  
  use interfaces, only: rngini, rngend, randnorm
  use covfun, only: calc_cov, calc_cov_pred
  implicit none 
  integer, intent(in) :: Nout, n, n0, p, icf, Ns
  double precision, intent(in) :: F(n,p), F0(n0,p), &
     dm(n,n), dmdm0(n,n0), betQ0(p,p), kappa, ssq(Nout), &
     z(n,Nout), phi, nsq, beta(p, Nout)
  double precision, intent(out) :: z0(n0, Ns, Nout)
  integer i, j, k
  logical lup(n,n)
  double precision T(n,n), TiF(n,p), FTF(p,p), TFFT(n,p), Ups(n,n), &
     TC(n,n0), FCTF(n0,p), z0_ups(n0), ldh_Ups, z0_mean(n0), z0_sd(n0), &
     z0_mfix(n0)

  do i = 1, n
    lup(:i-1,i) = .true.
    lup(i:,i) = .false. 
  end do

  call rngini
  ! Compute covariance
  call calc_cov (phi,nsq,dm,F,betQ0, &
     lup,kappa,icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
  call calc_cov_pred(z0_ups, TC, FCTF, phi, nsq, dmdm0, F, &
     F0, kappa, icf, T, n, n0, p)
  call dgemv ('t',n,n0,1d0,TC,n,z,1,0d0,z0_mfix,1) ! z0_mean = C'*T^{-1}*z
  do i = 1, Nout
    call rchkusr
    z0_mean = z0_mfix
    call dgemv ('n',n0,p,1d0,FCTF,n0,beta(:,i),1,1d0,z0_mean,1) ! + (F0-C'T^-1F)*b
    z0_sd = sqrt(ssq(i))*z0_ups
    do j = 1, Ns
      do k = 1, n0
        z0(k,j,i) = z0_sd(k)*randnorm() + z0_mean(k)
      end do
    end do
  end do
  call rngend
end subroutine mcspz0eb


subroutine DICmc (DIC, y, l, Ns, z, beta, ssq, phi, nsq, nu, F, F0, &
   kappa, icf, tsq, respdfh, dm, dm0, dmdm0, ifam, Nout, n, n0, p)

  ! tsq is tsqdfsc for transformed gaussian

  use interfaces, only: rngini, rngend, randnorm
  use covfun, only: calc_cov, calc_cov_pred, covmat
  use condyz

  implicit none 
  integer, intent(in) :: Nout, n, n0, p, icf, Ns, ifam
  double precision, intent(in) :: y(n0), l(n0), F(n,p), F0(n0,p), &
     dm(n,n), dm0(n0,n0), dmdm0(n,n0), kappa, ssq(Nout), &
     z(n,Nout), phi(Nout), nsq(Nout), beta(p,Nout), nu, tsq, respdfh
  double precision, intent(out) :: DIC
  integer i, j, k
  logical lup(n,n), lup0(n0,n0)
  double precision T(n,n), TC(n,n0), z0_mean(n0), u(n0), &
     z0(n0), Dy, T0(n0,n0), nsqp1, cnt, ssqrt, zz(n), alphz, alphz0
  
  do i = 1, n
    lup(:i-1,i) = .true.
    lup(i:,i) = .false. 
  end do

  do i = 1, n0
    lup0(:i-1,i) = .true.
    lup0(i:,i) = .false. 
  end do

  DIC = 0d0
  cnt = 0d0

  if (all (F .eq. 0d0)) then
    alphz = 0d0
  else
    alphz = -1d0
  end if
  if (all (F0 .eq. 0d0)) then
    alphz0 = 0d0
  else
    alphz0 = 1d0
  end if
  
  call rngini
  
  do i = 1, Nout
    call rchkusr
    ! Compute covariance
    nsqp1 = nsq(i) + 1d0
    T = dm
    call covmat (T, phi(i), kappa, n, n, lup, icf)
    do k = 1, n
      T(k, k) = nsqp1
    end do
    call dpotrf ('u', n, T, n, j)
    T0 = dm0
    call covmat (T0, phi(i), kappa, n0, n0, lup0, icf)
    do k = 1, n0
      T0(k, k) = nsqp1
    end do
    TC = dmdm0
    call covmat (TC, phi(i), kappa, n, n0, icf)
    call dtrsm ('l', 'u', 't', 'n', n, n0, 1d0, T, n, TC, n)
    call dsyrk ('u', 't', n0, n, -1d0, TC, n, 1d0, T0, n0)
    call dpotrf ('u', n0, T0, n0, j)

    ssqrt = sqrt(ssq(i))
    zz = z(:,i)
    call dgemv ('n', n, p, alphz, F, n, beta(:,i), 1, 1d0, zz, 1)
    call dtrmv ('u', 't', 'n', n, T, n, zz, 1)
    call dgemv ('n', n0, p, alphz0, F0, n0, beta(:,i), 1, 0d0, z0_mean, 1)
    call dgemv ('t', n, n0, 1d0, TC, n, zz, 1, 1d0, z0_mean, 1)

    do j = 1, Ns
      cnt = cnt + 1d0
      do k = 1, n0
        u(k) = randnorm()
        z0(k) = u(k)*ssqrt
      end do
      call dtrmv ('u', 't', 'n', n0, T0, n0, z0, 1)
      z0 = z0 + z0_mean
      select case (ifam)
      case (0) ! Transformed Gaussian
        Dy = -2d0*condyz_gt(n0, y, l, z0, nu, tsq, respdfh)
      case (1) ! Gaussian
        Dy = -2d0*condyz_ga(n0, y, l, z0, nu, tsq)
      case (2) ! Binomial
        Dy = -2d0*condyz_bi(n0, y, l, z0, nu, tsq)
      case (3) ! Poisson
        Dy = -2d0*condyz_po(n0, y, l, z0, nu, tsq)
      case (4) ! Gamma
        Dy = -2d0*condyz_gm(n0, y, l, z0, nu, tsq)
      case (5) ! Binomial Asymmetric
        Dy = -2d0*condyz_ba(n0, y, l, z0, nu, tsq)
      case (6) ! Binomial Asymmetric Decreasing
        Dy = -2d0*condyz_bd(n0, y, l, z0, nu, tsq)
      end select
      DIC = DIC + (Dy - DIC)/cnt
    end do
  end do
  
  call rngend
end subroutine DICmc


subroutine DICeb (DIC, y, l, Ns, z, beta, ssq, phi, nsq, nu, F, F0, &
   kappa, icf, tsq, respdfh, dm, dm0, dmdm0, ifam, Nout, n, n0, p)

  ! tsq is tsqdfsc for transformed gaussian

  use interfaces, only: rngini, rngend, randnorm
  use covfun, only: calc_cov, calc_cov_pred, covmat
  use condyz

  implicit none 
  integer, intent(in) :: Nout, n, n0, p, icf, Ns, ifam
  double precision, intent(in) :: y(n0), l(n0), F(n,p), F0(n0,p), &
     dm(n,n), dm0(n0,n0), dmdm0(n,n0), kappa, ssq(Nout), &
     z(n,Nout), phi, nsq, beta(p,Nout), nu, tsq, respdfh
  double precision, intent(out) :: DIC
  integer i, j, k
  logical lup(n,n), lup0(n0,n0)
  double precision T(n,n), TC(n,n0), z0_mean(n0), u(n0), &
     z0(n0), Dy, T0(n0,n0), nsqp1, cnt, zz(n), ssqrt, alphz, alphz0
  
  do i = 1, n
    lup(:i-1,i) = .true.
    lup(i:,i) = .false. 
  end do

  do i = 1, n0
    lup0(:i-1,i) = .true.
    lup0(i:,i) = .false. 
  end do

  DIC = 0d0
  cnt = 0d0

  if (all (F .eq. 0d0)) then
    alphz = 0d0
  else
    alphz = -1d0
  end if
  if (all (F0 .eq. 0d0)) then
    alphz0 = 0d0
  else
    alphz0 = 1d0
  end if

  call rngini
  
  ! Compute covariance
  nsqp1 = nsq + 1d0
  T = dm
  call covmat (T, phi, kappa, n, n, lup, icf)
  do k = 1, n
    T(k, k) = nsqp1
  end do
  call dpotrf ('u', n, T, n, j)
  T0 = dm0
  call covmat (T0, phi, kappa, n0, n0, lup0, icf)
  do k = 1, n0
    T0(k, k) = nsqp1
  end do
  TC = dmdm0
  call covmat (TC, phi, kappa, n, n0, icf)
  call dtrsm ('l', 'u', 't', 'n', n, n0, 1d0, T, n, TC, n)
  call dsyrk ('u', 't', n0, n, -1d0, TC, n, 1d0, T0, n0)
  call dpotrf ('u', n0, T0, n0, j)

  do i = 1, Nout
    call rchkusr
    ssqrt = sqrt(ssq(i))
    zz = z(:,i)
    call dgemv ('n', n, p, alphz, F, n, beta(:,i), 1, 1d0, zz, 1)
    call dtrmv ('u', 't', 'n', n, T, n, zz, 1)
    call dgemv ('n', n0, p, alphz0, F0, n0, beta(:,i), 1, 0d0, z0_mean, 1)
    call dgemv ('t', n, n0, 1d0, TC, n, zz, 1, 1d0, z0_mean, 1)
    do j = 1, Ns
      cnt = cnt + 1d0
      do k = 1, n0
        u(k) = randnorm()
        z0(k) = u(k)*ssqrt
      end do
      call dtrmv ('u', 't', 'n', n0, T0, n0, z0, 1)
      z0 = z0 + z0_mean
      select case (ifam)
      case (0) ! Transformed Gaussian
        Dy = -2d0*condyz_gt(n0, y, l, z0, nu, tsq, respdfh)
      case (1) ! Gaussian
        Dy = -2d0*condyz_ga(n0, y, l, z0, nu, tsq)
      case (2) ! Binomial
        Dy = -2d0*condyz_bi(n0, y, l, z0, nu, tsq)
      case (3) ! Poisson
        Dy = -2d0*condyz_po(n0, y, l, z0, nu, tsq)
      case (4) ! Gamma
        Dy = -2d0*condyz_gm(n0, y, l, z0, nu, tsq)
      case (5) ! Binomial Asymmetric
        Dy = -2d0*condyz_ba(n0, y, l, z0, nu, tsq)
      case (6) ! Binomial Asymmetric Decreasing
        Dy = -2d0*condyz_bd(n0, y, l, z0, nu, tsq)
      end select
      DIC = DIC + (Dy - DIC)/cnt
    end do
  end do
  
  call rngend
end subroutine DICeb

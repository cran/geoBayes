!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calcb_st (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, zsample, weights, n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam)

  use interfaces
  use flogsumexp
  use covfun
  use linkfcn
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), zsample(n, Ntot), weights(Ntot), &
     betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n, p), dm(n, n)
  double precision, intent(out) :: bfact(n_nu, n_cov)
  logical lup(n, n)
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh, zmxi(n, Ntot)
  integer i, j, k
  double precision logfy(n_nu, Ntot), lfy
  double precision T(n, n), TiF(n, p), FTF(p, p), TFFT(n, p), Ups(n, n), &
     ldh_Ups, llikw(n_nu, Ntot), mu, xi(n), Upsz(n), zUz

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
    zmxi = zsample
  else ! Normal prior
    modeldfh = .5d0*(n + ssqdf)
    xi = matmul(F,betm0)
    zmxi = spread(xi, 2, Ntot)
    zmxi = zsample - zmxi
  end if

  lup = .false.
  do i = 2, n
    lup(:i-1,i) = .true.
  end do

  call rchkusr

  ! Calculate log f(y|z,nu)
  select case (ifam)
  case (0)
    do j = 1, Ntot
      do k = 1, n_nu
        lfy = tsqdfsc
        do i = 1, n
          mu = invlink_ga(zsample(i,j),nu(k))
          mu = y(i) - mu
          lfy = lfy + l(i)*mu*mu
        end do
        logfy(k, j) = -respdfh*log(lfy)
      end do
      call rchkusr
    end do
  case (1)
    do j = 1, Ntot
      do k = 1, n_nu
        lfy = 0d0
        do i = 1, n
          mu = invlink_ga(zsample(i,j),nu(k))
          lfy = lfy + y(i)*mu - .5d0*l(i)*mu*mu
        end do
        logfy(k, j) = lfy/tsq
      end do
      call rchkusr
    end do
  case (2)
    do j = 1, Ntot
      do k = 1, n_nu
        lfy = 0d0
        do i = 1, n
          mu = invlink_bi(zsample(i,j),nu(k))
          lfy = lfy + y(i)*mu + l(i)*flog1mexp(mu)
        end do
        logfy(k, j) = lfy/tsq
      end do
      call rchkusr
    end do
  case (3)
    do j = 1, Ntot
      do k = 1, n_nu
        lfy = 0d0
        do i = 1, n
          mu = invlink_po(zsample(i,j),nu(k))
          lfy = lfy + y(i)*mu - l(i)*exp(mu)
        end do
        logfy(k, j) = lfy/tsq
      end do
      call rchkusr
    end do
  case (4)
    do j = 1, Ntot
      do k = 1, n_nu
        lfy = 0d0
        do i = 1, n
          mu = invlink_gm(zsample(i,j),nu(k))
          lfy = lfy + y(i)*mu - l(i)*log(-mu)
        end do
        logfy(k, j) = lfy/tsq
      end do
      call rchkusr
    end do
  end select

  do k = 1, n_cov
    call calc_cov (phi(k),nsq(k),dm,F,betQ0,&
       lup,kappa(k),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
    do j = 1, Ntot
      call dsymv ('u',n,1d0,Ups,n,zmxi(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
      zUz = dot_product(zmxi(:,j),Upsz) + ssqdfsc
      lfy = ldh_Ups - modeldfh*log(zUz)
      do i = 1, n_nu
        llikw(i, j) = logfy(i, j) + lfy - weights(j)
      end do
    end do
    bfact(:, k) = logrsumexp(llikw, n_nu, Ntot)
  end do
end subroutine calcb_st


subroutine calcb_cv (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, zsample, weights, zcv, n, p, kg, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam)

  use interfaces
  use flogsumexp
  use covfun
  use linkfcn
  implicit none
  integer, intent(in) :: n, p, kg, Ntot, n_cov, n_nu, ifam, icf
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), zsample(n, Ntot), weights(Ntot), zcv(Ntot, kg), &
     betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n, p), dm(n, n)
  double precision, intent(out) :: bfact(n_nu, n_cov)
  logical lup(n, n)
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh, zmxi(n, Ntot)
  integer i, j, k
  double precision, allocatable :: tmp(:)
  integer ltmp
  double precision qrtau(kg), zcvqr(Ntot, kg), logfy(n_nu, Ntot), lfy, dNtot
  double precision T(n, n), TiF(n, p), FTF(p, p), TFFT(n, p), Ups(n, n), &
     ldh_Ups, ycv(n_nu, Ntot), llikw, mu, xi(n), Upsz(n), zUz

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)
  dNtot = log(dble(Ntot))

  ! Determine flat or normal prior
  j=0
  do i = 1, p
    if (betQ0(i,i) /= 0d0) j = j + 1
  end do
  if (j == 0) then ! Flat prior
    modeldfh = .5d0*(n - p + ssqdf)
    xi = 0d0
    zmxi = zsample
  else ! Normal prior
    modeldfh = .5d0*(n + ssqdf)
    xi = matmul(F,betm0)
    zmxi = spread(xi, 2, Ntot)
    zmxi = zsample - zmxi
  end if

  lup = .false.
  do i = 2, n
    lup(:i-1,i) = .true.
  end do

  ! QR factorisation
  zcvqr = zcv
  ! Determine optimal work size for QR factorisation.
  ltmp = -1
  call dgeqrf (Ntot, kg, zcvqr, Ntot, qrtau, lfy, ltmp, i)
  if (i .ne. 0) then
    call rexit ('calcb - Cannot determine optimal space for regression')
  else
    ltmp = nint(lfy)
    allocate (tmp(ltmp))
  end if
  call dgeqrf (Ntot, kg, zcvqr, Ntot, qrtau, tmp, ltmp, i)
  if (i .ne. 0) then
    call rexit ('calcb - Error with the QR factorisation')
  end if

  call rchkusr

  ! Calculate log f(y|z,nu)
  select case (ifam)
  case (0)
    do j = 1, Ntot
      do k = 1, n_nu
        lfy = tsqdfsc
        do i = 1, n
          mu = invlink_ga(zsample(i,j),nu(k))
          mu = y(i) - mu
          lfy = lfy + l(i)*mu*mu
        end do
        logfy(k, j) = -respdfh*log(lfy)
      end do
      call rchkusr
    end do
  case (1)
    do j = 1, Ntot
      do k = 1, n_nu
        lfy = 0d0
        do i = 1, n
          mu = invlink_ga(zsample(i,j),nu(k))
          lfy = lfy + y(i)*mu - .5d0*l(i)*mu*mu
        end do
        logfy(k, j) = lfy/tsq
      end do
      call rchkusr
    end do
  case (2)
    do j = 1, Ntot
      do k = 1, n_nu
        lfy = 0d0
        do i = 1, n
          mu = invlink_bi(zsample(i,j),nu(k))
          lfy = lfy + y(i)*mu + l(i)*flog1mexp(mu)
        end do
        logfy(k, j) = lfy/tsq
      end do
      call rchkusr
    end do
  case (3)
    do j = 1, Ntot
      do k = 1, n_nu
        lfy = 0d0
        do i = 1, n
          mu = invlink_po(zsample(i,j),nu(k))
          lfy = lfy + y(i)*mu - l(i)*exp(mu)
        end do
        logfy(k, j) = lfy/tsq
      end do
      call rchkusr
    end do
  case (4)
    do j = 1, Ntot
      do k = 1, n_nu
        lfy = 0d0
        do i = 1, n
          mu = invlink_gm(zsample(i,j),nu(k))
          lfy = lfy + y(i)*mu - l(i)*log(-mu)
        end do
        logfy(k, j) = lfy/tsq
      end do
      call rchkusr
    end do
  end select

  do k = 1, n_cov
    call calc_cov (phi(k),nsq(k),dm,F,betQ0,&
       lup,kappa(k),icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
    do j = 1, Ntot
      call dsymv ('u',n,1d0,Ups,n,zmxi(:,j),1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
      zUz = dot_product(zmxi(:,j),Upsz) + ssqdfsc
      
      ! Calculate unnormalised log-likelihood at sampled points
      lfy = ldh_Ups - modeldfh*log(zUz)
      do i = 1, n_nu
        llikw = logfy(i, j) + lfy - weights(j)
        ycv(i, j) = exp(llikw + dNtot)
      end do
    end do

    ! Regression
    call dormqr ('r','n',n_nu,Ntot,kg,zcvqr,Ntot,qrtau,ycv,n_nu,tmp,ltmp,i)
    if (i .ne. 0) then
      call rexit ('calcb - Error in DORMQR')
    end if
    call dtrsm ('r','u','t','n',n_nu,kg,1d0,zcvqr,Ntot,ycv,n_nu)
    where (ycv(:, 1) > 0d0)
      bfact(:,k) = log(ycv(:, 1))
    elsewhere
      bfact(:,k) = -huge(1d0)
    end where
  end do
  deallocate(tmp)
end subroutine calcb_cv


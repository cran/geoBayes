!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calcbz_st (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, zsample, weights, n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam)

  use interfaces
  use flogsumexp
  use covfun
  use condyz
  use pdfz
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), zsample(n, Ntot), weights(Ntot), &
     betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n, p), dm(n, n)
  double precision, intent(out) :: bfact(n_nu, n_cov)
  logical lup(n, n), lmxi
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh
  integer i, j, k
  double precision logfy(n_nu, Ntot), lfz
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, llikw(n_nu, Ntot), xi(n)

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
  case (1)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_ga(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case (2)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_bi(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case (3)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_po(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case (4)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_gm(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case (5)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_ba(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case (6)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_bd(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case default
    call rexit ("Unrecognised family")
  end select

  do k = 1, n_cov
    call calc_cov (phi(k),nsq(k),dm,F,betQ0,&
       lup,kappa(k),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
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
end subroutine calcbz_st


subroutine calcbz_cv (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, zsample, weights, zcv, n, p, kg, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam)

  use interfaces
  use flogsumexp
  use covfun
  use condyz
  use pdfz
  implicit none
  integer, intent(in) :: n, p, kg, Ntot, n_cov, n_nu, ifam, icf
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), zsample(n, Ntot), weights(Ntot), zcv(Ntot, kg), &
     betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n, p), dm(n, n)
  double precision, intent(out) :: bfact(n_nu, n_cov)
  logical lup(n, n), lmxi
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh
  integer i, j, k
  double precision, allocatable :: tmp(:)
  integer ltmp
  double precision qrtau(kg), zcvqr(Ntot, kg), logfy(n_nu, Ntot), lfz, dNtot
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ycv(n_nu, Ntot), llikw, xi(n)

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
  case (1)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_ga(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case (2)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_bi(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case (3)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_po(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case (4)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_gm(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case (5)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_ba(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case (6)
    do j = 1, Ntot
      do k = 1, n_nu
        logfy(k, j) = condyz_bd(n, y, l, zsample(:, j), nu(k), tsq)
      end do
      call rchkusr
    end do
  case default
    call rexit ("Unrecognised family")
  end select

  ! QR factorisation
  zcvqr = zcv
  ! Determine optimal work size for QR factorisation.
  ltmp = -1
  call dgeqrf (Ntot, kg, zcvqr, Ntot, qrtau, lfz, ltmp, i)
  if (i .ne. 0) then
    call rexit ('calcb - Cannot determine optimal space for regression')
  else
    ltmp = nint(lfz)
    allocate (tmp(ltmp))
  end if
  call dgeqrf (Ntot, kg, zcvqr, Ntot, qrtau, tmp, ltmp, i)
  if (i .ne. 0) then
    call rexit ('calcb - Error with the QR factorisation')
  end if

  do k = 1, n_cov
    call calc_cov (phi(k),nsq(k),dm,F,betQ0,&
       lup,kappa(k),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
    do j = 1, Ntot
      ! Calculate unnormalised log-likelihood at sampled points
      lfz = logpdfz(n, zsample(:, j), Ups, ldh_Ups, xi, lmxi, &
         ssqdfsc, modeldfh)
      do i = 1, n_nu
        llikw = logfy(i, j) + lfz - weights(j)
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
end subroutine calcbz_cv



subroutine calcbmu_st (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, musample, weights, n, p, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam)

  use interfaces
  use flogsumexp
  use covfun
  use condymu
  use pdfmu
  implicit none
  integer, intent(in) :: n, p, Ntot, n_cov, n_nu, ifam, icf
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), musample(n, Ntot), weights(Ntot), &
     betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n, p), dm(n, n)
  double precision, intent(out) :: bfact(n_nu, n_cov)
  logical lup(n, n), lmxi
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh
  integer i, j, k
  double precision logfy(Ntot), lfmu
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, llikw(n_nu, Ntot), xi(n)

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

  call rchkusr

  select case (ifam)
  case (0)
    do j = 1, Ntot
      logfy(j) = condymu_gt(n, y, l, musample(:, j), tsqdfsc, respdfh)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_ga(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw(k, j) = logfy(j) + lfmu - weights(j)
        end do
      end do
      bfact(:, i) = logrsumexp(llikw, n_nu, Ntot)
    end do
  case (1)
    do j = 1, Ntot
      logfy(j) = condymu_ga(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_ga(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw(k, j) = logfy(j) + lfmu - weights(j)
        end do
      end do
      bfact(:, i) = logrsumexp(llikw, n_nu, Ntot)
    end do
  case (2)
    do j = 1, Ntot
      logfy(j) = condymu_bi(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_bi(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw(k, j) = logfy(j) + lfmu - weights(j)
        end do
      end do
      bfact(:, i) = logrsumexp(llikw, n_nu, Ntot)
    end do
  case (3)
    do j = 1, Ntot
      logfy(j) = condymu_po(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_po(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw(k, j) = logfy(j) + lfmu - weights(j)
        end do
      end do
      bfact(:, i) = logrsumexp(llikw, n_nu, Ntot)
    end do
  case (4)
    do j = 1, Ntot
      logfy(j) = condymu_gm(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_gm(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw(k, j) = logfy(j) + lfmu - weights(j)
        end do
      end do
      bfact(:, i) = logrsumexp(llikw, n_nu, Ntot)
    end do
  case (5)
    do j = 1, Ntot
      logfy(j) = condymu_ba(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_ba(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw(k, j) = logfy(j) + lfmu - weights(j)
        end do
      end do
      bfact(:, i) = logrsumexp(llikw, n_nu, Ntot)
    end do
  case (6)
    do j = 1, Ntot
      logfy(j) = condymu_bd(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_bd(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw(k, j) = logfy(j) + lfmu - weights(j)
        end do
      end do
      bfact(:, i) = logrsumexp(llikw, n_nu, Ntot)
    end do
  case default
    call rexit ("Unrecognised family")
  end select
end subroutine calcbmu_st


subroutine calcbmu_cv (bfact, phi, nu, nsq, kappa, icf, n_cov, n_nu, &
   Ntot, musample, weights, zcv, n, p, kg, betm0, betQ0, ssqdf, &
   ssqsc, tsqdf, tsq, y, l, F, dm, ifam)

  use interfaces
  use flogsumexp
  use covfun
  use condymu
  use pdfmu
  implicit none
  integer, intent(in) :: n, p, kg, Ntot, n_cov, n_nu, ifam, icf
  double precision, intent(in) :: phi(n_cov), nsq(n_cov), kappa(n_cov), &
     nu(n_nu), musample(n, Ntot), weights(Ntot), zcv(Ntot, kg), &
     betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n, p), dm(n, n)
  double precision, intent(out) :: bfact(n_nu, n_cov)
  logical lup(n, n), lmxi
  double precision ssqdfsc, tsqdfsc, respdfh, modeldfh
  integer i, j, k
  double precision, allocatable :: tmp(:)
  integer ltmp
  double precision qrtau(kg), zcvqr(Ntot, kg), logfy(Ntot), lfmu, dNtot
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ycv(n_nu, Ntot), llikw, xi(n)

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

  call rchkusr

  ! QR factorisation
  zcvqr = zcv
  ! Determine optimal work size for QR factorisation.
  ltmp = -1
  call dgeqrf (Ntot, kg, zcvqr, Ntot, qrtau, lfmu, ltmp, i)
  if (i .ne. 0) then
    call rexit ('calcb - Cannot determine optimal space for regression')
  else
    ltmp = nint(lfmu)
    allocate (tmp(ltmp))
  end if
  call dgeqrf (Ntot, kg, zcvqr, Ntot, qrtau, tmp, ltmp, i)
  if (i .ne. 0) then
    call rexit ('calcb - Error with the QR factorisation')
  end if

  select case (ifam)
  case (0)
    do j = 1, Ntot
      logfy(j) = condymu_gt(n, y, l, musample(:, j), tsqdfsc, respdfh)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_ga(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw = logfy(j) + lfmu - weights(j)
          ycv(k, j) = exp(llikw + dNtot)
        end do
      end do
      ! Regression
      call dormqr ('r','n',n_nu,Ntot,kg,zcvqr,Ntot,qrtau,ycv,n_nu,tmp,ltmp,j)
      if (j .ne. 0) then
        call rexit ('calcb - Error in DORMQR')
      end if
      call dtrsm ('r','u','t','n',n_nu,kg,1d0,zcvqr,Ntot,ycv,n_nu)
      where (ycv(:, 1) > 0d0)
        bfact(:,i) = log(ycv(:, 1))
      elsewhere
        bfact(:,i) = -huge(1d0)
      end where
    end do
  case (1)
    do j = 1, Ntot
      logfy(j) = condymu_ga(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_ga(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw = logfy(j) + lfmu - weights(j)
          ycv(k, j) = exp(llikw + dNtot)
        end do
      end do
      ! Regression
      call dormqr ('r','n',n_nu,Ntot,kg,zcvqr,Ntot,qrtau,ycv,n_nu,tmp,ltmp,j)
      if (j .ne. 0) then
        call rexit ('calcb - Error in DORMQR')
      end if
      call dtrsm ('r','u','t','n',n_nu,kg,1d0,zcvqr,Ntot,ycv,n_nu)
      where (ycv(:, 1) > 0d0)
        bfact(:,i) = log(ycv(:, 1))
      elsewhere
        bfact(:,i) = -huge(1d0)
      end where
    end do
  case (2)
    do j = 1, Ntot
      logfy(j) = condymu_bi(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_bi(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw = logfy(j) + lfmu - weights(j)
          ycv(k, j) = exp(llikw + dNtot)
        end do
      end do
      ! Regression
      call dormqr ('r','n',n_nu,Ntot,kg,zcvqr,Ntot,qrtau,ycv,n_nu,tmp,ltmp,j)
      if (j .ne. 0) then
        call rexit ('calcb - Error in DORMQR')
      end if
      call dtrsm ('r','u','t','n',n_nu,kg,1d0,zcvqr,Ntot,ycv,n_nu)
      where (ycv(:, 1) > 0d0)
        bfact(:,i) = log(ycv(:, 1))
      elsewhere
        bfact(:,i) = -huge(1d0)
      end where
    end do
  case (3)
    do j = 1, Ntot
      logfy(j) = condymu_po(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_po(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw = logfy(j) + lfmu - weights(j)
          ycv(k, j) = exp(llikw + dNtot)
        end do
      end do
      ! Regression
      call dormqr ('r','n',n_nu,Ntot,kg,zcvqr,Ntot,qrtau,ycv,n_nu,tmp,ltmp,j)
      if (j .ne. 0) then
        call rexit ('calcb - Error in DORMQR')
      end if
      call dtrsm ('r','u','t','n',n_nu,kg,1d0,zcvqr,Ntot,ycv,n_nu)
      where (ycv(:, 1) > 0d0)
        bfact(:,i) = log(ycv(:, 1))
      elsewhere
        bfact(:,i) = -huge(1d0)
      end where
    end do
  case (4)
    do j = 1, Ntot
      logfy(j) = condymu_gm(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_gm(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw = logfy(j) + lfmu - weights(j)
          ycv(k, j) = exp(llikw + dNtot)
        end do
      end do
      ! Regression
      call dormqr ('r','n',n_nu,Ntot,kg,zcvqr,Ntot,qrtau,ycv,n_nu,tmp,ltmp,j)
      if (j .ne. 0) then
        call rexit ('calcb - Error in DORMQR')
      end if
      call dtrsm ('r','u','t','n',n_nu,kg,1d0,zcvqr,Ntot,ycv,n_nu)
      where (ycv(:, 1) > 0d0)
        bfact(:,i) = log(ycv(:, 1))
      elsewhere
        bfact(:,i) = -huge(1d0)
      end where
    end do
  case (5)
    do j = 1, Ntot
      logfy(j) = condymu_ba(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_ba(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw = logfy(j) + lfmu - weights(j)
          ycv(k, j) = exp(llikw + dNtot)
        end do
      end do
      ! Regression
      call dormqr ('r','n',n_nu,Ntot,kg,zcvqr,Ntot,qrtau,ycv,n_nu,tmp,ltmp,j)
      if (j .ne. 0) then
        call rexit ('calcb - Error in DORMQR')
      end if
      call dtrsm ('r','u','t','n',n_nu,kg,1d0,zcvqr,Ntot,ycv,n_nu)
      where (ycv(:, 1) > 0d0)
        bfact(:,i) = log(ycv(:, 1))
      elsewhere
        bfact(:,i) = -huge(1d0)
      end where
    end do
  case (6)
    do j = 1, Ntot
      logfy(j) = condymu_bd(n, y, l, musample(:, j), tsq)
    end do
    do i = 1, n_cov
      call rchkusr
      call calc_cov (phi(i),nsq(i),dm,F,betQ0,&
         lup,kappa(i),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        do k = 1, n_nu
          lfmu = logpdfmu_bd(n, musample(:, j), Ups, ldh_Ups, &
             nu(k), xi, lmxi, ssqdfsc, modeldfh)
          llikw = logfy(j) + lfmu - weights(j)
          ycv(k, j) = exp(llikw + dNtot)
        end do
      end do
      ! Regression
      call dormqr ('r','n',n_nu,Ntot,kg,zcvqr,Ntot,qrtau,ycv,n_nu,tmp,ltmp,j)
      if (j .ne. 0) then
        call rexit ('calcb - Error in DORMQR')
      end if
      call dtrsm ('r','u','t','n',n_nu,kg,1d0,zcvqr,Ntot,ycv,n_nu)
      where (ycv(:, 1) > 0d0)
        bfact(:,i) = log(ycv(:, 1))
      elsewhere
        bfact(:,i) = -huge(1d0)
      end where
    end do
  case default
    call rexit ("Unrecognised family")
  end select

  deallocate(tmp)
end subroutine calcbmu_cv


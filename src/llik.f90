subroutine llikfcnz (lglk, philist, nsqlist, nulist, kappalist, &
   zsample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, kg, ifam)
  
  use covfun
  use jointyz
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, Ntot, icf
  double precision, intent(in) :: philist(kg), nsqlist(kg), nulist(kg), &
     zsample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot, kg)
  logical lup(n, n), lmxi
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, &
     tsqdfsc, respdfh, xi(n)
  integer i, ii, j

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  do i = 1, n
    lup(:i-1,i) = .true.
    lup(i:,i) = .false. 
  end do

  select case (ifam)
  case (0) ! Transformed Gaussian
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointyz_gt(n, zsample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
      end do
    end do
  case (1) ! Gaussian
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointyz_ga(n, zsample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (2) ! Binomial
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointyz_bi(n, zsample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (3) ! Poisson
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointyz_po(n, zsample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (4) ! Gamma
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointyz_gm(n, zsample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (5) ! Binomial Asymmetric
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointyz_ba(n, zsample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (6) ! Binomial Asymmetric Decreasing
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointyz_bd(n, zsample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (7) ! Binomial Wallace
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointyz_bw(n, zsample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case default
    call rexit ("Unrecognised family")
  end select
end subroutine llikfcnz




subroutine llikfcnmu (lglk, philist, nsqlist, nulist, kappalist, &
   musample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, kg, ifam)
  
  use covfun
  use jointymu
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, Ntot, icf
  double precision, intent(in) :: philist(kg), nsqlist(kg), nulist(kg), &
     musample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot, kg)
  logical lup(n, n), lmxi
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, &
     tsqdfsc, respdfh, xi(n)
  integer i, ii, j

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  do i = 1, n
    lup(:i-1,i) = .true.
    lup(i:,i) = .false. 
  end do

  select case (ifam)
  case (0) ! Transformed Gaussian
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointymu_gt(n, musample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
      end do
    end do
  case (1) ! Gaussian
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointymu_ga(n, musample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (2) ! Binomial
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointymu_bi(n, musample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (3) ! Poisson
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointymu_po(n, musample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (4) ! Gamma
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointymu_gm(n, musample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (5) ! Binomial Asymmetric
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointymu_ba(n, musample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (6) ! Binomial Asymmetric Decreasing
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointymu_bd(n, musample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case (7) ! Binomial Wallace
    do ii = 1, kg
      call calc_cov (philist(ii),nsqlist(ii),dm,F,betQ0,&
         lup,kappalist(ii),icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,ii) = jointymu_bw(n, musample(:, j), y, l, Ups, ldh_Ups, &
           nulist(ii), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  case default
    call rexit ("Unrecognised family")
  end select
end subroutine llikfcnmu


subroutine llikfcnrb (lglk, philist, nsqlist, nulist, kappalist, &
   wsample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, kg, ifam)
  
  use covfun
  use transfbinomial, only: jointyw_bi
  use jointyz, only: jointyz_bi
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, Ntot, icf
  double precision, intent(in) :: philist(kg), nsqlist(kg), nulist(kg), &
     wsample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot, kg)
  logical lup(n, n), lmxi
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, nu, phi, nsq, kappa, &
     tsqdfsc, respdfh, xi(n)
  integer i, ii, j

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  do i = 1, n
    lup(:i-1,i) = .true.
    lup(i:,i) = .false. 
  end do

  select case (ifam)
  case (2) ! Binomial
    do ii = 1, kg
      nu = nulist(ii)
      phi = philist(ii)
      nsq = nsqlist(ii)
      kappa = kappalist(ii)
      call calc_cov (phi,nsq,dm,F,betQ0,&
         lup,kappa,icf,n,p,T,TiF,FTF,Ups,ldh_Ups)
      if (nu .gt. 0d0) then
        do j = 1, Ntot
          call rchkusr
          lglk(j,ii) = jointyw_bi(n, wsample(:, j), y, l, Ups, ldh_Ups, &
             nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
        end do
      else
        do j = 1, Ntot
          call rchkusr
          lglk(j,ii) = jointyz_bi(n, wsample(:, j), y, l, Ups, ldh_Ups, &
             nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
        end do
      end if
    end do
  case default
    call rexit ("This function is unsed only for the binomial family.")
  end select
end subroutine llikfcnrb




! Local Variables:
! compile-command: "gfortran -c -fpic -Wunused-parameter -Wall \
!   -pedantic -o llik.o llik.f90"
! End:

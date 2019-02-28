!!! mcsp2.f90 ---
!!
!! Author: Evangelos Evangelou
!! Created: Tue, 15 Jul, 2014 16:59 (BST)
!! Last-Updated: Sun, 15 Nov, 2015 15:56 (GMT)
!!     Update #: 250
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mcspsample (lglk, z, z0, gmu, gmu0, &
   beta, ssq, phi, omg, kappa, acc, y, l, F, F0, &
   betm0, betQ0, ssqdf, ssqsc, phipars, omgpars, kappapars, &
   phisc, omgsc, kappasc, icf, &
   dft, tsq, dm, dmdm0, nch, Nmc, Nout, Nbi, Nthin, n, n0, p, ifam)

  use mcmcfcns
  use modelfcns
  use covfun, only: create_spcor
  implicit none
  integer, intent(in) :: nch, Nmc(nch), Nout, Nbi, Nthin, n, n0, p, ifam, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), F0(n0,p), &
     dft, tsq, dm(n,n), dmdm0(n,n0), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, phipars(4), omgpars(4), kappapars(2), phisc, omgsc, kappasc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), omg(Nout), &
     kappa(Nout)
  double precision, intent(out) :: z0(n0,Nout), beta(p,Nout), ssq(Nout), &
     lglk(Nout), gmu(n, Nout), gmu0(n0, Nout)
  integer, intent(out) :: acc(nch)
  integer i, j, k, ii
  double precision, parameter :: tsqdf = 0d0
  logical lnewcov
  double precision zUz, ldh_Ups, modeldfh, ssqdfsc, &
     respdf, tsqdfsc, tsqyy, betQm0(p), zmxi(n), z0_ups(n0), &
     T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), Upsz(n), TC(n,n0), FCTF(n0,p)

  acc = 0
  i = 0

  call create_model(ifam)
  call create_spcor(icf,n)

  do k = 1, nch
    i = i + 1

    call ini_mcmc(lglk(i),z(:,i),gmu(:,i),phi(i),omg(i),kappa(i),y,l,F,icf,dm,&
       betm0,betQ0,ssqdf,ssqsc,tsqdf,tsq,dft,n,p,ifam,betQm0,zmxi,T,&
       TiF,FTF,Ups,Upsz,zUz,ldh_Ups,modeldfh,ssqdfsc,respdf,tsqdfsc,tsqyy,&
       lnewcov)
    call rchkusr

    do j = 0, Nbi
      call sample_cov(lglk(i),phi(i),omg(i),kappa(i), &
         phipars,omgpars,kappapars,phisc,omgsc,kappasc,&
         dm,F,betQ0,n,p,acc(k),zmxi,T,TiF,FTF,Ups, &
         Upsz,lnewcov,zUz,ldh_Ups,modeldfh,ssqdfsc)
      call sample_ssq(ssq(i),modeldfh,zUz)
      call sample_z(lglk(i),z(:,i),gmu(:,i),y,l,dft,ssq(i),tsq,zmxi,Ups,&
         Upsz,zUz,modeldfh,n)
    end do
    call sample_beta(beta(:,i),z(:,i),ssq(i),n,p,betQm0,TiF,FTF)
    if (n0 .gt. 0) then
      call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),omg(i),&
         n0,n,p,dmdm0,F,F0,kappa(i),icf,T,z0_ups,TC,FCTF,lnewcov)
      gmu0(:,i) = invlink(z0(:,i), dft)
    end if
    call rchkusr
    do ii = 2, Nmc(k)
      i = i + 1
      lglk(i) = lglk(i-1)
      z(:,i) = z(:,i-1)
      gmu(:,i) = gmu(:,i-1)
      phi(i) = phi(i-1)
      omg(i) = omg(i-1)
      kappa(i) = kappa(i-1)
      do j = 1, Nthin
        call sample_cov(lglk(i),phi(i),omg(i),kappa(i),&
           phipars,omgpars,kappapars,phisc,omgsc,kappasc,&
           dm,F,betQ0,n,p,acc(k),zmxi,T,TiF,FTF,Ups, &
           Upsz,lnewcov,zUz,ldh_Ups,modeldfh,ssqdfsc)
        call sample_ssq(ssq(i),modeldfh,zUz)
        call sample_z(lglk(i),z(:,i),gmu(:,i),y,l,dft,ssq(i),tsq,zmxi,&
           Ups,Upsz,zUz,modeldfh,n)
      end do
      call sample_beta(beta(:,i),z(:,i),ssq(i),n,p,betQm0,TiF,FTF)
      if (n0 .gt. 0) then
        call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),omg(i),&
           n0,n,p,dmdm0,F,F0,kappa(i),icf,T,z0_ups,TC,FCTF,lnewcov)
        gmu0(:,i) = invlink(z0(:,i), dft)
      end if
      call rchkusr
    end do
    call end_mcmc
  end do
end subroutine mcspsample


subroutine samplemulti (lglk, z, gmu, phi, omg, y, l, F, &
   betm0, betQ0, ssqdf, ssqsc, kappa, icf, &
   nu, tsqdf, tsqin, dm, Ntot, Nout, Nbi, Nthin, n, p, kg, ifam)

  use mcmcfcns
  use modelfcns
  use covfun, only: create_spcor
  implicit none
  integer, intent(in) :: n, p, ifam, kg, icf, &
     Nbi(kg), Nthin(kg), Nout(kg), Ntot
  double precision, intent(in) :: y(n), l(n), F(n,p), &
     nu(kg), tsqdf, tsqin, dm(n,n), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, phi(kg), omg(kg), kappa(kg)
  double precision, intent(inout) :: z(n, Ntot)
  double precision, intent(out) :: lglk(Ntot), gmu(n, Ntot)
  double precision :: ssq, tsq
  integer i, ii, j, k
  logical lnewcov
  double precision zUz, ldh_Ups, modeldfh, ssqdfsc, &
     respdf, tsqdfsc, tsqyy, betQm0(p), zmxi(n), &
     T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), Upsz(n)

  i = 0

  call create_model (ifam)
  call create_spcor (icf,n)

  select case (ifam)
  case (0)
    do k = 1, kg
      i = i + 1
      call ini_mcmc(lglk(i),z(:,i),gmu(:,i),phi(k),omg(k),kappa(k),y,l,F,icf,&
         dm,betm0,betQ0,ssqdf,ssqsc,tsqdf,tsqin,nu(k),n,p,ifam,&
         betQm0,zmxi,T,TiF,FTF,Ups,Upsz,zUz,ldh_Ups,modeldfh,ssqdfsc,respdf,&
         tsqdfsc,tsqyy,lnewcov)
      call rchkusr
      do j = 0, Nbi(k)
        call sample_ssq(ssq,modeldfh,zUz)
        call sample_tsq(tsq,respdf,tsqyy)
        call samplez_gt(lglk(i),z(:,i),gmu(:,i),y,l,nu(k),ssq,zmxi,Ups,&
           Upsz,zUz,modeldfh,respdf,tsqyy,n)
      end do
      call rchkusr
      do ii = 2, Nout(k)
        i = i + 1
        lglk(i) = lglk(i-1)
        z(:,i) = z(:,i-1)
        gmu(:,i) = gmu(:,i-1)
        do j = 1, Nthin(k)
          call sample_ssq(ssq,modeldfh,zUz)
          call sample_tsq(tsq,respdf,tsqyy)
          call samplez_gt(lglk(i),z(:,i),gmu(:,i),y,l,nu(k),ssq,zmxi,&
             Ups,Upsz,zUz,modeldfh,respdf,tsqyy,n)
        end do
        call rchkusr
      end do
      call end_mcmc
    end do
  case default
    tsq = tsqin
    do k = 1, kg
      i = i + 1
      call ini_mcmc(lglk(i),z(:,i),gmu(:,i),phi(k),omg(k),kappa(k),y,l,F,icf,&
         dm,betm0,betQ0,ssqdf,ssqsc,tsqdf,tsqin,nu(k),n,p,ifam,&
         betQm0,zmxi,T,TiF,FTF,Ups,Upsz,zUz,ldh_Ups,modeldfh,ssqdfsc,respdf,&
         tsqdfsc,tsqyy,lnewcov)
      call rchkusr
      call create_model(ifam)
      do j = 0, Nbi(k)
        call sample_ssq(ssq,modeldfh,zUz)
        call sample_z(lglk(i),z(:,i),gmu(:,i),y,l,nu(k),ssq,tsq,zmxi,Ups,&
           Upsz,zUz,modeldfh,n)
      end do
      call rchkusr
      do ii = 2, Nout(k)
        i = i + 1
        lglk(i) = lglk(i-1)
        z(:,i) = z(:,i-1)
        gmu(:,i) = gmu(:,i-1)
        do j = 1, Nthin(k)
          call sample_ssq(ssq,modeldfh,zUz)
          call sample_z(lglk(i),z(:,i),gmu(:,i),y,l,nu(k),ssq,tsq,zmxi,&
             Ups,Upsz,zUz,modeldfh,n)
        end do
        call rchkusr
      end do
      call end_mcmc
    end do
  end select
end subroutine samplemulti


subroutine mcspsamtry (lglk, z, phi, omg, kappa, acc, y, l, F, &
   betm0, betQ0, ssqdf, ssqsc, phipars, omgpars, kappapars, &
   phisc, omgsc, kappasc, &
   icf, dft, tsq, dm, Nout, Npr, n, p, ifam)

  use mcmcfcns
  use modelfcns
  use msgmc
  use covfun, only: create_spcor
  implicit none
  integer, intent(in) :: Nout, n, p, ifam, Npr, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), &
     dft, tsq, dm(n,n), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, phipars(4), omgpars(4), kappapars(2), phisc, omgsc, kappasc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), omg(Nout), &
     kappa(Nout)
  double precision, intent(out) :: lglk(Nout)
  integer, intent(out) :: acc
  integer i, ia
  double precision lglk1, ssq, phi1, omg1, kappa1, z1(n), gmu(n)
!   integer, parameter :: nl = 19
!   character(len=nl) :: label
  integer iap
  double precision :: tsqdf
  logical lnewcov
  double precision zUz, ldh_Ups, modeldfh, ssqdfsc, &
     respdf, tsqdfsc, tsqyy, betQm0(p), zmxi(n), &
     T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), Upsz(n)

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
  omg1 = omg(1)
  kappa1 = kappa(1)

  call create_model (ifam)
  call create_spcor (icf,n)

  call ini_mcmc(lglk1,z1,gmu,phi1,omg1,kappa1,y,l,F,icf,&
     dm,betm0,betQ0,ssqdf,ssqsc,tsqdf,tsq,dft,n,p,ifam,&
     betQm0,zmxi,T,TiF,FTF,Ups,Upsz,zUz,ldh_Ups,modeldfh,ssqdfsc,respdf,&
     tsqdfsc,tsqyy,lnewcov)
  call rchkusr
  do i = 1, Nout
    call sample_cov(lglk1,phi1,omg1,kappa1,&
       phipars,omgpars,kappapars,phisc,omgsc,kappasc,dm,&
       F,betQ0,n,p,ia,zmxi,T,TiF,FTF,Ups, &
       Upsz,lnewcov,zUz,ldh_Ups,modeldfh,ssqdfsc)
    call sample_ssq(ssq,modeldfh,zUz)
    call sample_z(lglk1,z1,gmu,y,l,dft,ssq,tsq,zmxi,Ups,Upsz,zUz,modeldfh,n)
    lglk(i) = lglk1
    z(:,i) = z1
    phi(i) = phi1
    omg(i) = omg1
    kappa(i) = kappa1
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


subroutine trgasample (lglk, z, z0, gmu, gmu0, &
   beta, ssq, tsq, phi, omg, kappa, acc, y, l, F,&
   F0, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsqsc, phipars, omgpars, &
   kappapars, phisc, omgsc, kappasc, icf, dft, dm, dmdm0, nch, &
   Nmc, Nout, Nbi, Nthin, n, n0, p)

  use modelfcns
  use mcmcfcns
  use linkfcns, only: invlink_ga
  use covfun, only: create_spcor
  implicit none
  integer, intent(in) :: nch, Nmc(nch), Nout, Nbi, Nthin, n, n0, p, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), F0(n0,p), &
     dft, dm(n,n), dmdm0(n,n0), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, tsqdf, tsqsc, phipars(4), omgpars(4), kappapars(2), &
     phisc, omgsc, kappasc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), omg(Nout), &
     kappa(Nout)
  double precision, intent(out) :: z0(n0,Nout), beta(p,Nout), ssq(Nout), &
     tsq(Nout), lglk(Nout), gmu(n,Nout), gmu0(n0,Nout)
  integer, intent(out) :: acc(nch)
  integer i, j, k, ii
  integer, parameter :: ifam = 0
  logical lnewcov
  double precision zUz, ldh_Ups, modeldfh, ssqdfsc, &
     respdf, tsqdfsc, tsqyy, betQm0(p), zmxi(n), z0_ups(n0), &
     T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), Upsz(n), TC(n,n0), FCTF(n0,p)

  acc = 0
  i = 0

  call create_model(ifam)
  call create_spcor(icf,n)

  do k = 1, Nch
    i = i + 1

    call ini_mcmc(lglk(i),z(:,i),gmu(:,i),phi(i),omg(i),kappa(i),y,l,F,icf,&
       dm,betm0,betQ0,ssqdf,ssqsc,tsqdf,tsqsc,dft,n,p,ifam,&
       betQm0,zmxi,T,TiF,FTF,Ups,Upsz,zUz,ldh_Ups,modeldfh,ssqdfsc,respdf,&
       tsqdfsc,tsqyy,lnewcov)
    call rchkusr

    do j = 0, Nbi
      call sample_cov(lglk(i),phi(i),omg(i),kappa(i),&
         phipars,omgpars,kappapars,phisc,omgsc,kappasc,&
         dm,F,betQ0,n,p,acc(k),zmxi,T,TiF,FTF,Ups, &
         Upsz,lnewcov,zUz,ldh_Ups,modeldfh,ssqdfsc)
      call sample_ssq(ssq(i),modeldfh,zUz)
      call sample_tsq(tsq(i),respdf,tsqyy)
      call samplez_gt(lglk(i),z(:,i),gmu(:,i),y,l,dft,ssq(i),zmxi,Ups,&
         Upsz,zUz,modeldfh,respdf,tsqyy,n)
    end do
    call sample_beta(beta(:,i),z(:,i),ssq(i),n,p,betQm0,TiF,FTF)
    if (n0 .gt. 0) then
      call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),omg(i),&
         n0,n,p,dmdm0,F,F0,kappa(i),icf,T,z0_ups,TC,FCTF,lnewcov)
      gmu0(:,i) = invlink_ga(z0(:,i), dft)
    end if
    call rchkusr
    do ii = 2, Nmc(k)
      i = i + 1
      lglk(i) = lglk(i-1)
      z(:,i) = z(:,i-1)
      gmu(:,i) = gmu(:,i-1)
      phi(i) = phi(i-1)
      omg(i) = omg(i-1)
      kappa(i) = kappa(i-1)
      do j = 1, Nthin
        call sample_cov(lglk(i),phi(i),omg(i),kappa(i),&
           phipars,omgpars,kappapars,phisc,omgsc,kappasc,&
           dm,F,betQ0,n,p,acc(k),zmxi,T,TiF,FTF,Ups, &
           Upsz,lnewcov,zUz,ldh_Ups,modeldfh,ssqdfsc)
        call sample_ssq(ssq(i),modeldfh,zUz)
        call sample_tsq(tsq(i),respdf,tsqyy)
        call samplez_gt(lglk(i),z(:,i),gmu(:,i),y,l,dft,ssq(i),zmxi,&
           Ups,Upsz,zUz,modeldfh,respdf,tsqyy,n)
      end do
      call sample_beta(beta(:,i),z(:,i),ssq(i),n,p,betQm0,TiF,FTF)
      if (n0 .gt. 0) then
        call sample_z0(z0(:,i),z(:,i),beta(:,i),ssq(i),phi(i),omg(i),&
           n0,n,p,dmdm0,F,F0,kappa(i),icf,T,z0_ups,TC,FCTF,lnewcov)
        gmu0(:,i) = invlink_ga(z0(:,i), dft)
      end if
      call rchkusr
    end do

    call end_mcmc
  end do
end subroutine trgasample


subroutine trgasamtry (lglk, z, phi, omg, kappa, acc, y, l, F, &
   betm0, betQ0, ssqdf, ssqsc, tsqdf, tsqsc, phipars, omgpars, kappapars, &
   phisc, omgsc, kappasc, icf, dft, dm, Nout, Npr, n, p)

  use modelfcns
  use mcmcfcns
  use msgmc
  use covfun, only: create_spcor
  implicit none
  integer, intent(in) :: Nout, n, p, Npr, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), &
     dft, dm(n,n), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, tsqdf, tsqsc, phipars(4), omgpars(4), kappapars(2), &
     phisc, omgsc, kappasc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), omg(Nout), &
     kappa(Nout)
  double precision, intent(out) :: lglk(Nout)
  integer, intent(out) :: acc
  integer i, ia
  double precision ssq, tsq, phi1, omg1, z1(n), lglk1, gmu(n), kappa1
!   integer, parameter :: nl = 19
!   character(len=nl) :: label
  integer iap
  integer, parameter :: ifam = 0
  logical lnewcov
  double precision zUz, ldh_Ups, modeldfh, ssqdfsc, &
     respdf, tsqdfsc, tsqyy, betQm0(p), zmxi(n), &
     T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), Upsz(n)

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
  omg1 = omg(1)
  kappa1 = kappa(1)

  call create_model(ifam)
  call create_spcor(icf,n)

  call ini_mcmc(lglk1,z1,gmu,phi1,omg1,kappa1,y,l,F,icf,dm,&
     betm0,betQ0,ssqdf,ssqsc,tsqdf,tsqsc,dft,n,p,ifam,betQm0,zmxi,T,&
     TiF,FTF,Ups,Upsz,zUz,ldh_Ups,modeldfh,ssqdfsc,respdf,tsqdfsc,tsqyy,&
     lnewcov)
  call rchkusr

  do i = 1, Nout
    call sample_cov(lglk1,phi1,omg1,kappa1,&
       phipars,omgpars,kappapars,phisc,omgsc,kappasc,dm,&
       F,betQ0,n,p,ia,zmxi,T,TiF,FTF,Ups, &
         Upsz,lnewcov,zUz,ldh_Ups,modeldfh,ssqdfsc)
    call sample_ssq(ssq,modeldfh,zUz)
    call sample_tsq(tsq,respdf,tsqyy)
    call samplez_gt(lglk1,z1,gmu,y,l,dft,ssq,zmxi,Ups,Upsz,zUz,modeldfh,&
       respdf,tsqyy,n)
    lglk(i) = lglk1
    z(:,i) = z1
    phi(i) = phi1
    omg(i) = omg1
    kappa(i) = kappa1
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

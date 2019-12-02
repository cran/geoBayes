subroutine mcspsample_mala (lglk, z, z0, gmu, gmu0, &
   beta, ssq, phi, omg, kappa, acc, y, l, F, F0, &
   betm0, betQ0, ssqdf, ssqsc, phipars, omgpars, kappapars, &
   phisc, omgsc, kappasc, icf, &
   dft, tsq, dm, dmdm0, nch, Nmc, Nout, Nbi, Nthin, n, n0, p, ifam, z_eps, &
   acc_z)

  use mcmcfcns
  use modelfcns
  use covfun, only: create_spcor
  implicit none
  integer, intent(in) :: nch, Nmc(nch), Nout, Nbi, Nthin, n, n0, p, ifam, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), F0(n0,p), z_eps, &
     dft, tsq, dm(n,n), dmdm0(n,n0), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, phipars(4), omgpars(4), kappapars(2), phisc, omgsc, kappasc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), omg(Nout), &
     kappa(Nout)
  double precision, intent(out) :: z0(n0,Nout), beta(p,Nout), ssq(Nout), &
     lglk(Nout), gmu(n, Nout), gmu0(n0, Nout)
  integer, intent(out) :: acc(nch), acc_z(nch)
  integer i, j, k, ii
  double precision, parameter :: tsqdf = 0d0
  logical lnewcov
  double precision zUz, ldh_Ups, modeldfh, ssqdfsc, &
     respdf, tsqdfsc, tsqyy, betQm0(p), zmxi(n), z0_ups(n0), &
     T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), Upsz(n), TC(n,n0), FCTF(n0,p)

  acc = 0
  acc_z = 0
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
      call sample_z_mala(lglk(i),z(:,i),gmu(:,i),y,l,dft,ssq(i),tsq,zmxi,Ups,&
         Upsz,zUz,modeldfh,n,z_eps,acc_z(k))
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
        call sample_z_mala(lglk(i),z(:,i),gmu(:,i),y,l,dft,ssq(i),tsq,zmxi,&
           Ups,Upsz,zUz,modeldfh,n,z_eps,acc_z(k))
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
end subroutine mcspsample_mala

subroutine mcspsamtry_mala (lglk, z, phi, omg, kappa, acc, y, l, F, &
   betm0, betQ0, ssqdf, ssqsc, phipars, omgpars, kappapars, &
   phisc, omgsc, kappasc, &
   icf, dft, tsq, dm, Nout, Npr, n, p, ifam, z_eps, acc_z)

  use mcmcfcns
  use modelfcns
  use msgmc
  use covfun, only: create_spcor
  implicit none
  integer, intent(in) :: Nout, n, p, ifam, Npr, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), z_eps, &
     dft, tsq, dm(n,n), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, phipars(4), omgpars(4), kappapars(2), phisc, omgsc, kappasc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), omg(Nout), &
     kappa(Nout)
  double precision, intent(out) :: lglk(Nout)
  integer, intent(out) :: acc, acc_z
  integer i, ia, ib
  double precision lglk1, ssq, phi1, omg1, kappa1, z1(n), gmu(n)
!   integer, parameter :: nl = 19
!   character(len=nl) :: label
  integer iap, ibp
  double precision :: tsqdf
  logical lnewcov
  double precision zUz, ldh_Ups, modeldfh, ssqdfsc, &
     respdf, tsqdfsc, tsqyy, betQm0(p), zmxi(n), &
     T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), Upsz(n)

  tsqdf = 0d0

  call msgmca2
  call msgmcl2
  acc = 0
  acc_z = 0
  ia = 0
  ib = 0
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
    call sample_z_mala(lglk1,z1,gmu,y,l,dft,ssq,tsq,zmxi,Ups,Upsz,zUz,modeldfh,&
       n,z_eps,ib)
    lglk(i) = lglk1
    z(:,i) = z1
    phi(i) = phi1
    omg(i) = omg1
    kappa(i) = kappa1
    if (Npr .gt. 0 .and. mod(i,Npr) .eq. 0) then
      iap = 100*ia/Npr
      ibp = 100*ib/Npr
      call msgmci2 (i, iap, ibp)
      acc = acc + ia
      acc_z = acc_z + ib
      ia = 0
      ib = 0
      call rchkusr
    end if
  end do

  acc = acc + ia
  acc_z = acc_z + ib
  call end_mcmc
  call msgmcl2
  iap = 100*acc/Nout
  ibp = 100*acc_z/Nout
  call msgmce2 (iap, ibp)
  call msgmcl2
end subroutine mcspsamtry_mala


subroutine trgasample_mala (lglk, z, z0, gmu, gmu0, &
   beta, ssq, tsq, phi, omg, kappa, acc, y, l, F,&
   F0, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsqsc, phipars, omgpars, &
   kappapars, phisc, omgsc, kappasc, icf, dft, dm, dmdm0, nch, &
   Nmc, Nout, Nbi, Nthin, n, n0, p, z_eps, acc_z)

  use modelfcns
  use mcmcfcns
  use covfun, only: create_spcor
  implicit none
  integer, intent(in) :: nch, Nmc(nch), Nout, Nbi, Nthin, n, n0, p, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), F0(n0,p), z_eps, &
     dft, dm(n,n), dmdm0(n,n0), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, tsqdf, tsqsc, phipars(4), omgpars(4), kappapars(2), &
     phisc, omgsc, kappasc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), omg(Nout), &
     kappa(Nout)
  double precision, intent(out) :: z0(n0,Nout), beta(p,Nout), ssq(Nout), &
     tsq(Nout), lglk(Nout), gmu(n,Nout), gmu0(n0,Nout)
  integer, intent(out) :: acc(nch), acc_z(nch)
  integer i, j, k, ii
  integer, parameter :: ifam = 0
  logical lnewcov
  double precision zUz, ldh_Ups, modeldfh, ssqdfsc, &
     respdf, tsqdfsc, tsqyy, betQm0(p), zmxi(n), z0_ups(n0), &
     T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), Upsz(n), TC(n,n0), FCTF(n0,p)

  acc = 0
  acc_z = 0
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
      call samplez_gt_mala(lglk(i),z(:,i),gmu(:,i),y,l,dft,ssq(i),zmxi,Ups,&
         Upsz,zUz,modeldfh,respdf,tsqyy,n,z_eps,acc_z(k))
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
        call sample_tsq(tsq(i),respdf,tsqyy)
        call samplez_gt_mala(lglk(i),z(:,i),gmu(:,i),y,l,dft,ssq(i),zmxi,&
           Ups,Upsz,zUz,modeldfh,respdf,tsqyy,n,z_eps,acc_z(k))
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
end subroutine trgasample_mala


subroutine trgasamtry_mala (lglk, z, phi, omg, kappa, acc, y, l, F, &
   betm0, betQ0, ssqdf, ssqsc, tsqdf, tsqsc, phipars, omgpars, kappapars, &
   phisc, omgsc, kappasc, icf, dft, dm, Nout, Npr, n, p, z_eps, acc_z)

  use modelfcns
  use mcmcfcns
  use msgmc
  use covfun, only: create_spcor
  implicit none
  integer, intent(in) :: Nout, n, p, Npr, icf
  double precision, intent(in) :: y(n), l(n), F(n,p), z_eps, &
     dft, dm(n,n), betm0(p), betQ0(p,p), ssqdf, &
     ssqsc, tsqdf, tsqsc, phipars(4), omgpars(4), kappapars(2), &
     phisc, omgsc, kappasc
  double precision, intent(inout) :: z(n,Nout), phi(Nout), omg(Nout), &
     kappa(Nout)
  double precision, intent(out) :: lglk(Nout)
  integer, intent(out) :: acc, acc_z
  integer i, ia, ib
  double precision ssq, tsq, phi1, omg1, z1(n), lglk1, gmu(n), kappa1
!   integer, parameter :: nl = 19
!   character(len=nl) :: label
  integer iap, ibp
  integer, parameter :: ifam = 0
  logical lnewcov
  double precision zUz, ldh_Ups, modeldfh, ssqdfsc, &
     respdf, tsqdfsc, tsqyy, betQm0(p), zmxi(n), &
     T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), Upsz(n)

  call msgmca2
  call msgmcl2
  acc = 0
  ia = 0
  ib = 0
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
    call samplez_gt_mala(lglk1,z1,gmu,y,l,dft,ssq,zmxi,Ups,Upsz,zUz,modeldfh,&
       respdf,tsqyy,n,z_eps,ib)
    lglk(i) = lglk1
    z(:,i) = z1
    phi(i) = phi1
    omg(i) = omg1
    kappa(i) = kappa1
    if (Npr .gt. 0 .and. mod(i,Npr) .eq. 0) then
      iap = 100*ia/Npr
      ibp = 100*ib/Npr
      call msgmci2 (i, iap, ibp)
      acc = acc + ia
      acc_z = acc_z + ib
      ia = 0
      ib = 0
      call rchkusr
    end if
  end do
  acc = acc + ia
  call end_mcmc
  call msgmcl2
  iap = 100*acc/Nout
  ibp = 100*acc_z/Nout
  call msgmce2 (iap, ibp)
  call msgmcl2
end subroutine trgasamtry_mala

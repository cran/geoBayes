!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Commentary: Derivatives of BF
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcbd_no (bfact, bfdnu, bfdphi, bfdnsq, bfdkappa, &
   phi, nu, nsq, kappa, icf, Ntot, sample, weights, QRin, &
   n, p, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr, icv)

  use modelfcns, logpdfzf => logpdfz, logpdfydmu => logpdfydlnk
  use covfun
  use betaprior
  use calcbd_fcns
  implicit none
  integer, intent(in) :: n, p, Ntot, ifam, icf, itr(n), icv
  double precision, intent(in) :: phi, nsq, &
     kappa, nu, sample(n,Ntot), weights(Ntot), &
     betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n,p), dm(n,n)
  double precision, intent(in) :: QRin(Ntot)
  double precision, intent(out) :: bfact, bfdphi, bfdnu, bfdnsq, bfdkappa
  logical lmxi
  double precision ssqdfsc, respdfh, modeldfh, tsqval
  integer j
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, lglk(Ntot), xi(n)
  double precision zsam(n), msam(n), jsam(n), sam(n), dlogpy(n), dzdnu(n), &
     dmudnu(n), djsam(n), lpdfz, lpdfy
  double precision dy, dz, dj, dlglkdnu(Ntot), dlglkdphi(Ntot), &
     dlglkdnsq(Ntot), dlglkdkappa(Ntot)
  double precision zUz, zUDTUz, DT(n,n,3), trUpsDTh(3), Upsz(n), dfozUz
  double precision ymm(n), ymml(n)

  abstract interface
    function condymu_ (n, y1, y2, mu, tsqval, respdfh)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
         respdfh
      double precision condymu_
    end function condymu_
  end interface

  procedure(condymu_), pointer :: condymuf => null()

  abstract interface
    double precision function weigh_llik_ (llik, lw, q, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: llik(n), lw(n), q(n)
    end function weigh_llik_
  end interface

  abstract interface
    double precision function weigh_llik_deriv_ (dllik, llik, lw, q, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: dllik(n), llik(n), lw(n), q(n)
    end function weigh_llik_deriv_
  end interface

  procedure (weigh_llik_), pointer :: weigh_llik
  procedure (weigh_llik_deriv_), pointer :: weigh_llik_deriv

  select case (icv)
  case (0)
    weigh_llik => weigh_llik_st
    weigh_llik_deriv => weigh_llik_deriv_st
  case default
    weigh_llik => weigh_llik_cv
    weigh_llik_deriv => weigh_llik_deriv_cv
  end select

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
    condymuf => condymu_gt
  case default
    tsqval = tsq
    condymuf => condymu_sp
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  call calc_cov (phi,nsq,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
  do j = 1, 3
    DT(:,:,j) = cor_dcov(n, dm, phi, nsq, kappa, j)
    trUpsDTh(j) = .5d0*traceAB(Ups,DT(:,:,j),n)
  end do

  do j = 1, Ntot
    call rchkusr
    sam = sample(:,j)
    zsam = sam
    msam = invlink(zsam,nu)
    dlogpy = logpdfydmu(y,l,msam)/tsq
    dmudnu = invlinkdn(zsam,nu)
    if (lmxi) zsam = zsam - xi
    call dsymv ('u',n,1d0,Ups,n,zsam,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
    zUz = dot_product(zsam,Upsz) + ssqdfsc
    lpdfz = ldh_Ups - modeldfh*log(zUz)
    dfozUz = modeldfh/zUz
    if (ifam .eq. 0) then
      ymm = y - msam
      ymml = l*ymm
      lpdfy = condymuf(n,y,l,msam,tsqval,respdfh)
      dy = dot_product(ymm,ymml) + tsqval
      dlogpy = (2*respdfh/dy)*ymml
      dy = dot_product(dlogpy,dmudnu)
    else
      lpdfy = condymuf(n,y,l,msam,tsqval,respdfh)
      dy = dot_product(dlogpy,dmudnu)
    end if
    lglk(j) = lpdfz + lpdfy
    dlglkdnu(j) = dy
    zUDTUz = qform(Upsz,DT(:,:,1),n)
    dlglkdphi(j) = -trUpsDTh(1) + dfozUz*zUDTUz
    zUDTUz = qform(Upsz,DT(:,:,2),n)
    dlglkdnsq(j) = -trUpsDTh(2) + dfozUz*zUDTUz
    zUDTUz = qform(Upsz,DT(:,:,3),n)
    dlglkdkappa(j) = -trUpsDTh(3) + dfozUz*zUDTUz
  end do
  bfact = weigh_llik(lglk, weights, QRin, Ntot)
  bfdnu = weigh_llik_deriv(dlglkdnu, lglk, weights, QRin, Ntot)
  bfdphi = weigh_llik_deriv(dlglkdphi, lglk, weights, QRin, Ntot)
  bfdnsq = weigh_llik_deriv(dlglkdnsq, lglk, weights, QRin, Ntot)
  bfdkappa = weigh_llik_deriv(dlglkdkappa, lglk, weights, QRin, Ntot)
end subroutine calcbd_no

subroutine calcbd_mu (bfact, bfdnu, bfdphi, bfdnsq, bfdkappa, &
   phi, nu, nsq, kappa, icf, Ntot, sample, weights, QRin, &
   n, p, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr, icv)

  use modelfcns, logpdfzf => logpdfz, logpdfydmu => logpdfydlnk
  use covfun
  use betaprior
  use calcbd_fcns
  implicit none
  integer, intent(in) :: n, p, Ntot, ifam, icf, itr(n), icv
  double precision, intent(in) :: phi, nsq, &
     kappa, nu, sample(n,Ntot), weights(Ntot), &
     betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n,p), dm(n,n)
  double precision, intent(in) :: QRin(Ntot)
  double precision, intent(out) :: bfact, bfdphi, bfdnu, bfdnsq, bfdkappa
  logical lmxi
  double precision ssqdfsc, respdfh, modeldfh, tsqval
  integer j
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, lglk(Ntot), xi(n)
  double precision zsam(n), msam(n), jsam(n), sam(n), dlogpy(n), dzdnu(n), &
     dmudnu(n), djsam(n), lpdfz, lpdfy
  double precision dy, dz, dj, dlglkdnu(Ntot), dlglkdphi(Ntot), &
     dlglkdnsq(Ntot), dlglkdkappa(Ntot)
  double precision zUz, zUDTUz, DT(n,n,3), trUpsDTh(3), Upsz(n), dfozUz
  double precision ymm(n), ymml(n)

  abstract interface
    function condymu_ (n, y1, y2, mu, tsqval, respdfh)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
         respdfh
      double precision condymu_
    end function condymu_
  end interface

  procedure(condymu_), pointer :: condymuf => null()

  abstract interface
    double precision function weigh_llik_ (llik, lw, q, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: llik(n), lw(n), q(n)
    end function weigh_llik_
  end interface

  abstract interface
    double precision function weigh_llik_deriv_ (dllik, llik, lw, q, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: dllik(n), llik(n), lw(n), q(n)
    end function weigh_llik_deriv_
  end interface

  procedure (weigh_llik_), pointer :: weigh_llik
  procedure (weigh_llik_deriv_), pointer :: weigh_llik_deriv

  select case (icv)
  case (0)
    weigh_llik => weigh_llik_st
    weigh_llik_deriv => weigh_llik_deriv_st
  case default
    weigh_llik => weigh_llik_cv
    weigh_llik_deriv => weigh_llik_deriv_cv
  end select

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
    condymuf => condymu_gt
  case default
    tsqval = tsq
    condymuf => condymu_sp
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  call calc_cov (phi,nsq,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
  do j = 1, 3
    DT(:,:,j) = cor_dcov(n, dm, phi, nsq, kappa, j)
    trUpsDTh(j) = .5d0*traceAB(Ups,DT(:,:,j),n)
  end do

  do j = 1, Ntot
    call rchkusr
    sam = sample(:,j)
    msam = sam
    zsam = flink(msam,nu)
    jsam = logilinkdz(zsam,nu)
    dlogpy = 0d0
    dzdnu = -invlinkdn(zsam,nu)/invlinkdz(zsam,nu)
    dmudnu = 0d0
    djsam = (invlinkdzdn(zsam,nu) + invlinkhz(zsam,nu)*dzdnu) &
       /invlinkdz(zsam,nu)
    if (lmxi) zsam = zsam - xi
    call dsymv ('u',n,1d0,Ups,n,zsam,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
    zUz = dot_product(zsam,Upsz) + ssqdfsc
    lpdfz = ldh_Ups - modeldfh*log(zUz)
    dfozUz = modeldfh/zUz
    dz = -(2*dfozUz)*dot_product(Upsz,dzdnu)
    dj = -sum(djsam)
    lpdfy = 0d0
    dy = 0d0
    lglk(j) = lpdfz + lpdfy - sum(jsam)
    dlglkdnu(j) = dy + dz + dj
    zUDTUz = qform(Upsz,DT(:,:,1),n)
    dlglkdphi(j) = -trUpsDTh(1) + dfozUz*zUDTUz
    zUDTUz = qform(Upsz,DT(:,:,2),n)
    dlglkdnsq(j) = -trUpsDTh(2) + dfozUz*zUDTUz
    zUDTUz = qform(Upsz,DT(:,:,3),n)
    dlglkdkappa(j) = -trUpsDTh(3) + dfozUz*zUDTUz
  end do
  bfact = weigh_llik(lglk, weights, QRin, Ntot)
  bfdnu = weigh_llik_deriv(dlglkdnu, lglk, weights, QRin, Ntot)
  bfdphi = weigh_llik_deriv(dlglkdphi, lglk, weights, QRin, Ntot)
  bfdnsq = weigh_llik_deriv(dlglkdnsq, lglk, weights, QRin, Ntot)
  bfdkappa = weigh_llik_deriv(dlglkdkappa, lglk, weights, QRin, Ntot)
end subroutine calcbd_mu


subroutine calcbd_wo (bfact, bfdnu, bfdphi, bfdnsq, bfdkappa, &
   phi, nu, nsq, kappa, icf, Ntot, sample, weights, QRin, &
   n, p, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr, icv)

  use modelfcns, logpdfzf => logpdfz, logpdfydmu => logpdfydlnk
  use covfun
  use betaprior
  use calcbd_fcns
  implicit none
  integer, intent(in) :: n, p, Ntot, ifam, icf, itr(n), icv
  double precision, intent(in) :: phi, nsq, &
     kappa, nu, sample(n,Ntot), weights(Ntot), &
     betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n,p), dm(n,n)
  double precision, intent(in) :: QRin(Ntot)
  double precision, intent(out) :: bfact, bfdphi, bfdnu, bfdnsq, bfdkappa
  logical lmxi
  double precision ssqdfsc, respdfh, modeldfh, tsqval
  integer j
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, lglk(Ntot), xi(n)
  double precision zsam(n), msam(n), jsam(n), sam(n), dlogpy(n), dzdnu(n), &
     dmudnu(n), djsam(n), lpdfz, lpdfy
  double precision dy, dz, dj, dlglkdnu(Ntot), dlglkdphi(Ntot), &
     dlglkdnsq(Ntot), dlglkdkappa(Ntot)
  double precision zUz, zUDTUz, DT(n,n,3), trUpsDTh(3), Upsz(n), dfozUz
  double precision ymm(n), ymml(n)

  abstract interface
    function condymu_ (n, y1, y2, mu, tsqval, respdfh)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
         respdfh
      double precision condymu_
    end function condymu_
  end interface

  procedure(condymu_), pointer :: condymuf => null()

  abstract interface
    double precision function weigh_llik_ (llik, lw, q, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: llik(n), lw(n), q(n)
    end function weigh_llik_
  end interface

  abstract interface
    double precision function weigh_llik_deriv_ (dllik, llik, lw, q, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: dllik(n), llik(n), lw(n), q(n)
    end function weigh_llik_deriv_
  end interface

  procedure (weigh_llik_), pointer :: weigh_llik
  procedure (weigh_llik_deriv_), pointer :: weigh_llik_deriv

  select case (icv)
  case (0)
    weigh_llik => weigh_llik_st
    weigh_llik_deriv => weigh_llik_deriv_st
  case default
    weigh_llik => weigh_llik_cv
    weigh_llik_deriv => weigh_llik_deriv_cv
  end select

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
    condymuf => condymu_gt
  case default
    tsqval = tsq
    condymuf => condymu_sp
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  call calc_cov (phi,nsq,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
  do j = 1, 3
    DT(:,:,j) = cor_dcov(n, dm, phi, nsq, kappa, j)
    trUpsDTh(j) = .5d0*traceAB(Ups,DT(:,:,j),n)
  end do

  do j = 1, Ntot
    call rchkusr
    sam = sample(:,j)
    zsam = transfw(sam,nu)
    msam = invlink(zsam,nu)
    jsam = logitrwdz(zsam,nu)
    dlogpy = logpdfydmu(y,l,msam)/tsq
    dzdnu = -invtrwdn(zsam,nu)/invtrwdz(zsam,nu)
    dmudnu = invlinkdn(zsam,nu) + invlinkdz(zsam,nu)*dzdnu
    djsam = (invtrwdzdn(zsam,nu) + invtrwhz(zsam,nu)*dzdnu)/invtrwdz(zsam,nu)
    if (lmxi) zsam = zsam - xi
    call dsymv ('u',n,1d0,Ups,n,zsam,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
    zUz = dot_product(zsam,Upsz) + ssqdfsc
    lpdfz = ldh_Ups - modeldfh*log(zUz)
    dfozUz = modeldfh/zUz
    dz = -(2*dfozUz)*dot_product(Upsz,dzdnu)
    dj = -sum(djsam)
    if (ifam .eq. 0) then
      ymm = y - msam
      ymml = l*ymm
      lpdfy = condymuf(n,y,l,msam,tsqval,respdfh)
      dy = dot_product(ymm,ymml) + tsqval
      dlogpy = (2*respdfh/dy)*ymml
      dy = dot_product(dlogpy,dmudnu)
    else
      lpdfy = condymuf(n,y,l,msam,tsqval,respdfh)
      dy = dot_product(dlogpy,dmudnu)
    end if
    lglk(j) = lpdfz + lpdfy - sum(jsam)
    dlglkdnu(j) = dy + dz + dj
    zUDTUz = qform(Upsz,DT(:,:,1),n)
    dlglkdphi(j) = -trUpsDTh(1) + dfozUz*zUDTUz
    zUDTUz = qform(Upsz,DT(:,:,2),n)
    dlglkdnsq(j) = -trUpsDTh(2) + dfozUz*zUDTUz
    zUDTUz = qform(Upsz,DT(:,:,3),n)
    dlglkdkappa(j) = -trUpsDTh(3) + dfozUz*zUDTUz
  end do
  bfact = weigh_llik(lglk, weights, QRin, Ntot)
  bfdnu = weigh_llik_deriv(dlglkdnu, lglk, weights, QRin, Ntot)
  bfdphi = weigh_llik_deriv(dlglkdphi, lglk, weights, QRin, Ntot)
  bfdnsq = weigh_llik_deriv(dlglkdnsq, lglk, weights, QRin, Ntot)
  bfdkappa = weigh_llik_deriv(dlglkdkappa, lglk, weights, QRin, Ntot)
end subroutine calcbd_wo


subroutine calcbd_tr (bfact, bfdnu, bfdphi, bfdnsq, bfdkappa, &
   phi, nu, nsq, kappa, icf, Ntot, sample, weights, QRin, &
   n, p, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, y, l, F, dm, ifam, itr, icv)

  use modelfcns, logpdfzf => logpdfz, logpdfydmu => logpdfydlnk
  use covfun
  use betaprior
  use calcbd_fcns
  implicit none
  integer, intent(in) :: n, p, Ntot, ifam, icf, itr(n), icv
  double precision, intent(in) :: phi, nsq, &
     kappa, nu, sample(n,Ntot), weights(Ntot), &
     betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, y(n), l(n), &
     F(n,p), dm(n,n)
  double precision, intent(in) :: QRin(Ntot)
  double precision, intent(out) :: bfact, bfdphi, bfdnu, bfdnsq, bfdkappa
  logical lmxi
  double precision ssqdfsc, respdfh, modeldfh, tsqval
  integer j
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, lglk(Ntot), xi(n)
  double precision zsam(n), msam(n), jsam(n), sam(n), dlogpy(n), dzdnu(n), &
     dmudnu(n), djsam(n), lpdfz, lpdfy
  double precision dy, dz, dj, dlglkdnu(Ntot), dlglkdphi(Ntot), &
     dlglkdnsq(Ntot), dlglkdkappa(Ntot)
  double precision zUz, zUDTUz, DT(n,n,3), trUpsDTh(3), Upsz(n), dfozUz
  double precision ymm(n), ymml(n)
  logical ltr0(n), ltr1(n), ltr2(n), ltrmu

  abstract interface
    function condymu_ (n, y1, y2, mu, tsqval, respdfh)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
         respdfh
      double precision condymu_
    end function condymu_
  end interface

  procedure(condymu_), pointer :: condymuf => null()

  abstract interface
    double precision function weigh_llik_ (llik, lw, q, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: llik(n), lw(n), q(n)
    end function weigh_llik_
  end interface

  abstract interface
    double precision function weigh_llik_deriv_ (dllik, llik, lw, q, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: dllik(n), llik(n), lw(n), q(n)
    end function weigh_llik_deriv_
  end interface

  procedure (weigh_llik_), pointer :: weigh_llik
  procedure (weigh_llik_deriv_), pointer :: weigh_llik_deriv

  select case (icv)
  case (0)
    weigh_llik => weigh_llik_st
    weigh_llik_deriv => weigh_llik_deriv_st
  case default
    weigh_llik => weigh_llik_cv
    weigh_llik_deriv => weigh_llik_deriv_cv
  end select

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
    condymuf => condymu_gt
  case default
    tsqval = tsq
    condymuf => condymu_sp
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  call calc_cov (phi,nsq,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
  do j = 1, 3
    DT(:,:,j) = cor_dcov(n, dm, phi, nsq, kappa, j)
    trUpsDTh(j) = .5d0*traceAB(Ups,DT(:,:,j),n)
  end do

  ltr0 = itr .eq. 0
  ltr1 = itr .eq. 1
  ltr2 = itr .eq. 2
  ltrmu = all(ltr1)

  do j = 1, Ntot
    call rchkusr
    sam = sample(:,j)
    where (ltr0)
      zsam = sam
      msam = invlink(zsam,nu)
      jsam = 0d0
      dlogpy = logpdfydmu(y,l,msam)/tsq
      dzdnu = 0d0
      dmudnu = invlinkdn(zsam,nu)
      djsam = 0d0
    elsewhere (ltr1)
      msam = sam
      zsam = flink(msam,nu)
      jsam = logilinkdz(zsam,nu)
      dlogpy = 0d0
      dzdnu = -invlinkdn(zsam,nu)/invlinkdz(zsam,nu)
      dmudnu = 0d0
      djsam = (invlinkdzdn(zsam,nu) + invlinkhz(zsam,nu)*dzdnu) &
         /invlinkdz(zsam,nu)
    elsewhere (ltr2)
      zsam = transfw(sam,nu)
      msam = invlink(zsam,nu)
      jsam = logitrwdz(zsam,nu)
      dlogpy = logpdfydmu(y,l,msam)/tsq
      dzdnu = -invtrwdn(zsam,nu)/invtrwdz(zsam,nu)
      dmudnu = invlinkdn(zsam,nu) + invlinkdz(zsam,nu)*dzdnu
      djsam = (invtrwdzdn(zsam,nu) + invtrwhz(zsam,nu)*dzdnu)/invtrwdz(zsam,nu)
    end where
    if (lmxi) zsam = zsam - xi
    call dsymv ('u',n,1d0,Ups,n,zsam,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
    zUz = dot_product(zsam,Upsz) + ssqdfsc
    lpdfz = ldh_Ups - modeldfh*log(zUz)
    dfozUz = modeldfh/zUz
    dz = -(2*dfozUz)*dot_product(Upsz,dzdnu)
    dj = -sum(djsam)
    if (ltrmu) then
      lpdfy = 0d0
      dy = 0d0
    else
      if (ifam .eq. 0) then
        ymm = y - msam
        ymml = l*ymm
        lpdfy = condymuf(n,y,l,msam,tsqval,respdfh)
        dy = dot_product(ymm,ymml) + tsqval
        dlogpy = (2*respdfh/dy)*ymml
        dy = dot_product(dlogpy,dmudnu)
      else
        lpdfy = condymuf(n,y,l,msam,tsqval,respdfh)
        dy = dot_product(dlogpy,dmudnu)
      end if
    end if
    lglk(j) = lpdfz + lpdfy - sum(jsam)
    dlglkdnu(j) = dy + dz + dj
    zUDTUz = qform(Upsz,DT(:,:,1),n)
    dlglkdphi(j) = -trUpsDTh(1) + dfozUz*zUDTUz
    zUDTUz = qform(Upsz,DT(:,:,2),n)
    dlglkdnsq(j) = -trUpsDTh(2) + dfozUz*zUDTUz
    zUDTUz = qform(Upsz,DT(:,:,3),n)
    dlglkdkappa(j) = -trUpsDTh(3) + dfozUz*zUDTUz
  end do
  bfact = weigh_llik(lglk, weights, QRin, Ntot)
  bfdnu = weigh_llik_deriv(dlglkdnu, lglk, weights, QRin, Ntot)
  bfdphi = weigh_llik_deriv(dlglkdphi, lglk, weights, QRin, Ntot)
  bfdnsq = weigh_llik_deriv(dlglkdnsq, lglk, weights, QRin, Ntot)
  bfdkappa = weigh_llik_deriv(dlglkdkappa, lglk, weights, QRin, Ntot)
end subroutine calcbd_tr

subroutine llikfcn_dcov_tr (dlglk, phi, omg, nu, kappa, &
   sample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, ifam, itr)
!! Log-likelihood function derivative w.r.t covariance parameters.
!! itr is the type of transformation used: 0 = no, 1 = mu, 2 = wo.

  use modelfcns
  use covfun
  use betaprior
  implicit none
  integer, intent(in) :: n, p, ifam, Ntot, icf, itr(n)
  double precision, intent(in) :: phi, omg, nu, &
     sample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappa
  double precision, intent(out) :: dlglk(3,Ntot)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, tsqval, respdfh, xi(n)
  integer i, j
  double precision zsam(n), zU(n), zUz, zUDTUz, DT(n,n,3), trUpsDT(3)

  call create_model (ifam)
  call create_spcor (icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)

  do i = 1, 3
    DT(:,:,i) = cor_dcov(n, dm, phi, omg, kappa, i)
    trUpsDT(i) = traceAB(Ups,DT(:,:,i),n)
  end do

  do j = 1, Ntot
    call rchkusr
    where (itr .eq. 0)
      zsam = sample(:,j)
    elsewhere (itr .eq. 1)
      zsam = flink(sample(:,j),nu)
    elsewhere (itr .eq. 2)
      zsam = transfw(sample(:,j),nu)
    end where
    if (lmxi) zsam = zsam - xi
    call dsymv ('u',n,1d0,Ups,n,zsam,1,0d0,zU,1)
    zUz = dot_product(zsam,zU) + ssqdfsc
    do i = 1, 3
      zUDTUz = qform(zU,DT(:,:,i),n)
      dlglk(i,j) = -.5d0*trUpsDT(i) + modeldfh*zUDTUz/zUz
    end do
  end do

contains

  double precision function qform(v,A,n)
    ! Compute v'*A*v where A is symmetric with only the upper triangle
    ! specified.
    integer, intent(in) :: n
    double precision, intent(in) :: v(n), A(n,n)
    double precision Av(n)
    call dsymv ('u',n,1d0,A,n,v,1,0d0,Av,1)
    qform = dot_product(Av,v)
  end function qform

  double precision function traceAB(A,B,n)
    ! Compute the trace(matmul(A,B)) where A,B are symmetric matrices with
    ! only their upper triangular elements specified.
    integer, intent(in) :: n
    double precision, intent(in) :: A(n,n), B(n,n)
    integer j
    double precision c
    traceAB = A(1,1)*B(1,1)
    do j = 2, n
      c = dot_product(A(:j-1,j),B(:j-1,j))
      traceAB = traceAB + c + c + A(j,j)*B(j,j)
    end do
  end function traceAB

  double precision function cor_dcov (n, dm, phi, omg, kappa, id) result (DT)
    ! id = which derivative do you want? phi, omg, kappa
    ! Return dT where T = omg*I + R(phi,kappa)
    use covfun, only: spcor_dh, spcor_dk, fill_symmetric_matrix, upper_tri
    implicit none
    integer, intent(in) :: n, id
    double precision, intent(in) :: dm(n,n), phi, omg, kappa
    !double precision, intent(out) :: DT(n,n)
    dimension :: DT(n,n)
    integer i
    double precision H(n,n)
    select case (id)
    case (1)
      if (phi .eq. 0d0) then
        DT = 0d0
        return
      else
        where (upper_tri())
          H = dm/phi
          DT = -(H/phi)*spcor_dh(H,kappa)
        end where
      end if
    case (2)
      DT = 0d0
      do i = 1, n
        DT(i,i) = 1d0
      end do
    case (3)
      if (phi .eq. 0d0) then
        DT = 0d0
        return
      else
        where (upper_tri())
          H = dm/phi
          DT = spcor_dk(H,kappa)
        end where
      end if
    end select
  end function cor_dcov
end subroutine llikfcn_dcov_tr

subroutine llikfcn_hcov_tr (hlglk, phi, omg, nu, kappa, &
   sample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, ifam, itr)
!! Log-likelihood function derivative w.r.t covariance parameters.
!! itr is the type of transformation used: 0 = no, 1 = mu, 2 = wo.

  use modelfcns
  use covfun
  use betaprior
  implicit none
  integer, intent(in) :: n, p, ifam, Ntot, icf, itr(n)
  double precision, intent(in) :: phi, omg, nu, &
     sample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappa
  double precision, intent(out) :: hlglk(3,3,Ntot)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, tsqval, respdfh, xi(n)
  integer i, j, k
  double precision zsam(n), zU(n), zUz, zUDTUz(3), zUHTUz, zUDTUDTUz
  double precision DT(n,n,3), HT(n,n,3,3), trUpsHT(3,3), trUpsDTUpsDT(3,3)
  double precision UpsDT(n,n), DTUpsDT(n,n)

  call create_model (ifam)
  call create_spcor (icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
  call fill_symmetric_matrix(Ups,n)

  do i = 1, 3
    DT(:,:,i) = cor_dcov(n, dm, phi, omg, kappa, i)
  end do

  do j = 1, 3
    do i = 1, j
      HT(:,:,i,j) = cor_hcov(n, dm, phi, omg, kappa, i, j)
      trUpsHT(i,j) = traceAB(Ups,HT(:,:,i,j),n)
      trUpsDTUpsDT(i,j) = traceABAC(Ups,DT(:,:,i),DT(:,:,j),n)
    end do
  end do

  do k = 1, Ntot
    call rchkusr
    where (itr .eq. 0)
      zsam = sample(:,k)
    elsewhere (itr .eq. 1)
      zsam = flink(sample(:,k),nu)
    elsewhere (itr .eq. 2)
      zsam = transfw(sample(:,k),nu)
    end where
    if (lmxi) zsam = zsam - xi
    call dsymv ('u',n,1d0,Ups,n,zsam,1,0d0,zU,1)
    zUz = dot_product(zsam,zU) + ssqdfsc
    do j = 1, 3
      zUDTUz(j) = qform(zU,DT(:,:,j),n)/zUz
    end do
    do j = 1, 3
      call dsymm ('r','u',n,n,1d0,DT(:,:,j),n,Ups,n,0d0,UpsDT,n)
      do i = 1, j
        zUHTUz = qform(zU,HT(:,:,i,j),n)
        call dsymm ('l','u',n,n,1d0,DT(:,:,i),n,UpsDT,n,0d0,DTUpsDT,n)
        DTUpsDT = .5d0*(DTUpsDT + transpose(DTUpsDT))
        zUDTUDTUz = qform(zU,DTUpsDT,n)
        hlglk(i,j,k) = -.5d0*trUpsHT(i,j) + .5d0*trUpsDTUpsDT(i,j) &
           + modeldfh*zUDTUz(i)*zUDTUz(j) &
           - modeldfh*(zUDTUDTUz + zUDTUDTUz - zUHTUz)/zUz
      end do
    end do
  end do

contains

  double precision function qform(v,A,n)
    ! Compute v'*A*v where A is symmetric with only the upper triangle
    ! specified.
    integer, intent(in) :: n
    double precision, intent(in) :: v(n), A(n,n)
    double precision Av(n)
    call dsymv ('u',n,1d0,A,n,v,1,0d0,Av,1)
    qform = dot_product(Av,v)
  end function qform

  double precision function traceAB(A,B,n)
    ! Compute the trace(matmul(A,B)) where A,B are symmetric matrices with
    ! only their upper triangular elements specified.
    integer, intent(in) :: n
    double precision, intent(in) :: A(n,n), B(n,n)
    integer j
    double precision c
    traceAB = A(1,1)*B(1,1)
    do j = 2, n
      c = dot_product(A(:j-1,j),B(:j-1,j))
      traceAB = traceAB + c + c + A(j,j)*B(j,j)
    end do
  end function traceAB

  double precision function traceABAC(A,B,C,n)
    ! Compute the trace(matmul(A,B,A,C)) where A is some matrix and B,C are
    ! symmetric matrices with only their upper triangular elements
    ! specified.
    integer, intent(in) :: n
    double precision, intent(in) :: A(n,n), B(n,n), C(n,n)
    !integer i,j
    double precision AB(n,n), AC(n,n)
    call dsymm('r','u',n,n,1d0,B,n,A,n,0d0,AB,n)
    call dsymm('r','u',n,n,1d0,C,n,A,n,0d0,AC,n)
    traceABAC = sum(AB*transpose(AC))
  end function traceABAC

  double precision function cor_dcov (n, dm, phi, omg, kappa, id) result (DT)
    ! id = which derivative do you want? phi, omg, kappa
    ! Return dT where T = omg*I + R(phi,kappa)
    use covfun, only: spcor_dh, spcor_dk, fill_symmetric_matrix, upper_tri
    implicit none
    integer, intent(in) :: n, id
    double precision, intent(in) :: dm(n,n), phi, omg, kappa
    !double precision, intent(out) :: DT(n,n)
    dimension :: DT(n,n)
    integer i
    double precision H(n,n)
    select case (id)
    case (1)
      if (phi .eq. 0d0) then
        DT = 0d0
        return
      else
        where (upper_tri())
          H = dm/phi
          DT = -(H/phi)*spcor_dh(H,kappa)
        end where
      end if
    case (2)
      DT = 0d0
      do i = 1, n
        DT(i,i) = 1d0
      end do
    case (3)
      if (phi .eq. 0d0) then
        DT = 0d0
        return
      else
        where (upper_tri())
          H = dm/phi
          DT = spcor_dk(H,kappa)
        end where
      end if
    end select
  end function cor_dcov

  double precision function cor_hcov (n, dm, phi, omg, kappa, id, jd) result (HT)
    ! id,jd = which derivative do you want? phi, omg, kappa
    ! Return d^2 T where T = omg*I + R(phi,kappa)
    use covfun, only: spcor_dh, spcor_dk, spcor_hh, spcor_hk, spcor_dhdk, &
       fill_symmetric_matrix, upper_tri
    implicit none
    integer, intent(in) :: n, id, jd
    double precision, intent(in) :: dm(n,n), phi, omg, kappa
    dimension :: HT(n,n)
    integer ij
    double precision H(n,n), H1(n,n)
    ij = id + 3*(jd-1)
    select case (ij)
    case (1) ! (phi,phi)
      if (phi .eq. 0d0) then
        HT = 0d0
        return
      else
        where (upper_tri())
          H = dm/phi
          H1 = -H/phi
          HT = H1*H1*spcor_hh(H,kappa) - (2/phi)*H1*spcor_dh(H,kappa)
        end where
      end if
    case (2,4) ! (phi,omg)
      HT = 0d0
    case (3,7) ! (phi,kappa)
      H = dm/phi
      H1 = -H/phi
      HT = H1*spcor_dhdk(H,kappa)
    case (5) ! (omg,omg)
      HT = 0d0
    case (6,8) ! (omg,kappa)
      HT = 0d0
    case (9) ! (kappa,kappa)
      if (phi .eq. 0d0) then
        HT = 0d0
        return
      else
        where (upper_tri())
          H = dm/phi
          HT = spcor_hk(H,kappa)
        end where
      end if
    end select
  end function cor_hcov
end subroutine llikfcn_hcov_tr



subroutine llikfcn_dlnk_tr (dlglk, phi, omg, nu, kappa, &
   sample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, ifam, itr)
!! Log-likelihood function derivative w.r.t link parameter.
!! itr is the type of transformation used: 0 = no, 1 = mu, 2 = wo.

  use modelfcns, logpdfydmu => logpdfydlnk
  use covfun
  use betaprior
  implicit none
  integer, intent(in) :: n, p, ifam, Ntot, icf, itr(n)
  double precision, intent(in) :: phi, omg, nu, &
     sample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappa
  double precision, intent(out) :: dlglk(Ntot)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, tsqval, respdfh, xi(n)
  integer j
  double precision zsam(n), msam(n), jsam(n), sam(n)
  double precision dlogpy(n), dmudnu(n), dy, dj, dz, zU(n), dzdnu(n)

  call create_model (ifam)
  call create_spcor (icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)
  call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)

  do j = 1, Ntot
    call rchkusr
    sam = sample(:,j)
    where (itr == 0)
      zsam = sam
      msam = invlink(zsam,nu)
      dlogpy = logpdfydmu(y,l,msam)
      dmudnu = invlinkdn(zsam,nu)
      dzdnu = 0d0
      jsam = 0d0
    elsewhere (itr == 1)
      msam = sam
      zsam = flink(msam,nu)
      dlogpy = 0d0
      dmudnu = 0d0
      dzdnu = -invlinkdn(zsam,nu)/invlinkdz(zsam,nu)
      jsam = (invlinkdzdn(zsam,nu) + invlinkhz(zsam,nu)*dzdnu) &
         /invlinkdz(zsam,nu)
      !dzdnu = invlinkdn(zsam,nu)
    elsewhere (itr == 2)
      zsam = transfw(sam,nu)
      msam = invlink(zsam,nu)
      dzdnu = -invtrwdn(zsam,nu)/invtrwdz(zsam,nu)
      dlogpy = logpdfydmu(y,l,msam)
      dmudnu = invlinkdn(zsam,nu) + invlinkdz(zsam,nu)*dzdnu
      jsam = (invtrwdzdn(zsam,nu) + invtrwhz(zsam,nu)*dzdnu)/invtrwdz(zsam,nu)
      !dzdnu = invtrwdn(zsam,nu)
    end where
    dy = dot_product(dlogpy,dmudnu)
    zU = logpdfz_dz(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh)
    dz = dot_product(zU,dzdnu)
    dj = -sum(jsam)
    dlglk(j) = dy + dz + dj
  end do
end subroutine llikfcn_dlnk_tr



subroutine llikfcn_hlnk_tr (hlglk, phi, omg, nu, kappa, &
   sample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, ifam, itr)
!! Log-likelihood function.
!! itr is the type of transformation used: 0 = no, 1 = mu, 2 = wo.

  use modelfcns, logpdfydmu => logpdfydlnk, logpdfyhmu => logpdfyhlnk
  use covfun
  use betaprior
  implicit none
  integer, intent(in) :: n, p, ifam, Ntot, icf, itr(n)
  double precision, intent(in) :: phi, omg, nu, &
     sample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappa
  double precision, intent(out) :: hlglk(Ntot)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, tsqval, respdfh, xi(n)
  integer j
  double precision zsam(n), msam(n), djsam(n), sam(n), hjsam(n)
  double precision dlogpy(n), dmudnu(n), dy, dj, dz, dlogpz(n), dzdnu(n)
  double precision hlogpy(n), hmuhnu(n), hy, hj, hz, hzhnu(n), hlogpz(n,n)

  call create_model (ifam)
  call create_spcor (icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
  do j = 1, Ntot
    call rchkusr
    sam = sample(:,j)
    where (itr == 0)
      zsam = sam
      msam = invlink(zsam,nu)
      dlogpy = logpdfydmu(y,l,msam)
      hlogpy = logpdfyhmu(y,l,msam)
      dmudnu = invlinkdn(zsam,nu)
      hmuhnu = invlinkhn(zsam,nu)
      djsam = 0d0
      dzdnu = 0d0
      hzhnu = 0d0
      hjsam = 0d0
    elsewhere (itr == 1)
      msam = sam
      zsam = flink(msam,nu)
      dlogpy = 0d0
      hlogpy = 0d0
      dzdnu = -invlinkdn(zsam,nu)/invlinkdz(zsam,nu)
      djsam = (invlinkdzdn(zsam,nu) + invlinkhz(zsam,nu)*dzdnu) &
         /invlinkdz(zsam,nu)
      hzhnu = -(invlinkhn(zsam,nu) + 2*invlinkdzdn(zsam,nu)*dzdnu &
         + invlinkhz(zsam,nu)*dzdnu*dzdnu)/invlinkdz(zsam,nu)
      hjsam = (invlink3z(zsam,nu)*dzdnu*dzdnu + invlinkhz(zsam,nu)*hzhnu &
         + 2*invlinkhzdn(zsam,nu)*dzdnu + invlinkdzhn(zsam,nu))&
         /invlinkdz(zsam,nu) - djsam*djsam
      dmudnu = 0d0
      hmuhnu = 0d0
    elsewhere (itr == 2)
      zsam = transfw(sam,nu)
      msam = invlink(zsam,nu)
      dlogpy = logpdfydmu(y,l,msam)
      hlogpy = logpdfyhmu(y,l,msam)
      dzdnu = -invtrwdn(zsam,nu)/invtrwdz(zsam,nu)
      hzhnu = -(invtrwhn(zsam,nu) + 2*invtrwdzdn(zsam,nu)*dzdnu &
         + invtrwhz(zsam,nu)*dzdnu*dzdnu)/invtrwdz(zsam,nu)
      dmudnu = invlinkdn(zsam,nu) + invlinkdz(zsam,nu)*dzdnu
      hmuhnu = invlinkhn(zsam,nu) + invlinkdz(zsam,nu)*hzhnu &
         + 2*invlinkdzdn(zsam,nu)*dzdnu + invlinkhz(zsam,nu)*dzdnu*dzdnu
      djsam = (invtrwdzdn(zsam,nu) + invtrwhz(zsam,nu)*dzdnu)/invtrwdz(zsam,nu)
      hjsam = (invtrw3z(zsam,nu)*dzdnu*dzdnu + invtrwhz(zsam,nu)*hzhnu &
         + 2*invtrwhzdn(zsam,nu)*dzdnu + invtrwdzhn(zsam,nu))&
         /invtrwdz(zsam,nu) - djsam*djsam
    end where
    dy = dot_product(dlogpy,dmudnu)
    dlogpz = logpdfz_dz(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh)
    hlogpz = logpdfz_hz(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh)
    dz = dot_product(dlogpz,dzdnu)
    dj = -sum(djsam)
    hy = dot_product(hlogpy,dmudnu*dmudnu) + dot_product(dlogpy,hmuhnu)
    hz = qform(dzdnu,hlogpz,n) + dot_product(dlogpz,hzhnu)
    hj = -sum(hjsam)
    hlglk(j) = hy + hz + hj
  end do

contains

  double precision function qform(v,A,n)
    ! Compute v'*A*v where A is symmetric with only the upper triangle
    ! specified.
    integer, intent(in) :: n
    double precision, intent(in) :: v(n), A(n,n)
    double precision Av(n)
    call dsymv ('u',n,1d0,A,n,v,1,0d0,Av,1)
    qform = dot_product(Av,v)
  end function qform
end subroutine llikfcn_hlnk_tr


subroutine llikfcn_dlnkdcov_tr (dlglk, phi, omg, nu, kappa, &
   sample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, ifam, itr)
!! Log-likelihood function derivative w.r.t link and covariance parameters.
!! itr is the type of transformation used: 0 = no, 1 = mu, 2 = wo.

  use modelfcns
  use covfun
  use betaprior
  implicit none
  integer, intent(in) :: n, p, ifam, Ntot, icf, itr(n)
  double precision, intent(in) :: phi, omg, nu, &
     sample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappa
  double precision, intent(out) :: dlglk(3,Ntot)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, tsqval, respdfh, xi(n)
  integer i, j
  double precision zsam(n), zU(n), zUz, zUDTUz, DT(n,n,3), Dz(n), DzU(n), &
     DzUDTUz, DzUz, DTzU(n)

  call create_model (ifam)
  call create_spcor (icf,n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifam)
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)

  do i = 1, 3
    DT(:,:,i) = cor_dcov(n, dm, phi, omg, kappa, i)
  end do

  do j = 1, Ntot
    call rchkusr
    where (itr .eq. 0)
      zsam = sample(:,j)
      Dz = 0d0
    elsewhere (itr .eq. 1)
      zsam = flink(sample(:,j),nu)
      Dz = -invlinkdn(zsam,nu)/invlinkdz(zsam,nu)
    elsewhere (itr .eq. 2)
      zsam = transfw(sample(:,j),nu)
      Dz = -invtrwdn(zsam,nu)/invtrwdz(zsam,nu)
    end where
    if (lmxi) zsam = zsam - xi
    call dsymv ('u',n,1d0,Ups,n,zsam,1,0d0,zU,1)
    zUz = dot_product(zsam,zU) + ssqdfsc
    call dsymv ('u',n,1d0,Ups,n,Dz,1,0d0,DzU,1)
    DzUz = dot_product(Dz,zU)
    do i = 1, 3
      call dsymv ('u',n,1d0,DT(:,:,i),n,zU,1,0d0,DTzU,1)
      zUDTUz = dot_product(zU,DTzU)
      DzUDTUz = dot_product(DzU,DTzU)
      dlglk(i,j) = (2*modeldfh)*(DzUDTUz/zUz - zUDTUz*DzUz/(zUz*zUz))
    end do
  end do

contains

  double precision function qform(v,A,n)
    ! Compute v'*A*v where A is symmetric with only the upper triangle
    ! specified.
    integer, intent(in) :: n
    double precision, intent(in) :: v(n), A(n,n)
    double precision Av(n)
    call dsymv ('u',n,1d0,A,n,v,1,0d0,Av,1)
    qform = dot_product(Av,v)
  end function qform

  double precision function cor_dcov (n, dm, phi, omg, kappa, id) result (DT)
    ! id = which derivative do you want? phi, omg, kappa
    ! Return dT where T = omg*I + R(phi,kappa)
    use covfun, only: spcor_dh, spcor_dk, fill_symmetric_matrix, upper_tri
    implicit none
    integer, intent(in) :: n, id
    double precision, intent(in) :: dm(n,n), phi, omg, kappa
    !double precision, intent(out) :: DT(n,n)
    dimension :: DT(n,n)
    integer i
    double precision H(n,n)
    select case (id)
    case (1)
      if (phi .eq. 0d0) then
        DT = 0d0
        return
      else
        where (upper_tri())
          H = dm/phi
          DT = -(H/phi)*spcor_dh(H,kappa)
        end where
      end if
    case (2)
      DT = 0d0
      do i = 1, n
        DT(i,i) = 1d0
      end do
    case (3)
      if (phi .eq. 0d0) then
        DT = 0d0
        return
      else
        where (upper_tri())
          H = dm/phi
          DT = spcor_dk(H,kappa)
        end where
      end if
    end select
  end function cor_dcov
end subroutine llikfcn_dlnkdcov_tr

subroutine llikfcnmc_11 (lglk, nulist, kappalist, &
   sample, phisample, omgsample, phipars, omgpars, Ntot, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icflist, n, p, kg, ifamlist, itr)
!! Log-likelihood function.
!! itr is the type of transformation used: 0 = no, 1 = mu, 2 = wo.

  use modelfcns, logpdfzf => logpdfz, condymu_mf => condymu
  use covfun
  use betaprior
  use condymu, only: condymu_gt
  implicit none
  integer, intent(in) :: n, p, kg, ifamlist(kg), Ntot, icflist(kg), itr(n)
  double precision, intent(in) :: nulist(kg), phipars(4,kg), omgpars(4,kg), &
     sample(n,Ntot,kg), phisample(Ntot), omgsample(Ntot), &
     y(n), l(n), F(n,p), &
     dm(n,n), betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot,kg)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, nu, phi, omg, kappa, &
     tsqval, respdfh, xi(n)
  integer i, j, ifam, icf
  double precision zsam(n), msam(n), jsam(n), sam(n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifamlist(1))
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  do i = 1, kg
    ifam = ifamlist(i)
    icf = icflist(i)
    call create_model (ifam)
    call create_spcor (icf,n)
    nu = nulist(i)
    kappa = kappalist(i)
    do j = 1, Ntot
      call rchkusr
      phi = phisample(j)
      omg = omgsample(j)
      call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
      sam = sample(:,j,i)
      where (itr == 0)
        zsam = sam
        msam = invlink(zsam,nu)
        jsam = 0d0
      elsewhere (itr == 1)
        msam = sam
        zsam = flink(msam,nu)
        jsam = loginvlinkdz(zsam,nu)
      elsewhere (itr == 2)
        zsam = transfw(sam,nu)
        msam = invlink(zsam,nu)
        jsam = loginvtrwdz(zsam,nu)
      end where
      lglk(j,i) = logpdfzf(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh) &
         + condymuf(ifam,n,y,l,msam,tsqval,respdfh) - sum(jsam) &
         + lgengam(phi,phipars(:,i)) + lgengam(omg,omgpars(:,i))
    end do
  end do

contains
  pure function lgengam (x, q)
    double precision, intent(in) :: x, q(4)
    double precision lgengam
    lgengam = (x - q(4))/q(1)
    lgengam = log(lgengam)
    lgengam = (q(2)-1)*lgengam - exp(lgengam*q(3))
  end function lgengam

  pure function condymu_sp (n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymu_sp
    condymu_sp = condymu_mf(n,y1,y2,mu,tsqval)
  end function condymu_sp

  pure function condymuf (ifam, n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: ifam
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymuf
    select case (ifam)
    case (0)
      condymuf = condymu_gt(n, y1, y2, mu, tsqval, respdfh)
    case default
      condymuf = condymu_sp(n, y1, y2, mu, tsqval, respdfh)
    end select
  end function condymuf
end subroutine llikfcnmc_11

subroutine llikfcnmc_01 (lglk, nulist, philist, kappalist, &
   sample, omgsample, omgpars, Ntot, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icflist, n, p, kg, ifamlist, itr)
!! Log-likelihood function.
!! itr is the type of transformation used: 0 = no, 1 = mu, 2 = wo.

  use modelfcns, logpdfzf => logpdfz, condymu_mf => condymu
  use covfun
  use betaprior
  use condymu, only: condymu_gt
  implicit none
  integer, intent(in) :: n, p, kg, ifamlist(kg), Ntot, icflist(kg), itr(n)
  double precision, intent(in) :: nulist(kg), omgpars(4,kg), &
     sample(n,Ntot,kg), philist(kg), omgsample(Ntot), y(n), l(n), F(n,p), &
     dm(n,n), betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot,kg)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, nu, phi, omg, kappa, &
     tsqval, respdfh, xi(n)
  integer i, j, ifam, icf
  double precision zsam(n), msam(n), jsam(n), sam(n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifamlist(1))
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  do i = 1, kg
    ifam = ifamlist(i)
    icf = icflist(i)
    call create_model (ifam)
    call create_spcor(icf,n)
    nu = nulist(i)
    kappa = kappalist(i)
    phi = philist(i)
    do j = 1, Ntot
      call rchkusr
      omg = omgsample(j)
      call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
      sam = sample(:,j,i)
      where (itr == 0)
        zsam = sam
        msam = invlink(zsam,nu)
        jsam = 0d0
      elsewhere (itr == 1)
        msam = sam
        zsam = flink(msam,nu)
        jsam = loginvlinkdz(zsam,nu)
      elsewhere (itr == 2)
        zsam = transfw(sam,nu)
        msam = invlink(zsam,nu)
        jsam = loginvtrwdz(zsam,nu)
      end where
      lglk(j,i) = logpdfzf(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh) &
         + condymuf(ifam,n,y,l,msam,tsqval,respdfh) - sum(jsam) &
         + lgengam(omg,omgpars(:,i))
    end do
  end do

contains
  pure function lgengam (x, q)
    double precision, intent(in) :: x, q(4)
    double precision lgengam
    lgengam = (x - q(4))/q(1)
    lgengam = log(lgengam)
    lgengam = (q(2)-1)*lgengam - exp(lgengam*q(3))
  end function lgengam

  pure function condymu_sp (n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymu_sp
    condymu_sp = condymu_mf(n,y1,y2,mu,tsqval)
  end function condymu_sp

  pure function condymuf (ifam, n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: ifam
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymuf
    select case (ifam)
    case (0)
      condymuf = condymu_gt(n, y1, y2, mu, tsqval, respdfh)
    case default
      condymuf = condymu_sp(n, y1, y2, mu, tsqval, respdfh)
    end select
  end function condymuf
end subroutine llikfcnmc_01

subroutine llikfcnmc_10 (lglk, nulist, omglist, kappalist, &
   sample, phisample, phipars, Ntot, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icflist, n, p, kg, ifamlist, itr)
!! Log-likelihood function.
!! itr is the type of transformation used: 0 = no, 1 = mu, 2 = wo.

  use modelfcns, logpdfzf => logpdfz, condymu_mf => condymu
  use covfun
  use betaprior
  use condymu, only: condymu_gt
  implicit none
  integer, intent(in) :: n, p, kg, ifamlist(kg), Ntot, icflist(kg), itr(n)
  double precision, intent(in) :: nulist(kg), phipars(4,kg), &
     sample(n,Ntot,kg), omglist(kg), phisample(Ntot), y(n), l(n), F(n,p), &
     dm(n,n), betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot,kg)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, nu, phi, omg, kappa, &
     tsqval, respdfh, xi(n)
  integer i, j, ifam, icf
  double precision zsam(n), msam(n), jsam(n), sam(n)

  ssqdfsc = ssqdf*ssqsc
  select case (ifamlist(1))
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  do i = 1, kg
    ifam = ifamlist(i)
    icf = icflist(i)
    call create_model (ifam)
    call create_spcor(icf,n)
    nu = nulist(i)
    kappa = kappalist(i)
    omg = omglist(i)
    do j = 1, Ntot
      call rchkusr
      phi = phisample(j)
      call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
      sam = sample(:,j,i)
      where (itr == 0)
        zsam = sam
        msam = invlink(zsam,nu)
        jsam = 0d0
      elsewhere (itr == 1)
        msam = sam
        zsam = flink(msam,nu)
        jsam = loginvlinkdz(zsam,nu)
      elsewhere (itr == 2)
        zsam = transfw(sam,nu)
        msam = invlink(zsam,nu)
        jsam = loginvtrwdz(zsam,nu)
      end where
      lglk(j,i) = logpdfzf(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh) &
         + condymuf(ifam,n,y,l,msam,tsqval,respdfh) - sum(jsam) &
         + lgengam(phi,phipars(:,i))
    end do
  end do

contains
  pure function lgengam (x, q)
    double precision, intent(in) :: x, q(4)
    double precision lgengam
    lgengam = (x - q(4))/q(1)
    lgengam = log(lgengam)
    lgengam = (q(2)-1)*lgengam - exp(lgengam*q(3))
  end function lgengam

  pure function condymu_sp (n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymu_sp
    condymu_sp = condymu_mf(n,y1,y2,mu,tsqval)
  end function condymu_sp

  pure function condymuf (ifam, n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: ifam
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymuf
    select case (ifam)
    case (0)
      condymuf = condymu_gt(n, y1, y2, mu, tsqval, respdfh)
    case default
      condymuf = condymu_sp(n, y1, y2, mu, tsqval, respdfh)
    end select
  end function condymuf
end subroutine llikfcnmc_10

subroutine llikfcnmc_00 (lglk, nulist, philist, omglist, kappalist, &
   sample, Ntot, &
   y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icflist, n, p, kg, ifamlist, itr)
!! Log-likelihood function.
!! itr is the type of transformation used: 0 = no, 1 = mu, 2 = wo.

  use modelfcns, logpdfzf => logpdfz, condymu_mf => condymu
  use covfun
  use betaprior
  use condymu, only: condymu_gt
  implicit none
  integer, intent(in) :: n, p, kg, ifamlist(kg), Ntot, icflist(kg), itr(n)
  double precision, intent(in) :: nulist(kg), philist(kg), &
     sample(n,Ntot,kg), omglist(kg), y(n), l(n), F(n,p), &
     dm(n,n), betm0(p), betQ0(p,p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot,kg)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n,n), &
     ldh_Ups, ssqdfsc, modeldfh, nu, phi, omg, kappa, &
     tsqval, respdfh, xi(n)
  integer i, j, ifam, icf
  double precision zsam(n), msam(n), jsam(n), sam(n)
  double precision ll1, ll2, ll3

  ssqdfsc = ssqdf*ssqsc
  select case (ifamlist(1))
  case (0)
    tsqval = tsqdf*tsq
    respdfh = .5d0*(n + tsqdf)
  case default
    tsqval = tsq
  end select

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  call rchkusr

  do i = 1, kg
    ifam = ifamlist(i)
    icf = icflist(i)
    call create_model (ifam)
    call create_spcor(icf,n)
    nu = nulist(i)
    kappa = kappalist(i)
    phi = philist(i)
    omg = omglist(i)
    call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
    do j = 1, Ntot
      call rchkusr
      sam = sample(:,j,i)
      where (itr == 0)
        zsam = sam
        msam = invlink(zsam,nu)
        jsam = 0d0
      elsewhere (itr == 1)
        msam = sam
        zsam = flink(msam,nu)
        jsam = loginvlinkdz(zsam,nu)
      elsewhere (itr == 2)
        zsam = transfw(sam,nu)
        msam = invlink(zsam,nu)
        jsam = loginvtrwdz(zsam,nu)
      end where
      ll1 = logpdfzf(n,zsam,Ups,ldh_Ups,xi,lmxi,ssqdfsc,modeldfh)
      ll2 = condymuf(ifam,n,y,l,msam,tsqval,respdfh)
      ll3 = sum(jsam)
      lglk(j,i) = ll1 + ll2 - ll3
    end do
  end do
  
contains
  pure function condymu_sp (n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymu_sp
    condymu_sp = condymu_mf(n,y1,y2,mu,tsqval)
  end function condymu_sp

  pure function condymuf (ifam, n, y1, y2, mu, tsqval, respdfh)
    implicit none
    integer, intent(in) :: ifam
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymuf
    select case (ifam)
    case (0)
      condymuf = condymu_gt(n, y1, y2, mu, tsqval, respdfh)
    case default
      condymuf = condymu_sp(n, y1, y2, mu, tsqval, respdfh)
    end select
  end function condymuf
end subroutine llikfcnmc_00

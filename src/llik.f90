subroutine llikfcn_no (lglk, philist, omglist, nulist, kappalist, &
   zsample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, kg, ifam, itr)

  use modelfcns, jointyz_sp => jointyz
  use covfun
  use jointyz, only: jointyz_gt
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, Ntot, icf, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     zsample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot, kg)
  logical lmxi
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, &
     tsqdfsc, respdfh, xi(n)
  integer i, j

  call create_model (ifam)
  call create_spcor (icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  select case (ifam)
  case (0)
    do i = 1, kg
      call calc_cov (philist(i),omglist(i),dm,F,betQ0,&
         kappalist(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,i) = jointyz_gt(n, zsample(:,j), y, l, Ups, ldh_Ups, &
           nulist(i), xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
      end do
    end do
  case default
    do i = 1, kg
      call calc_cov (philist(i),omglist(i),dm,F,betQ0,&
         kappalist(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,i) = jointyz_sp(n, zsample(:,j), y, l, Ups, ldh_Ups, &
           nulist(i), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  end select
end subroutine llikfcn_no




subroutine llikfcn_mu (lglk, philist, omglist, nulist, kappalist, &
   musample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, kg, ifam, itr)

  use modelfcns, jointymu_sp => jointymu, logpdfmu_sp => logpdfmu
  use covfun
  use jointymu, only: jointymu_gt
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, Ntot, icf, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     musample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot, kg)
  logical lmxi
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, &
     tsqdfsc, respdfh, xi(n)
  integer i, j

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  select case (ifam)
  case (0)
    do i = 1, kg
      call calc_cov (philist(i),omglist(i),dm,F,betQ0,&
         kappalist(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,i) = jointymu_gt(n, musample(:,j), y, l, Ups, ldh_Ups, &
           nulist(i), xi, lmxi, ssqdfsc, tsqdfsc, modeldfh, respdfh)
      end do
    end do
  case default
    do i = 1, kg
      call calc_cov (philist(i),omglist(i),dm,F,betQ0,&
         kappalist(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        lglk(j,i) = jointymu_sp(n, musample(:,j), y, l, Ups, ldh_Ups, &
           nulist(i), xi, lmxi, ssqdfsc, tsq, modeldfh)
      end do
    end do
  end select
end subroutine llikfcn_mu

subroutine lpdffcn_mu (lglk, philist, omglist, nulist, kappalist, &
   musample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, kg, ifam, itr)
  ! Like llikfcn_mu but not multiplying by p(y|mu)

  use modelfcns, jointymu_sp => jointymu, logpdfmu_sp => logpdfmu
  use covfun
  use jointymu, only: jointymu_gt
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, Ntot, icf, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     musample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot, kg)
  logical lmxi
  double precision T(n, n), TiF(n, p), FTF(p, p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, &
     tsqdfsc, respdfh, xi(n)
  integer i, j

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  do i = 1, kg
    call calc_cov (philist(i),omglist(i),dm,F,betQ0,&
       kappalist(i),n,p,T,TiF,FTF,Ups,ldh_Ups)
    do j = 1, Ntot
      call rchkusr
      lglk(j,i) = logpdfmu_sp(n, musample(:,j), Ups, ldh_Ups, &
         nulist(i), xi, lmxi, ssqdfsc, modeldfh)
    end do
  end do
end subroutine lpdffcn_mu


subroutine llikfcn_wo (lglk, philist, omglist, nulist, kappalist, &
   wsample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, kg, ifam, itr)
  !! Work-around for the binomial family.

  use modelfcns, jointyz_sp => jointyz
  use covfun
  use betaprior
  implicit none
  integer, intent(in) :: n, p, kg, ifam, Ntot, icf, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     wsample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot, kg)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, nu, phi, omg, kappa, &
     tsqdfsc, respdfh, xi(n)
  integer i, j, m
  double precision zsam(n)

  call create_model (ifam)
  call create_spcor(icf,n)

  ssqdfsc = ssqdf*ssqsc
  tsqdfsc = tsqdf*tsq
  respdfh = .5d0*(n + tsqdf)

  ! Determine flat or normal prior
  call betapriorz (modeldfh, xi, lmxi, betm0, betQ0, F, n, p, ssqdf)

  select case (ifam)
  case (0)
    call rexit ("This method has not been implemented.")
  case default
    do i = 1, kg
      nu = nulist(i)
      phi = philist(i)
      omg = omglist(i)
      kappa = kappalist(i)
      call calc_cov (phi,omg,dm,F,betQ0,&
         kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
      do j = 1, Ntot
        call rchkusr
        zsam = transfw(wsample(:,j),nu)
        lglk(j,i) = jointyz_sp(n, zsam, y, l, Ups, ldh_Ups, &
           nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
        do m = 1, n
          lglk(j,i) = lglk(j,i) - loginvtrwdz(zsam(m),nu)
        end do
      end do
    end do
  end select
end subroutine llikfcn_wo

subroutine llikfcn_tr (lglk, philist, omglist, nulist, kappalist, &
   sample, Ntot, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsq, &
   icf, n, p, kg, ifam, itr)
!! Log-likelihood function.
!! itr is the type of transformation used: 0 = no, 1 = mu, 2 = wo.

  use modelfcns, logpdfzf => logpdfz, condymu_mf => condymu
  use covfun
  use betaprior
  use condymu, only: condymu_gt
  implicit none
  integer, intent(in) :: n, p, kg, ifam, Ntot, icf, itr(n)
  double precision, intent(in) :: philist(kg), omglist(kg), nulist(kg), &
     sample(n, Ntot), y(n), l(n), F(n, p), &
     dm(n, n), betm0(p), betQ0(p, p), ssqdf, ssqsc, tsqdf, tsq, kappalist(kg)
  double precision, intent(out) :: lglk(Ntot,kg)
  logical lmxi
  double precision T(n,n), TiF(n,p), FTF(p,p), Ups(n, n), &
     ldh_Ups, ssqdfsc, modeldfh, nu, phi, omg, kappa, &
     tsqval, respdfh, xi(n)
  integer i, j
  double precision zsam(n), msam(n), jsam(n), sam(n)

  call create_model (ifam)
  call create_spcor(icf,n)

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

  do i = 1, kg
    nu = nulist(i)
    phi = philist(i)
    omg = omglist(i)
    kappa = kappalist(i)
    call calc_cov (phi,omg,dm,F,betQ0,kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
    do j = 1, Ntot
      call rchkusr
      sam = sample(:,j)
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
         + condymuf(ifam,n,y,l,msam,tsqval,respdfh) - sum(jsam)
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
end subroutine llikfcn_tr



! Local Variables:
! compile-command: "gfortran -c -fpic -Wunused-parameter -Wall \
!   -pedantic -o llik.o llik.f90"
! End:

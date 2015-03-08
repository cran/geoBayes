!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module covfun
  implicit none 
  public
  private covmat_a, covmat_l

  interface covmat
    module procedure covmat_a, covmat_l
  end interface covmat

contains
  subroutine calc_cov (phi,nsq,dm,F,betQ0, &
     lup,kappa,icf,n,p,T,TiF,FTF,TFFT,Ups,ldh_Ups)
    implicit none
    logical, intent(in) :: lup(n,n)
    integer, intent(in) :: n, p, icf
    double precision, intent(in) :: phi, nsq, dm(n,n), &
       F(n,p), betQ0(p,p), kappa
    double precision, intent(out) :: T(n,n), FTF(p,p), TFFT(n,p), Ups(n,n), &
       ldh_Ups, TiF(n,p)
    integer i, j
    double precision nsqp1, Tih(n,n), ldh_T, ldh_FTF

    nsqp1 = nsq + 1d0
    T = dm
    call covmat_l(T,phi,kappa,n,n,lup,icf)
    do i = 1, n
      T(i,i) = nsqp1
    end do
    Tih = T
    ldh_T = oppdf(Tih,n) ! Tih is uppertri s.t. Tih*Tih' = T^{-1}
    TiF = F
    call dtrmm ('l','u','t','n',n,p,1d0,Tih,n,TiF,n) ! TiF = T^{-1/2}*F
    FTF = betQ0
    call dsyrk ('u','t',p,n,1d0,TiF,n,1d0,FTF,p) ! FTF = F'*T^{-1}*F + Q0
    call dtrmm ('l','u','n','n',n,p,1d0,Tih,n,TiF,n) ! TiF = T^{-1}*F
    ldh_FTF = oppdf(FTF,p) ! FTF is uptri s.t. FTF*FTF' = (F'*T^{-1}*F + Q0)^{-1}
    TFFT = TiF
    call dtrmm ('r','u','n','n',n,p,1d0,FTF,p,TFFT,n) ! TFFT =
    ! T^{-1}*F*(F'*T^{-1}*F + Q0)^{-1/2}
    do i = 1, n
      Ups(1:i,i) = 0d0
      do j = i, 1, -1
        Ups(1:j,j) = Ups(1:j,j) + Tih(1:j,i)*Tih(j,i) ! Ups = T^{-1}
      end do
    end do
    T = Ups ! T is now T^{-1}
    call dsyrk ('u','n',n,p,-1d0,TFFT,n,1d0,Ups,n) ! Ups = T^{-1} -
    ! T^{-1}*F*(F'*T^{-1}*F  Q0)^{-1}*F'*T^{-1}
    ldh_Ups = -ldh_T - ldh_FTF
  end subroutine calc_cov

  subroutine calc_cov_pred (z0_ups, TC, FCTF, &
     phi, nsq, dmdm0, F, F0, kappa, icf, T, n, n0, p)
    implicit none
    integer, intent(in) :: n, n0, p, icf
    double precision, intent(in) :: phi, nsq, dmdm0(n,n0), F(n,p), &
       F0(n0,p), kappa, T(n,n)
    double precision, intent(out) :: z0_ups(n0), TC(n,n0), FCTF(n0,p)
    double precision C(n,n0), nsqp1

    nsqp1 = nsq + 1d0
    C = dmdm0
    call covmat_a (C,phi,kappa,n,n0,icf)
    call dsymm ('l','u',n,n0,1d0,T,n,C,n,0d0,TC,n)
    z0_ups = sqrt(nsqp1 - sum(TC*C,1))
    FCTF = F0
    call dgemm ('t','n',n0,p,n,-1d0,TC,n,F,n,1d0,FCTF,n0)
  end subroutine calc_cov_pred


  subroutine covlist (kg,philist,nsqlist,n,p,betQ0,F,dm,kappa,icf, &
     Ulist,ldh_Ulist)
!!! Compute list of relevant covariance matrices

    implicit none

    integer, intent(in) :: kg, n, p, icf
    double precision, intent(in) :: philist(kg), betQ0(p,p), nsqlist(kg), &
       F(n,p), dm(n,n), kappa
    double precision, intent(out) :: ldh_Ulist(kg), Ulist(n,n,kg)
    integer i
    double precision FTF(p,p), TFFT(n,p), T(n,n), TiF(n,p)
    logical lup(n,n)

    do i = 1, n
      lup(:i-1,i) = .true.
      lup(i:,i) = .false. 
    end do

    do i = 1, kg
      call calc_cov(philist(i),nsqlist(i),dm,F,betQ0,&
         lup,kappa,icf,n,p,T,TiF,FTF,TFFT,Ulist(:,:,i),ldh_Ulist(i))
    end do
  end subroutine covlist


  subroutine covmat_a (dm,phi,kappa,n1,n2,icf)
    use interfaces, only: matern
    implicit none
    integer, intent(in) :: n1, n2, icf
    double precision, intent(in) :: phi, kappa
    double precision, intent(inout) :: dm(n1,n2)
    if (phi .eq. 0d0) then
      where (dm .eq. 0d0)
        dm = 1d0
      elsewhere
        dm = 0d0
      end where
    else if (phi .gt. 0d0) then
      select case (icf)
      case (1) ! Matern
        if (kappa .eq. .5d0) then
          dm = dm/phi
          dm=exp(-dm)
        else if (kappa .eq. 1.5d0) then
          dm = dm/phi
          dm=exp(-dm)*(1d0+dm)
        else if (kappa .eq. 2.5d0) then
          dm = dm/phi
          dm=exp(-dm)*(1d0 + dm + dm*dm/3d0)
        else if (kappa .gt. 0d0) then
          dm = dm/phi
          dm = matern(dm,kappa)
        else
          call rexit ('covmat - Invalid kappa')
        end if
      case (2) ! Spherical
        where (dm .lt. phi)
          dm = dm/phi
          dm = 1d0 - 1.5d0*dm + .5d0*dm*dm*dm
        elsewhere
          dm = 0d0
        end where
      case (3) ! Power-exponential
        if (kappa .le. 2d0 .and. kappa .ge. 0d0) then
          if (kappa .eq. 1d0) then ! exponential
            dm = dm/phi
            dm = exp(-dm)
          else if (kappa .eq. 2d0) then ! Gaussian
            dm = dm/phi
            dm = dm*dm
            dm = exp(-dm)
          else
            dm = dm/phi
            dm = dm**(kappa)
            dm = exp(-dm)
          end if
        else
          call rexit ('covmat - Invalid kappa')
        end if
      end select
    else
      call rexit ('covmat - Negative phi')
    end if
  end subroutine covmat_a


  subroutine covmat_l (dm,phi,kappa,n1,n2,ldm,icf)
    use interfaces, only: matern
    implicit none
    integer, intent(in) :: n1, n2, icf
    double precision, intent(in) :: phi, kappa
    double precision, intent(inout) :: dm(n1,n2)
    logical, intent(in) :: ldm(n1,n2)
    if (phi .eq. 0d0) then
      where (ldm .and. (dm .eq. 0d0))
        dm = 1d0
      elsewhere (ldm)
        dm = 0d0
      end where
    else if (phi .gt. 0d0) then
      select case (icf)
      case (1) ! Matern
        if (kappa .eq. .5d0) then
          where (ldm)
            dm = dm/phi
            dm=exp(-dm)
          end where
        else if (kappa .eq. 1.5d0) then
          where (ldm)
            dm = dm/phi
            dm=exp(-dm)*(1d0+dm)
          end where
        else if (kappa .eq. 2.5d0) then
          where (ldm)
            dm = dm/phi
            dm=exp(-dm)*(1d0 + dm + dm*dm/3d0)
          end where
        else if (kappa .gt. 0d0) then
          where (ldm)
            dm = dm/phi
            dm = matern(dm,kappa)
          end where
        else
          call rexit('covmat - Invalid kappa')
        end if
      case (2) ! Spherical
        where (ldm .and. dm .lt. phi)
          dm = dm/phi
          dm = 1d0 - 1.5d0*dm + .5d0*dm*dm*dm
        elsewhere (ldm)
          dm = 0d0
        end where
      case (3)
        if (kappa .le. 2d0 .and. kappa .ge. 0d0) then
          if (kappa .eq. 1d0) then ! exponential
            where (ldm)
              dm = dm/phi
              dm = exp(-dm)
            end where
          else if (kappa .eq. 2d0) then
            where (ldm)
              dm = dm/phi
              dm = dm*dm
              dm = exp(-dm)
            end where
          else
            where (ldm)
              dm = dm/phi
              dm = dm**(kappa)
              dm = exp(-dm)
            end where
          end if
        else
          call rexit ('covmat - Invalid kappa')
        end if
      end select
    else
      call rexit ('covmat - Negative phi')
    end if
  end subroutine covmat_l

  function oppdf (A,n) result (ldh)
    ! INPUT:
    ! A - A positive definite matrix
    ! OUTPUT:
    ! A - An upper triangular matrix U s.t. U*U' = A^{-1}
    ! ldh - .5* log det (A) where A is the original (input) matrix
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: A(n,n)
    double precision ldh
    integer i
    ! First call cholesky decomposition
    call dpotrf ('u',n,A,n,i)
    if (i /= 0) then
      call rexit ("Error: oppdf - Matrix not positive definite")
    end if
    ldh = 0d0
    do i = 1, n
      ldh = ldh + log(A(i,i))
    end do
    ! Next invert
    call dtrtri ('u','n',n,A,n,i)
    if (i /= 0) then
      call rexit ("Error: oppdf - Matrix not invertible")
    end if
  end function oppdf
end module covfun


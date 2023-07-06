!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Commentary:
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module covfun
  use cor_fcns
  implicit none
  public
  private :: covmat_a, covmat_l, covmat_v
  public :: spcor, spcor_dh, spcor_dk, spcor_hh, spcor_hk, &
     spcor_dhdk, upper_tri
  integer, parameter, private :: CORCODES(5) = (/1,2,3,4,5/)
  logical, private :: CORDEF = .false.
  integer, private :: CORIS = 0, DIMSP = 0
  logical, private, pointer :: lup(:,:) => null()


  interface covmat
    module procedure covmat_a, covmat_l, covmat_v
  end interface covmat

contains

  subroutine create_spcor (icf,n)
    use cor_fcns
    implicit none
    integer, intent(in) :: icf
    integer, intent(in) :: n
    integer i

    if (.not.(CORDEF .and. (CORIS .eq. icf))) then
      if (.not. (any(icf .eq. CORCODES))) then
        call rexit ("Unrecognised correlation.")
      end if
      CORIS = icf
      CORDEF = .true.
    end if

    if (n .gt. 0) then
      if (DIMSP .ne. n) then
        nullify(lup)
        allocate(lup(n,n))
        DIMSP = n
        do i = 1, n
          lup(1:i,i) = .true.
          lup(i+1:n,i) = .false.
        end do
      end if
    end if
  end subroutine create_spcor

  function upper_tri () result (lup_)
    logical, pointer :: lup_(:,:)
    allocate(lup_(DIMSP,DIMSP))
    lup_ = lup
  end function upper_tri

  pure elemental double precision function spcor (h,k)
    implicit none
    double precision, intent(in) :: h,k
    select case (CORIS)
    case (1) ! Matern
      spcor = cor_matern(h,k)
    case (2) ! Spherical
      spcor = cor_spher(h,k)
    case (3) ! Power exponential
      spcor = cor_powexp(h,k)
    case (4) ! Exponential
      spcor = cor_exp(h,k)
    case (5) ! Gaussian
      spcor = cor_gaussian(h,k)
    end select
  end function spcor

  pure elemental double precision function spcor_dh (h,k)
    implicit none
    double precision, intent(in) :: h,k
    select case (CORIS)
    case (1) ! Matern
      spcor_dh = cor_dh_matern(h,k)
    case (2) ! Spherical
      spcor_dh = cor_dh_spher(h,k)
    case (3) ! Power exponential
      spcor_dh = cor_dh_powexp(h,k)
    case (4) ! Exponential
      spcor_dh = cor_dh_exp(h,k)
    case (5) ! Gaussian
      spcor_dh = cor_dh_gaussian(h,k)
    end select
  end function spcor_dh

  pure elemental double precision function spcor_dk (h,k)
    implicit none
    double precision, intent(in) :: h,k
    select case (CORIS)
    case (1) ! Matern
      spcor_dk = cor_dk_matern(h,k)
    case (2) ! Spherical
      spcor_dk = cor_dk_spher(h,k)
    case (3) ! Power exponential
      spcor_dk = cor_dk_powexp(h,k)
    case (4) ! Exponential
      spcor_dk = cor_dk_exp(h,k)
    case (5) ! Gaussian
      spcor_dk = cor_dk_gaussian(h,k)
    end select
  end function spcor_dk

  pure elemental double precision function spcor_hh (h,k)
    implicit none
    double precision, intent(in) :: h,k
    select case (CORIS)
    case (1) ! Matern
      spcor_hh = cor_hh_matern(h,k)
    case (2) ! Spherical
      spcor_hh = cor_hh_spher(h,k)
    case (3) ! Power exponential
      spcor_hh = cor_hh_powexp(h,k)
    case (4) ! Exponential
      spcor_hh = cor_hh_exp(h,k)
    case (5) ! Gaussian
      spcor_hh = cor_hh_gaussian(h,k)
    end select
  end function spcor_hh

  pure elemental double precision function spcor_hk (h,k)
    implicit none
    double precision, intent(in) :: h,k
    select case (CORIS)
    case (1) ! Matern
      spcor_hk = cor_hk_matern(h,k)
    case (2) ! Spherical
      spcor_hk = cor_hk_spher(h,k)
    case (3) ! Power exponential
      spcor_hk = cor_hk_powexp(h,k)
    case (4) ! Exponential
      spcor_hk = cor_hk_exp(h,k)
    case (5) ! Gaussian
      spcor_hk = cor_hk_gaussian(h,k)
    end select
  end function spcor_hk

  pure elemental double precision function spcor_dhdk (h,k)
    implicit none
    double precision, intent(in) :: h,k
    select case (CORIS)
    case (1) ! Matern
      spcor_dhdk = cor_dhdk_matern(h,k)
    case (2) ! Spherical
      spcor_dhdk = cor_dhdk_spher(h,k)
    case (3) ! Power exponential
      spcor_dhdk = cor_dhdk_powexp(h,k)
    case (4) ! Exponential
      spcor_dhdk = cor_dhdk_exp(h,k)
    case (5) ! Gaussian
      spcor_dhdk = cor_dhdk_gaussian(h,k)
    end select
  end function spcor_dhdk



  subroutine calc_cov (phi,omg,dm,F,betQ0, &
     kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
    implicit none
    integer, intent(in) :: n, p
    double precision, intent(in) :: phi, omg, dm(n,n), &
       F(n,p), betQ0(p,p), kappa
    double precision, intent(out) :: T(n,n), FTF(p,p), Ups(n,n), &
       ldh_Ups, TiF(n,p)
    integer i, j
    double precision omgp1, Tih(n,n), TFFT(n,p), ldh_T, ldh_FTF

    omgp1 = omg + 1d0
    T = dm
    call covmat_l(T,phi,kappa,n,n,lup)
    do i = 1, n
      T(i,i) = omgp1
    end do
    Tih = T
    call oppdf(n,Tih,ldh_T) ! Tih is uppertri s.t. Tih*Tih' = T^{-1}
    TiF = F
    call dtrmm ('l','u','t','n',n,p,1d0,Tih,n,TiF,n) ! TiF = T^{-1/2}*F
    FTF = betQ0
    call dsyrk ('u','t',p,n,1d0,TiF,n,1d0,FTF,p) ! FTF = F'*T^{-1}*F + Q0
    call dtrmm ('l','u','n','n',n,p,1d0,Tih,n,TiF,n) ! TiF = T^{-1}*F
    call oppdf(p,FTF,ldh_FTF) ! FTF is uptri s.t. FTF*FTF' = (F'*T^{-1}*F + Q0)^{-1}
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
    ! T^{-1}*F*(F'*T^{-1}*F + Q0)^{-1}*F'*T^{-1}
    ldh_Ups = -ldh_T - ldh_FTF
  end subroutine calc_cov

  subroutine calc_cov_pred (z0_ups, TC, FCTF, &
     phi, omg, dmdm0, F, F0, kappa,  T, n, n0, p)
    implicit none
    integer, intent(in) :: n, n0, p
    double precision, intent(in) :: phi, omg, dmdm0(n,n0), F(n,p), &
       F0(n0,p), kappa, T(n,n)
    double precision, intent(out) :: z0_ups(n0), TC(n,n0), FCTF(n0,p)
    double precision C(n,n0), omgp1

    omgp1 = omg + 1d0
    C = dmdm0
    call covmat_a (C,phi,kappa,n,n0)
    call dsymm ('l','u',n,n0,1d0,T,n,C,n,0d0,TC,n)
    z0_ups = sqrt(omgp1 - sum(TC*C,1))
    FCTF = F0
    call dgemm ('t','n',n0,p,n,-1d0,TC,n,F,n,1d0,FCTF,n0)
  end subroutine calc_cov_pred


  subroutine covlist (kg,philist,omglist,n,p,betQ0,F,dm,kappa, &
     Ulist,ldh_Ulist)
!!! Compute list of relevant covariance matrices

    implicit none

    integer, intent(in) :: kg, n, p
    double precision, intent(in) :: philist(kg), betQ0(p,p), omglist(kg), &
       F(n,p), dm(n,n), kappa
    double precision, intent(out) :: ldh_Ulist(kg), Ulist(n,n,kg)
    integer i
    double precision FTF(p,p), T(n,n), TiF(n,p)

    do i = 1, kg
      call calc_cov(philist(i),omglist(i),dm,F,betQ0,&
         kappa,n,p,T,TiF,FTF,Ulist(:,:,i),ldh_Ulist(i))
    end do
  end subroutine covlist

  subroutine covmat_v (dm,phi,kappa,n)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: phi, kappa
    double precision, intent(inout) :: dm(n)
    call covmat_a (dm,phi,kappa,n,1)
  end subroutine covmat_v

  subroutine covmat_a (dm,phi,kappa,n1,n2)
    implicit none
    integer, intent(in) :: n1, n2
    double precision, intent(in) :: phi, kappa
    double precision, intent(inout) :: dm(n1,n2)
    if (phi .eq. 0d0) then
      where (dm .eq. 0d0)
        dm = 1d0
      elsewhere
        dm = 0d0
      end where
    else if (phi .gt. 0d0) then
      dm = dm/phi
      dm = spcor(dm,kappa)
    else
      call rexit ('covmat - Negative phi')
    end if
  end subroutine covmat_a


  subroutine covmat_l (dm,phi,kappa,n1,n2,ldm)
    implicit none
    integer, intent(in) :: n1, n2
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
      where (ldm) dm = dm/phi
      where (ldm) dm = spcor(dm,kappa)
    else
      call rexit ('covmat - Negative phi')
    end if
  end subroutine covmat_l


  subroutine oppdf (n,A,ldh)
    ! INPUT:
    ! A - A positive definite matrix
    ! OUTPUT:
    ! A - An upper triangular matrix U s.t. U*U' = A^{-1}
    ! ldh - .5* log det (A) where A is the original (input) matrix
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: A(n,n)
    double precision, intent(out) :: ldh
    integer i
    ! First call cholesky decomposition
    call dpotrf ('u',n,A,n,i)
    if (i .ne. 0) then
      call rexit ("oppdf - Matrix not positive definite")
    end if
    ldh = 0d0
    do i = 1, n
      ldh = ldh + log(A(i,i))
    end do
    ! Next invert
    call dtrtri ('u','n',n,A,n,i)
    if (i .ne. 0) then
      call rexit ("oppdf - Matrix not invertible")
    end if
  end subroutine oppdf

  subroutine fill_symmetric_matrix (A,n)
    ! If A is symmetric, but only the upper triangle is specified, this
    ! routine returns the complete matrix A.
    integer, intent(in) :: n
    double precision, intent(inout) :: A(n,n)
    integer i
    do i = 1, n-1
      A(i+1:,i) = A(i,i+1:)
    end do
  end subroutine fill_symmetric_matrix
end module covfun

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
  private :: spcor_, spcor_dh_, spcor_dk_, spcor_hh_, spcor_hk_, &
     spcor_dhdk_
  private :: spcor__, spcor_dh__, spcor_dk__, spcor_hh__, spcor_hk__, &
     spcor_dhdk__
  integer, parameter, private :: CORCODES(3) = (/1,2,3/)
  logical, private :: CORDEF = .false.
  integer, private :: CORIS = 0, DIMSP = 0
  logical, private, pointer :: lup(:,:) => null()


  interface covmat
    module procedure covmat_a, covmat_l, covmat_v
  end interface covmat

  abstract interface
    pure double precision function spcor_ (h,k)
      double precision, intent(in) :: h,k
    end function spcor_
  end interface

  abstract interface
    pure double precision function spcor_dh_ (h,k)
      double precision, intent(in) :: h,k
    end function spcor_dh_
  end interface

  abstract interface
    pure double precision function spcor_dk_ (h,k)
      double precision, intent(in) :: h,k
    end function spcor_dk_
  end interface

  abstract interface
    pure double precision function spcor_hh_ (h,k)
      double precision, intent(in) :: h,k
    end function spcor_hh_
  end interface

  abstract interface
    pure double precision function spcor_hk_ (h,k)
      double precision, intent(in) :: h,k
    end function spcor_hk_
  end interface

  abstract interface
    pure double precision function spcor_dhdk_ (h,k)
      double precision, intent(in) :: h,k
    end function spcor_dhdk_
  end interface

  procedure (spcor_     ), pointer :: spcor__      => null()
  procedure (spcor_dh_  ), pointer :: spcor_dh__   => null()
  procedure (spcor_dk_  ), pointer :: spcor_dk__   => null()
  procedure (spcor_hh_  ), pointer :: spcor_hh__   => null()
  procedure (spcor_hk_  ), pointer :: spcor_hk__   => null()
  procedure (spcor_dhdk_), pointer :: spcor_dhdk__ => null()

contains

  function upper_tri () result (lup_)
    logical, allocatable :: lup_(:,:)
    allocate(lup_(DIMSP,DIMSP))
    lup_ = lup
  end function upper_tri

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

      select case (icf)
      case (1) ! Matern
        spcor__      => cor_matern
        spcor_dh__   => cor_dh_matern
        spcor_dk__   => cor_dk_matern
        spcor_hh__   => cor_hh_matern
        spcor_hk__   => cor_hk_matern
        spcor_dhdk__ => cor_dhdk_matern
      case (2) ! Spherical
        spcor__      => cor_spher
        spcor_dh__   => cor_dh_spher
        spcor_dk__   => cor_dk_spher
        spcor_hh__   => cor_hh_spher
        spcor_hk__   => cor_hk_spher
        spcor_dhdk__ => cor_dhdk_spher
      case (3) ! Power exponential
        spcor__      => cor_powexp
        spcor_dh__   => cor_dh_powexp
        spcor_dk__   => cor_dk_powexp
        spcor_hh__   => cor_hh_powexp
        spcor_hk__   => cor_hk_powexp
        spcor_dhdk__ => cor_dhdk_powexp
      case default
        call rexit("Correlation code not used.")
      end select

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


  elemental double precision function spcor (h,k)
    double precision, intent(in) :: h,k
    spcor = spcor__(h,k)
  end function spcor

  elemental double precision function spcor_dh (h,k)
    double precision, intent(in) :: h,k
    spcor_dh = spcor_dh__(h,k)
  end function spcor_dh

  elemental double precision function spcor_dk (h,k)
    double precision, intent(in) :: h,k
    spcor_dk = spcor_dk__(h,k)
  end function spcor_dk

  elemental double precision function spcor_hh (h,k)
    double precision, intent(in) :: h,k
    spcor_hh = spcor_hh__(h,k)
  end function spcor_hh

  elemental double precision function spcor_hk (h,k)
    double precision, intent(in) :: h,k
    spcor_hk = spcor_hk__(h,k)
  end function spcor_hk

  elemental double precision function spcor_dhdk (h,k)
    double precision, intent(in) :: h,k
    spcor_dhdk = spcor_dhdk__(h,k)
  end function spcor_dhdk



  subroutine calc_cov (phi,nsq,dm,F,betQ0, &
     kappa,n,p,T,TiF,FTF,Ups,ldh_Ups)
    implicit none
    integer, intent(in) :: n, p
    double precision, intent(in) :: phi, nsq, dm(n,n), &
       F(n,p), betQ0(p,p), kappa
    double precision, intent(out) :: T(n,n), FTF(p,p), Ups(n,n), &
       ldh_Ups, TiF(n,p)
    integer i, j
    double precision nsqp1, Tih(n,n), TFFT(n,p), ldh_T, ldh_FTF

    nsqp1 = nsq + 1d0
    T = dm
    call covmat_l(T,phi,kappa,n,n,lup)
    do i = 1, n
      T(i,i) = nsqp1
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
     phi, nsq, dmdm0, F, F0, kappa,  T, n, n0, p)
    implicit none
    integer, intent(in) :: n, n0, p
    double precision, intent(in) :: phi, nsq, dmdm0(n,n0), F(n,p), &
       F0(n0,p), kappa, T(n,n)
    double precision, intent(out) :: z0_ups(n0), TC(n,n0), FCTF(n0,p)
    double precision C(n,n0), nsqp1

    nsqp1 = nsq + 1d0
    C = dmdm0
    call covmat_a (C,phi,kappa,n,n0)
    call dsymm ('l','u',n,n0,1d0,T,n,C,n,0d0,TC,n)
    z0_ups = sqrt(nsqp1 - sum(TC*C,1))
    FCTF = F0
    call dgemm ('t','n',n0,p,n,-1d0,TC,n,F,n,1d0,FCTF,n0)
  end subroutine calc_cov_pred


  subroutine covlist (kg,philist,nsqlist,n,p,betQ0,F,dm,kappa, &
     Ulist,ldh_Ulist)
!!! Compute list of relevant covariance matrices

    implicit none

    integer, intent(in) :: kg, n, p
    double precision, intent(in) :: philist(kg), betQ0(p,p), nsqlist(kg), &
       F(n,p), dm(n,n), kappa
    double precision, intent(out) :: ldh_Ulist(kg), Ulist(n,n,kg)
    integer i
    double precision FTF(p,p), T(n,n), TiF(n,p)

    do i = 1, kg
      call calc_cov(philist(i),nsqlist(i),dm,F,betQ0,&
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

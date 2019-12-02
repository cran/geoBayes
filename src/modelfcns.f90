!! This module contains functions specific to each model. The functions
!! are pointers to model-specific functions.

module modelfcns

  integer, private, parameter :: &                 ! Models
     MODELCODES(15) = (/1,2,-2,3,4,5,6,7,-7,8,9,10,11,12,-12/)
  logical, private :: MODELDEF = .false.
  integer, private :: MODELIS = 0

contains

  subroutine create_model (ifam)
    integer, intent(in) :: ifam
    if (MODELDEF .and. (MODELIS .eq. ifam)) return ! Model is already defined
    if (.not. ((ifam .eq. 0) .or. any(ifam .eq. MODELCODES))) &
       call rexit ("Unrecognised family.")
    MODELIS = ifam
    MODELDEF = .true.
  end subroutine create_model

  elemental double precision function logitrwhz (z,d)
    double precision, intent(in) :: z, d
    double precision, parameter :: bigneg = -huge(1d0)
    logitrwhz = invtrwhz(z, d)
    if (logitrwhz .gt. 0d0) then
      logitrwhz = log(logitrwhz)
    else
      logitrwhz = bigneg
    end if
  end function logitrwhz

  function logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
!! log-pdf of z after integrating out beta and ssq.
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, modeldfh, xi(n)
    double precision logpdfz
    double precision Upsz(n), zUz, zmxi(n)
    if (lmxi) then
      zmxi = z - xi
      call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
      zUz = dot_product(zmxi,Upsz) + ssqdfsc
    else
      call dsymv ('u',n,1d0,Ups,n,z,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
      zUz = dot_product(z,Upsz) + ssqdfsc
    end if
    logpdfz = ldh_Ups - modeldfh*log(zUz)
  end function logpdfz

  function logpdfz_dz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, modeldfh, xi(n)
    double precision logpdfz_dz(n)
    double precision Upsz(n), zUz, zmxi(n)
    if (lmxi) then
      zmxi = z - xi
      call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
      zUz = dot_product(zmxi,Upsz) + ssqdfsc
    else
      call dsymv ('u',n,1d0,Ups,n,z,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
      zUz = dot_product(z,Upsz) + ssqdfsc
    end if
    logpdfz_dz = -(2*modeldfh/zUz)*Upsz
  end function logpdfz_dz

  function logpdfz_hz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), Ups(n,n), &
       ldh_Ups, ssqdfsc, modeldfh, xi(n)
    double precision logpdfz_hz(n,n)
    double precision Upsz(n), zUz, zmxi(n), UzzU(n,n), alpha
    if (lmxi) then
      zmxi = z - xi
    else
      zmxi = z
    end if
    call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,Upsz,1) ! Upsz = Ups*(z-xi)
    zUz = dot_product(zmxi,Upsz) + ssqdfsc
    alpha = -2d0/zUz
    UzzU = Ups
    call dsyr('u',n,alpha,Upsz,1,UzzU,n)
    logpdfz_hz = alpha*modeldfh*UzzU
  end function logpdfz_hz

  function logpdfmu (n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, modeldfh
    double precision logpdfmu
    integer i
    double precision z(n)
    ! Linear predictor
    do i = 1, n
      z(i) = flink(mu(i), nu)
    end do
    ! log-likelihood for z
    logpdfmu = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    ! Jacobian
    do i = 1, n
      logpdfmu = logpdfmu - loginvlinkdz(z(i),nu)
    end do
  end function logpdfmu


  function condyz (n, y1, y2, z, nu, tsq)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision condyz
    integer i
    double precision w, lfy
    lfy = 0d0
    do i = 1, n
      w = invlink(z(i),nu)
      lfy = lfy + logpdfy(y1(i), y2(i), w)
    end do
    condyz = lfy/tsq
  end function condyz

  function jointyz (n, z, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: z(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, ssqdfsc, tsq, nu, modeldfh, xi(n)
    double precision jointyz
    double precision lfz, lfy
    lfz = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condyz(n, y, l, z, nu, tsq)
    jointyz = lfz + lfy
  end function jointyz

  pure function condymu (n, y1, y2, mu, tsq)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsq
    double precision condymu
    integer i
    double precision lfy
    lfy = 0d0
    do i = 1, n
      lfy = lfy + logpdfy(y1(i), y2(i), mu(i))
    end do
    condymu = lfy/tsq
  end function condymu

  function jointymu (n, mu, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, tsq, modeldfh
    double precision jointymu
    double precision lfmu, lfy
    lfmu = logpdfmu(n, mu, Ups, ldh_Ups, nu, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condymu(n, y, l, mu, tsq)
    jointymu = lfmu + lfy
  end function jointymu


  ! Distribution functions

  elemental double precision function logpdfy (y1, y2, par)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: y1, y2, par
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      logpdfy = logpdfy_gt(y1,y2,par)
    case (1) ! Gaussian response
      logpdfy = logpdfy_ga(y1,y2,par)
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      logpdfy = logpdfy_bi(y1,y2,par)
    case (6,7,-7) ! Poisson response
      logpdfy = logpdfy_po(y1,y2,par)
    case (8,9) ! Gamma repsponse
      logpdfy = logpdfy_gm(y1,y2,par)
    end select
  end function logpdfy

  elemental double precision function logdffy (y1, y2, p1, p2)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: y1, y2, p1, p2
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      logdffy = logdffy_gt(y1,y2,p1,p2)
    case (1) ! Gaussian response
      logdffy = logdffy_ga(y1,y2,p1,p2)
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      logdffy = logdffy_bi(y1,y2,p1,p2)
    case (6,7,-7) ! Poisson response
      logdffy = logdffy_po(y1,y2,p1,p2)
    case (8,9) ! Gamma repsponse
      logdffy = logdffy_gm(y1,y2,p1,p2)
    end select
  end function logdffy

  elemental double precision function logpdfydlnk (y1, y2, par)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: y1, y2, par
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      logpdfydlnk = logpdfydlnk_gt(y1,y2,par)
    case (1) ! Gaussian response
      logpdfydlnk = logpdfydlnk_ga(y1,y2,par)
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      logpdfydlnk = logpdfydlnk_bi(y1,y2,par)
    case (6,7,-7) ! Poisson response
      logpdfydlnk = logpdfydlnk_po(y1,y2,par)
    case (8,9) ! Gamma repsponse
      logpdfydlnk = logpdfydlnk_gm(y1,y2,par)
    end select
  end function logpdfydlnk

  elemental double precision function logpdfyhlnk (y1, y2, par)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: y1, y2, par
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      logpdfyhlnk = logpdfyhlnk_gt(y1,y2,par)
    case (1) ! Gaussian response
      logpdfyhlnk = logpdfyhlnk_ga(y1,y2,par)
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      logpdfyhlnk = logpdfyhlnk_bi(y1,y2,par)
    case (6,7,-7) ! Poisson response
      logpdfyhlnk = logpdfyhlnk_po(y1,y2,par)
    case (8,9) ! Gamma repsponse
      logpdfyhlnk = logpdfyhlnk_gm(y1,y2,par)
    end select
  end function logpdfyhlnk

  elemental double precision function logpdfy3lnk (y1, y2, par)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: y1, y2, par
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      logpdfy3lnk = logpdfy3lnk_gt(y1,y2,par)
    case (1) ! Gaussian response
      logpdfy3lnk = logpdfy3lnk_ga(y1,y2,par)
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      logpdfy3lnk = logpdfy3lnk_bi(y1,y2,par)
    case (6,7,-7) ! Poisson response
      logpdfy3lnk = logpdfy3lnk_po(y1,y2,par)
    case (8,9) ! Gamma repsponse
      logpdfy3lnk = logpdfy3lnk_gm(y1,y2,par)
    end select
  end function logpdfy3lnk

  elemental double precision function mustart (y1, y2)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: y1, y2
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      mustart = mustart_gt(y1,y2)
    case (1) ! Gaussian response
      mustart = mustart_ga(y1,y2)
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      mustart = mustart_bi(y1,y2)
    case (6,7,-7) ! Poisson response
      mustart = mustart_po(y1,y2)
    case (8,9) ! Gamma repsponse
      mustart = mustart_gm(y1,y2)
    end select
  end function mustart

  elemental double precision function fcntruemu (x)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: x
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      fcntruemu = x
    case (1) ! Gaussian response
      fcntruemu = x
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      fcntruemu = exp(x)
    case (6,7,-7) ! Poisson response
      fcntruemu = exp(x)
    case (8,9) ! Gamma repsponse
      fcntruemu = exp(x)
    end select
  end function fcntruemu

  elemental double precision function invtruemu (x)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: x
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invtruemu = x
    case (1) ! Gaussian response
      invtruemu = x
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      invtruemu = log(x)
    case (6,7,-7) ! Poisson response
      invtruemu = log(x)
    case (8,9) ! Gamma repsponse
      invtruemu = log(x)
    end select
  end function invtruemu

  elemental double precision function fcncum (x)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: x
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      fcncum = fcncum_gt(x)
    case (1) ! Gaussian response
      fcncum = fcncum_ga(x)
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      fcncum = fcncum_bi(x)
    case (6,7,-7) ! Poisson response
      fcncum = fcncum_po(x)
    case (8,9) ! Gamma repsponse
      fcncum = fcncum_gm(x)
    end select
  end function fcncum

  elemental double precision function fcncumd2 (x)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: x
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      fcncumd2 = fcncumd2_gt(x)
    case (1) ! Gaussian response
      fcncumd2 = fcncumd2_ga(x)
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      fcncumd2 = fcncumd2_bi(x)
    case (6,7,-7) ! Poisson response
      fcncumd2 = fcncumd2_po(x)
    case (8,9) ! Gamma repsponse
      fcncumd2 = fcncumd2_gm(x)
    end select
  end function fcncumd2

  elemental double precision function fcncumd3 (x)
    use modelfcns_pdfy
    implicit none
    double precision, intent(in) :: x
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      fcncumd3 = fcncumd3_gt(x)
    case (1) ! Gaussian response
      fcncumd3 = fcncumd3_ga(x)
    case (2,-2,3,4,5,10,11,-12,12) ! Binomial response
      fcncumd3 = fcncumd3_bi(x)
    case (6,7,-7) ! Poisson response
      fcncumd3 = fcncumd3_po(x)
    case (8,9) ! Gamma repsponse
      fcncumd3 = fcncumd3_gm(x)
    end select
  end function fcncumd3


  ! Link functions

  elemental double precision function flink (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      flink = flink_ga(w,d)
    case (1) ! Gaussian boxcox
      flink = flink_ga(w,d)
    case (2,-2) ! Robit
      flink = flink_robit(w,d)
    case (3) ! Logit
      flink = flink_logit(w,d)
    case (4) ! Probit
      flink = flink_probit(w,d)
    case (5) ! Wallace
      flink = flink_wallace(w,d)
    case (6,8) ! Modified boxcox
      flink = flink_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      flink = flink_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      flink = flink_modgev(w,d)
    case (12,-12) ! Binomial GEV
      flink = flink_gev(w,d)
    case (11) ! Binomial modifiedGEV
      flink = flink_modgevns(w,d)
    end select
  end function flink

  elemental double precision function invlink (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invlink = invlink_ga(w,d)
    case (1) ! Gaussian boxcox
      invlink = invlink_ga(w,d)
    case (2,-2) ! Robit
      invlink = invlink_robit(w,d)
    case (3) ! Logit
      invlink = invlink_logit(w,d)
    case (4) ! Probit
      invlink = invlink_probit(w,d)
    case (5) ! Wallace
      invlink = invlink_wallace(w,d)
    case (6,8) ! Modified boxcox
      invlink = invlink_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      invlink = invlink_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      invlink = invlink_modgev(w,d)
    case (12,-12) ! Binomial GEV
      invlink = invlink_gev(w,d)
    case (11) ! Binomial modifiedGEV
      invlink = invlink_modgevns(w,d)
    end select
  end function invlink

  elemental double precision function invlinkdz (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invlinkdz = invlinkdz_ga(w,d)
    case (1) ! Gaussian boxcox
      invlinkdz = invlinkdz_ga(w,d)
    case (2,-2) ! Robit
      invlinkdz = invlinkdz_robit(w,d)
    case (3) ! Logit
      invlinkdz = invlinkdz_logit(w,d)
    case (4) ! Probit
      invlinkdz = invlinkdz_probit(w,d)
    case (5) ! Wallace
      invlinkdz = invlinkdz_wallace(w,d)
    case (6,8) ! Modified boxcox
      invlinkdz = invlinkdz_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      invlinkdz = invlinkdz_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      invlinkdz = invlinkdz_modgev(w,d)
    case (12,-12) ! Binomial GEV
      invlinkdz = invlinkdz_gev(w,d)
    case (11) ! Binomial modifiedGEV
      invlinkdz = invlinkdz_modgevns(w,d)
    end select
  end function invlinkdz

  elemental double precision function loginvlinkdz (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      loginvlinkdz = loginvlinkdz_ga(w,d)
    case (1) ! Gaussian boxcox
      loginvlinkdz = loginvlinkdz_ga(w,d)
    case (2,-2) ! Robit
      loginvlinkdz = loginvlinkdz_robit(w,d)
    case (3) ! Logit
      loginvlinkdz = loginvlinkdz_logit(w,d)
    case (4) ! Probit
      loginvlinkdz = loginvlinkdz_probit(w,d)
    case (5) ! Wallace
      loginvlinkdz = loginvlinkdz_wallace(w,d)
    case (6,8) ! Modified boxcox
      loginvlinkdz = loginvlinkdz_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      loginvlinkdz = loginvlinkdz_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      loginvlinkdz = loginvlinkdz_modgev(w,d)
    case (12,-12) ! Binomial GEV
      loginvlinkdz = loginvlinkdz_gev(w,d)
    case (11) ! Binomial modifiedGEV
      loginvlinkdz = loginvlinkdz_modgevns(w,d)
    end select
  end function loginvlinkdz

  elemental double precision function invlinkhz (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invlinkhz = invlinkhz_ga(w,d)
    case (1) ! Gaussian boxcox
      invlinkhz = invlinkhz_ga(w,d)
    case (2,-2) ! Robit
      invlinkhz = invlinkhz_robit(w,d)
    case (3) ! Logit
      invlinkhz = invlinkhz_logit(w,d)
    case (4) ! Probit
      invlinkhz = invlinkhz_probit(w,d)
    case (5) ! Wallace
      invlinkhz = invlinkhz_wallace(w,d)
    case (6,8) ! Modified boxcox
      invlinkhz = invlinkhz_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      invlinkhz = invlinkhz_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      invlinkhz = invlinkhz_modgev(w,d)
    case (12,-12) ! Binomial GEV
      invlinkhz = invlinkhz_gev(w,d)
    case (11) ! Binomial modifiedGEV
      invlinkhz = invlinkhz_modgevns(w,d)
    end select
  end function invlinkhz

  elemental double precision function invlink3z (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invlink3z = invlink3z_ga(w,d)
    case (1) ! Gaussian boxcox
      invlink3z = invlink3z_ga(w,d)
    case (2,-2) ! Robit
      invlink3z = invlink3z_robit(w,d)
    case (3) ! Logit
      invlink3z = invlink3z_logit(w,d)
    case (4) ! Probit
      invlink3z = invlink3z_probit(w,d)
    case (5) ! Wallace
      invlink3z = invlink3z_wallace(w,d)
    case (6,8) ! Modified boxcox
      invlink3z = invlink3z_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      invlink3z = invlink3z_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      invlink3z = invlink3z_modgev(w,d)
    case (12,-12) ! Binomial GEV
      invlink3z = invlink3z_gev(w,d)
    case (11) ! Binomial modifiedGEV
      invlink3z = invlink3z_modgevns(w,d)
    end select
  end function invlink3z

  elemental double precision function invlinkdn (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invlinkdn = invlinkdn_ga(w,d)
    case (1) ! Gaussian boxcox
      invlinkdn = invlinkdn_ga(w,d)
    case (2,-2) ! Robit
      invlinkdn = invlinkdn_robit(w,d)
    case (3) ! Logit
      invlinkdn = invlinkdn_logit(w,d)
    case (4) ! Probit
      invlinkdn = invlinkdn_probit(w,d)
    case (5) ! Wallace
      invlinkdn = invlinkdn_wallace(w,d)
    case (6,8) ! Modified boxcox
      invlinkdn = invlinkdn_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      invlinkdn = invlinkdn_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      invlinkdn = invlinkdn_modgev(w,d)
    case (12,-12) ! Binomial GEV
      invlinkdn = invlinkdn_gev(w,d)
    case (11) ! Binomial modifiedGEV
      invlinkdn = invlinkdn_modgevns(w,d)
    end select
  end function invlinkdn

  elemental double precision function invlinkhn (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invlinkhn = invlinkhn_ga(w,d)
    case (1) ! Gaussian boxcox
      invlinkhn = invlinkhn_ga(w,d)
    case (2,-2) ! Robit
      invlinkhn = invlinkhn_robit(w,d)
    case (3) ! Logit
      invlinkhn = invlinkhn_logit(w,d)
    case (4) ! Probit
      invlinkhn = invlinkhn_probit(w,d)
    case (5) ! Wallace
      invlinkhn = invlinkhn_wallace(w,d)
    case (6,8) ! Modified boxcox
      invlinkhn = invlinkhn_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      invlinkhn = invlinkhn_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      invlinkhn = invlinkhn_modgev(w,d)
    case (12,-12) ! Binomial GEV
      invlinkhn = invlinkhn_gev(w,d)
    case (11) ! Binomial modifiedGEV
      invlinkhn = invlinkhn_modgevns(w,d)
    end select
  end function invlinkhn

  elemental double precision function invlinkdzdn (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invlinkdzdn = invlinkdzdn_ga(w,d)
    case (1) ! Gaussian boxcox
      invlinkdzdn = invlinkdzdn_ga(w,d)
    case (2,-2) ! Robit
      invlinkdzdn = invlinkdzdn_robit(w,d)
    case (3) ! Logit
      invlinkdzdn = invlinkdzdn_logit(w,d)
    case (4) ! Probit
      invlinkdzdn = invlinkdzdn_probit(w,d)
    case (5) ! Wallace
      invlinkdzdn = invlinkdzdn_wallace(w,d)
    case (6,8) ! Modified boxcox
      invlinkdzdn = invlinkdzdn_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      invlinkdzdn = invlinkdzdn_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      invlinkdzdn = invlinkdzdn_modgev(w,d)
    case (12,-12) ! Binomial GEV
      invlinkdzdn = invlinkdzdn_gev(w,d)
    case (11) ! Binomial modifiedGEV
      invlinkdzdn = invlinkdzdn_modgevns(w,d)
    end select
  end function invlinkdzdn

  elemental double precision function invlinkhzdn (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invlinkhzdn = invlinkhzdn_ga(w,d)
    case (1) ! Gaussian boxcox
      invlinkhzdn = invlinkhzdn_ga(w,d)
    case (2,-2) ! Robit
      invlinkhzdn = invlinkhzdn_robit(w,d)
    case (3) ! Logit
      invlinkhzdn = invlinkhzdn_logit(w,d)
    case (4) ! Probit
      invlinkhzdn = invlinkhzdn_probit(w,d)
    case (5) ! Wallace
      invlinkhzdn = invlinkhzdn_wallace(w,d)
    case (6,8) ! Modified boxcox
      invlinkhzdn = invlinkhzdn_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      invlinkhzdn = invlinkhzdn_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      invlinkhzdn = invlinkhzdn_modgev(w,d)
    case (12,-12) ! Binomial GEV
      invlinkhzdn = invlinkhzdn_gev(w,d)
    case (11) ! Binomial modifiedGEV
      invlinkhzdn = invlinkhzdn_modgevns(w,d)
    end select
  end function invlinkhzdn

  elemental double precision function invlinkdzhn (w, d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invlinkdzhn = invlinkdzhn_ga(w,d)
    case (1) ! Gaussian boxcox
      invlinkdzhn = invlinkdzhn_ga(w,d)
    case (2,-2) ! Robit
      invlinkdzhn = invlinkdzhn_robit(w,d)
    case (3) ! Logit
      invlinkdzhn = invlinkdzhn_logit(w,d)
    case (4) ! Probit
      invlinkdzhn = invlinkdzhn_probit(w,d)
    case (5) ! Wallace
      invlinkdzhn = invlinkdzhn_wallace(w,d)
    case (6,8) ! Modified boxcox
      invlinkdzhn = invlinkdzhn_modbc(w,d)
    case (7,-7,9) ! Standard boxcox
      invlinkdzhn = invlinkdzhn_boxcox(w,d)
    case (10) ! Binomial modifiedGEV
      invlinkdzhn = invlinkdzhn_modgev(w,d)
    case (12,-12) ! Binomial GEV
      invlinkdzhn = invlinkdzhn_gev(w,d)
    case (11) ! Binomial modifiedGEV
      invlinkdzhn = invlinkdzhn_modgevns(w,d)
    end select
  end function invlinkdzhn


  ! Transformation

  elemental double precision function transfw (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      transfw = w
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      transfw = w
    case (-2) ! Binomial robit workaround
      transfw = flink_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      transfw = flink_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      transfw = flink_modgev(w,d)
    end select
  end function transfw

  elemental double precision function invtrw (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invtrw = w
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      invtrw = w
    case (-2) ! Binomial robit workaround
      invtrw = invlink_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      invtrw = invlink_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      invtrw = invlink_modgev(w,d)
    end select
  end function invtrw

  elemental double precision function invtrwdz (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invtrwdz = 1d0
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      invtrwdz = 1d0
    case (-2) ! Binomial robit workaround
      invtrwdz = invlinkdz_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      invtrwdz = invlinkdz_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      invtrwdz = invlinkdz_modgev(w,d)
    end select
  end function invtrwdz

  elemental double precision function loginvtrwdz (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      loginvtrwdz = 0d0
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      loginvtrwdz = 0d0
    case (-2) ! Binomial robit workaround
      loginvtrwdz = loginvlinkdz_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      loginvtrwdz = loginvlinkdz_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      loginvtrwdz = loginvlinkdz_modgev(w,d)
    end select
  end function loginvtrwdz

  elemental double precision function invtrwhz (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invtrwhz = 0d0
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      invtrwhz = 0d0
    case (-2) ! Binomial robit workaround
      invtrwhz = invlinkhz_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      invtrwhz = invlinkhz_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      invtrwhz = invlinkhz_modgev(w,d)
    end select
  end function invtrwhz

  elemental double precision function invtrw3z (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invtrw3z = 0d0
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      invtrw3z = 0d0
    case (-2) ! Binomial robit workaround
      invtrw3z = invlink3z_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      invtrw3z = invlink3z_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      invtrw3z = invlink3z_modgev(w,d)
    end select
  end function invtrw3z

  elemental double precision function invtrwdn (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invtrwdn = 0d0
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      invtrwdn = 0d0
    case (-2) ! Binomial robit workaround
      invtrwdn = invlinkdn_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      invtrwdn = invlinkdn_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      invtrwdn = invlinkdn_modgev(w,d)
    end select
  end function invtrwdn

  elemental double precision function invtrwhn (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invtrwhn = 0d0
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      invtrwhn = 0d0
    case (-2) ! Binomial robit workaround
      invtrwhn = invlinkhn_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      invtrwhn = invlinkhn_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      invtrwhn = invlinkhn_modgev(w,d)
    end select
  end function invtrwhn

  elemental double precision function invtrwdzdn (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invtrwdzdn = 0d0
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      invtrwdzdn = 0d0
    case (-2) ! Binomial robit workaround
      invtrwdzdn = invlinkdzdn_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      invtrwdzdn = invlinkdzdn_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      invtrwdzdn = invlinkdzdn_modgev(w,d)
    end select
  end function invtrwdzdn

  elemental double precision function invtrwhzdn (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invtrwhzdn = 0d0
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      invtrwhzdn = 0d0
    case (-2) ! Binomial robit workaround
      invtrwhzdn = invlinkhzdn_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      invtrwhzdn = invlinkhzdn_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      invtrwhzdn = invlinkhzdn_modgev(w,d)
    end select
  end function invtrwhzdn

  elemental double precision function invtrwdzhn (w,d)
    use modelfcns_link
    implicit none
    double precision, intent(in) :: w, d
    select case (MODELIS)
    case (0) ! Transformed Gaussian
      invtrwdzhn = 0d0
    case (1,2,3,4,5,6,7,8,9,10,11,12) ! No transformation
      invtrwdzhn = 0d0
    case (-2) ! Binomial robit workaround
      invtrwdzhn = invlinkdzhn_wallace(w,d)
    case (-7) ! Poisson boxcox workaround
      invtrwdzhn = invlinkdzhn_modbc(w,d)
    case (-12) ! Binomial GEV workaround
      invtrwdzhn = invlinkdzhn_modgev(w,d)
    end select
  end function invtrwdzhn
end module modelfcns

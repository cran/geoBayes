!! This module contains functions specific to each model. The functions
!! are pointers to model-specific functions.

module modelfcns

  integer, parameter :: &                             ! Models
     MODELCODES(14) = (/1,2,-2,3,4,5,6,7,-7,8,9,10,11,-11/)

  private :: logpdfy_ , logdffy_ , logpdfydlnk_ , logpdfyhlnk_ ,&
     mustart_ , fcntruemu_ , invtruemu_ , flink_ , invlink_ , &
     fcncum_, fcncumd2_, fcncumd3_, &
     invlinkdz_ , invlinkhz_ , invlink3z_ , &
     invlinkdn_ , invlinkhn_ , &
     invlinkdzdn_ , invlinkhzdn_ , invlinkdzhn_ , &
     transfw_ , invtrw_ , invtrwdz_, invtrwhz_, invtrw3z_ , invtrwdn_, &
     invtrwhn_, invtrwdzdn_, invtrwhzdn_ , invtrwdzhn_
  private :: logpdfy__ , logdffy__ , logpdfydlnk__ , logpdfyhlnk__ ,&
     mustart__ , fcntruemu__ , invtruemu__ , flink__ , invlink__ , &
     fcncum__, fcncumd2__, fcncumd3__, &
     invlinkdz__, invlinkhz__, invlink3z__, &
     invlinkdn__, invlinkhn__, &
     invlinkdzdn__, invlinkhzdn__, invlinkdzhn__, &
     transfw__ , invtrw__ , invtrwdz__, invtrwhz__, invtrw3z__, invtrwdn__, &
     invtrwhn__, invtrwdzdn__, invtrwhzdn__, invtrwdzhn__
  private :: self, slog, sexp
  private :: MODELCODES

  logical, private :: MODELDEF = .false.
  integer, private :: MODELIS = 0

  abstract interface
    pure double precision function logpdfy_ (y1, y2, par)
      double precision, intent(in) :: y1, y2, par
    end function logpdfy_
  end interface

  abstract interface
    pure double precision function logdffy_ (y1, y2, p1, p2)
      double precision, intent(in) :: y1, y2, p1, p2
    end function logdffy_
  end interface

  abstract interface
    pure double precision function logpdfydlnk_ (y1, y2, par)
      double precision, intent(in) :: y1, y2, par
    end function logpdfydlnk_
  end interface

  abstract interface
    pure double precision function logpdfyhlnk_ (y1, y2, par)
      double precision, intent(in) :: y1, y2, par
    end function logpdfyhlnk_
  end interface

  abstract interface
    pure double precision function mustart_ (y1, y2)
      double precision, intent(in) :: y1, y2
    end function mustart_
  end interface

  abstract interface
    pure double precision function fcntruemu_ (w)
      double precision, intent(in) :: w
    end function fcntruemu_
  end interface

  abstract interface
    pure double precision function invtruemu_ (mu)
      double precision, intent(in) :: mu
    end function invtruemu_
  end interface

  abstract interface
    pure double precision function fcncum_ (mu)
      double precision, intent(in) :: mu
    end function fcncum_
  end interface

  abstract interface
    pure double precision function fcncumd2_ (mu)
      double precision, intent(in) :: mu
    end function fcncumd2_
  end interface

  abstract interface
    pure double precision function fcncumd3_ (mu)
      double precision, intent(in) :: mu
    end function fcncumd3_
  end interface

  abstract interface
    pure double precision function flink_ (w, d)
      double precision, intent(in) :: w, d
    end function flink_
  end interface

  abstract interface
    pure double precision function invlink_ (z,d)
      double precision, intent(in) :: z, d
    end function invlink_
  end interface

  abstract interface
    pure double precision function invlinkdz_ (z,d)
      double precision, intent(in) :: z, d
    end function invlinkdz_
  end interface

  abstract interface
    pure double precision function invlinkhz_ (z,d)
      double precision, intent(in) :: z, d
    end function invlinkhz_
  end interface

  abstract interface
    pure double precision function invlink3z_ (z,d)
      double precision, intent(in) :: z, d
    end function invlink3z_
  end interface

  abstract interface
    pure double precision function invlinkdn_ (z,d)
      double precision, intent(in) :: z, d
    end function invlinkdn_
  end interface

  abstract interface
    pure double precision function invlinkhn_ (z,d)
      double precision, intent(in) :: z, d
    end function invlinkhn_
  end interface

  abstract interface
    pure double precision function invlinkdzdn_ (z,d)
      double precision, intent(in) :: z, d
    end function invlinkdzdn_
  end interface

  abstract interface
    pure double precision function invlinkhzdn_ (z,d)
      double precision, intent(in) :: z, d
    end function invlinkhzdn_
  end interface

  abstract interface
    pure double precision function invlinkdzhn_ (z,d)
      double precision, intent(in) :: z, d
    end function invlinkdzhn_
  end interface

  abstract interface
    pure double precision function transfw_ (w,d)
      double precision, intent(in) :: w, d
    end function transfw_
  end interface

  abstract interface
    pure double precision function invtrw_ (z,d)
      double precision, intent(in) :: z, d
    end function invtrw_
  end interface

  abstract interface
    pure double precision function invtrwdz_ (z,d)
      double precision, intent(in) :: z, d
    end function invtrwdz_
  end interface

  abstract interface
    pure double precision function invtrwhz_ (z,d)
      double precision, intent(in) :: z, d
    end function invtrwhz_
  end interface

  abstract interface
    pure double precision function invtrw3z_ (z,d)
      double precision, intent(in) :: z, d
    end function invtrw3z_
  end interface

  abstract interface
    pure double precision function invtrwdn_ (z,d)
      double precision, intent(in) :: z, d
    end function invtrwdn_
  end interface

  abstract interface
    pure double precision function invtrwhn_ (z,d)
      double precision, intent(in) :: z, d
    end function invtrwhn_
  end interface

  abstract interface
    pure double precision function invtrwdzdn_ (z,d)
      double precision, intent(in) :: z, d
    end function invtrwdzdn_
  end interface

  abstract interface
    pure double precision function invtrwhzdn_ (z,d)
      double precision, intent(in) :: z, d
    end function invtrwhzdn_
  end interface

  abstract interface
    pure double precision function invtrwdzhn_ (z,d)
      double precision, intent(in) :: z, d
    end function invtrwdzhn_
  end interface

  procedure (logpdfy_    ), pointer :: logpdfy__     => null()
  procedure (logdffy_    ), pointer :: logdffy__     => null()
  procedure (logpdfydlnk_), pointer :: logpdfydlnk__ => null()
  procedure (logpdfyhlnk_), pointer :: logpdfyhlnk__ => null()
  procedure (mustart_    ), pointer :: mustart__     => null()
  procedure (fcntruemu_  ), pointer :: fcntruemu__   => null()
  procedure (invtruemu_  ), pointer :: invtruemu__   => null()
  procedure (fcncum_     ), pointer :: fcncum__      => null()
  procedure (fcncumd2_   ), pointer :: fcncumd2__    => null()
  procedure (fcncumd3_   ), pointer :: fcncumd3__    => null()
  procedure (flink_      ), pointer :: flink__       => null()
  procedure (invlink_    ), pointer :: invlink__     => null()
  procedure (invlinkdz_  ), pointer :: invlinkdz__   => null()
  procedure (invlinkhz_  ), pointer :: invlinkhz__   => null()
  procedure (invlink3z_  ), pointer :: invlink3z__   => null()
  procedure (invlinkdn_  ), pointer :: invlinkdn__   => null()
  procedure (invlinkhn_  ), pointer :: invlinkhn__   => null()
  procedure (invlinkdzdn_), pointer :: invlinkdzdn__ => null()
  procedure (invlinkhzdn_), pointer :: invlinkhzdn__ => null()
  procedure (invlinkdzhn_), pointer :: invlinkdzhn__ => null()
  procedure (transfw_    ), pointer :: transfw__     => null()
  procedure (invtrw_     ), pointer :: invtrw__      => null()
  procedure (invtrwdz_   ), pointer :: invtrwdz__    => null()
  procedure (invtrwhz_   ), pointer :: invtrwhz__    => null()
  procedure (invtrw3z_   ), pointer :: invtrw3z__    => null()
  procedure (invtrwdn_   ), pointer :: invtrwdn__    => null()
  procedure (invtrwhn_   ), pointer :: invtrwhn__    => null()
  procedure (invtrwdzdn_ ), pointer :: invtrwdzdn__  => null()
  procedure (invtrwhzdn_ ), pointer :: invtrwhzdn__  => null()
  procedure (invtrwdzhn_ ), pointer :: invtrwdzhn__  => null()

contains

  subroutine create_model (ifam)
    ! Assign pointer functions to the family and link.
    ! Sets the following pointers
    ! logpdfy__
    ! logdffy__
    ! logpdfydlnk__
    ! logpdfyhlnk__
    ! mustart__
    ! fcntruemu__
    ! invtruemu__
    ! fcncum__
    ! fcncumd2__
    ! fcncumd3__
    ! flink__
    ! invlink__
    ! invlinkdz__
    ! invlinkhz__
    ! invlink3z__
    ! invlinkdn__
    ! invlinkhn__
    ! invlinkdzdn__
    ! invlinkhzdn__
    ! invlinkdzhn__
    ! transfw__
    ! invtrw__
    ! invtrwdz__
    ! invtrwhz__
    ! invtrwdn__
    ! invtrwhn__
    ! invtrwdzdn__
    ! invtrwdzhn__

    use modelfcns_pdfy
    use modelfcns_link
    integer, intent(in) :: ifam

    if (MODELDEF .and. (MODELIS .eq. ifam)) return ! Model is already defined

    if (.not. ((ifam .eq. 0) .or. any(ifam .eq. MODELCODES))) then
      call rexit ("Unrecognised family.")
    end if

    select case (ifam) ! Choice of distribution functions
    case (0) ! Transformed Gaussian, not really used except mustart
      logpdfy__     => logpdfy_ga
      logdffy__     => logdffy_ga
      logpdfydlnk__ => logpdfydlnk_ga
      logpdfyhlnk__ => logpdfyhlnk_ga
      mustart__     => mustart_ga
      fcntruemu__   => self
      invtruemu__   => self
      fcncum__      => fcncum_ga
      fcncumd2__    => fcncumd2_ga
      fcncumd3__    => fcncumd3_ga
    case (1) ! Gaussian response
      logpdfy__     => logpdfy_ga
      logdffy__     => logdffy_ga
      logpdfydlnk__ => logpdfydlnk_ga
      logpdfyhlnk__ => logpdfyhlnk_ga
      mustart__     => mustart_ga
      fcntruemu__   => self
      invtruemu__   => self
      fcncum__      => fcncum_ga
      fcncumd2__    => fcncumd2_ga
      fcncumd3__    => fcncumd3_ga
    case (2,-2,3,4,5,10,11,-11) ! Binomial response
      logpdfy__     => logpdfy_bi
      logdffy__     => logdffy_bi
      logpdfydlnk__ => logpdfydlnk_bi
      logpdfyhlnk__ => logpdfyhlnk_bi
      mustart__     => mustart_bi
      fcntruemu__   => sexp
      invtruemu__   => slog
      fcncum__      => fcncum_bi
      fcncumd2__    => fcncumd2_bi
      fcncumd3__    => fcncumd3_bi
    case (6,7,-7) ! Poisson response
      logpdfy__     => logpdfy_po
      logdffy__     => logdffy_po
      logpdfydlnk__ => logpdfydlnk_po
      logpdfyhlnk__ => logpdfyhlnk_po
      mustart__     => mustart_po
      fcntruemu__   => sexp
      invtruemu__   => slog
      fcncum__      => fcncum_po
      fcncumd2__    => fcncumd2_po
      fcncumd3__    => fcncumd3_po
    case (8,9) ! Gamma repsponse
      logpdfy__     => logpdfy_gm
      logdffy__     => logdffy_gm
      logpdfydlnk__ => logpdfydlnk_gm
      logpdfyhlnk__ => logpdfyhlnk_gm
      fcntruemu__   => sexp
      invtruemu__   => slog
      fcncum__      => fcncum_gm
      fcncumd2__    => fcncumd2_gm
      fcncumd3__    => fcncumd3_gm
    case default
      call rexit("Model code not used.")
    end select

    select case (ifam) ! Link functions
    case (0) ! Transformed Gaussian
      flink__       => flink_ga
      invlink__     => invlink_ga
      invlinkdz__   => invlinkdz_ga
      invlinkhz__   => invlinkhz_ga
      invlink3z__   => invlink3z_ga
      invlinkdn__   => invlinkdn_ga
      invlinkhn__   => invlinkhn_ga
      invlinkdzdn__ => invlinkdzdn_ga
      invlinkhzdn__ => invlinkhzdn_ga
      invlinkdzhn__ => invlinkdzhn_ga
    case (1) ! Gaussian boxcox
      flink__       => flink_ga
      invlink__     => invlink_ga
      invlinkdz__   => invlinkdz_ga
      invlinkhz__   => invlinkhz_ga
      invlink3z__   => invlink3z_ga
      invlinkdn__   => invlinkdn_ga
      invlinkhn__   => invlinkhn_ga
      invlinkdzdn__ => invlinkdzdn_ga
      invlinkhzdn__ => invlinkhzdn_ga
      invlinkdzhn__ => invlinkdzhn_ga
    case (2,-2) ! Robit
      flink__       => flink_robit
      invlink__     => invlink_robit
      invlinkdz__   => invlinkdz_robit
      invlinkhz__   => invlinkhz_robit
      invlink3z__   => invlink3z_robit
      invlinkdn__   => invlinkdn_robit
      invlinkhn__   => invlinkhn_robit
      invlinkdzdn__ => invlinkdzdn_robit
      invlinkhzdn__ => invlinkhzdn_robit
      invlinkdzhn__ => invlinkdzhn_robit
    case (3) ! Logit
      flink__       => flink_logit
      invlink__     => invlink_logit
      invlinkdz__   => invlinkdz_logit
      invlinkhz__   => invlinkhz_logit
      invlink3z__   => invlink3z_logit
      invlinkdn__   => invlinkdn_logit
      invlinkhn__   => invlinkhn_logit
      invlinkdzdn__ => invlinkdzdn_logit
      invlinkhzdn__ => invlinkhzdn_logit
      invlinkdzhn__ => invlinkdzhn_logit
    case (4) ! Probit
      flink__       => flink_probit
      invlink__     => invlink_probit
      invlinkdz__   => invlinkdz_probit
      invlinkhz__   => invlinkhz_probit
      invlink3z__   => invlink3z_probit
      invlinkdn__   => invlinkdn_probit
      invlinkhn__   => invlinkhn_probit
      invlinkdzdn__ => invlinkdzdn_probit
      invlinkhzdn__ => invlinkhzdn_probit
      invlinkdzhn__ => invlinkdzhn_probit
    case (5) ! Wallace
      flink__       => flink_wallace
      invlink__     => invlink_wallace
      invlinkdz__   => invlinkdz_wallace
      invlinkhz__   => invlinkhz_wallace
      invlink3z__   => invlink3z_wallace
      invlinkdn__   => invlinkdn_wallace
      invlinkhn__   => invlinkhn_wallace
      invlinkdzdn__ => invlinkdzdn_wallace
      invlinkhzdn__ => invlinkhzdn_wallace
      invlinkdzhn__ => invlinkdzhn_wallace
    case (6,8) ! Modified boxcox
      flink__       => flink_modbc
      invlink__     => invlink_modbc
      invlinkdz__   => invlinkdz_modbc
      invlinkhz__   => invlinkhz_modbc
      invlink3z__   => invlink3z_modbc
      invlinkdn__   => invlinkdn_modbc
      invlinkhn__   => invlinkhn_modbc
      invlinkdzdn__ => invlinkdzdn_modbc
      invlinkhzdn__ => invlinkhzdn_modbc
      invlinkdzhn__ => invlinkdzhn_modbc
    case (7,-7,9) ! Standard boxcox
      flink__       => flink_boxcox
      invlink__     => invlink_boxcox
      invlinkdz__   => invlinkdz_boxcox
      invlinkhz__   => invlinkhz_boxcox
      invlink3z__   => invlink3z_boxcox
      invlinkdn__   => invlinkdn_boxcox
      invlinkhn__   => invlinkhn_boxcox
      invlinkdzdn__ => invlinkdzdn_boxcox
      invlinkhzdn__ => invlinkhzdn_boxcox
      invlinkdzhn__ => invlinkdzhn_boxcox
    case (10) ! Binomial modifiedGEV
      flink__       => flink_modgev
      invlink__     => invlink_modgev
      invlinkdz__   => invlinkdz_modgev
      invlinkhz__   => invlinkhz_modgev
      invlink3z__   => invlink3z_modgev
      invlinkdn__   => invlinkdn_modgev
      invlinkhn__   => invlinkhn_modgev
      invlinkdzdn__ => invlinkdzdn_modgev
      invlinkhzdn__ => invlinkhzdn_modgev
      invlinkdzhn__ => invlinkdzhn_modgev
    case (11,-11) ! Binomial GEV
      flink__       => flink_gev
      invlink__     => invlink_gev
      invlinkdz__   => invlinkdz_gev
      invlinkhz__   => invlinkhz_gev
      invlink3z__   => invlink3z_gev
      invlinkdn__   => invlinkdn_gev
      invlinkhn__   => invlinkhn_gev
      invlinkdzdn__ => invlinkdzdn_gev
      invlinkhzdn__ => invlinkhzdn_gev
      invlinkdzhn__ => invlinkdzhn_gev
    case default
      call rexit("Model code not used.")
    end select

    select case (ifam) ! Transformation workaround
    case (0) ! Transformed Gaussian
      transfw__     => identity
      invtrw__      => identity
      invtrwdz__    => constant
      invtrwhz__    => zero
      invtrw3z__    => zero
      invtrwdn__    => zero
      invtrwhn__    => zero
      invtrwdzdn__  => zero
      invtrwhzdn__  => zero
      invtrwdzhn__  => zero
    case (1,2,3,4,5,6,7,8,9,10,11) ! No transformation
      transfw__     => identity
      invtrw__      => identity
      invtrwdz__    => constant
      invtrwhz__    => zero
      invtrw3z__    => zero
      invtrwdn__    => zero
      invtrwhn__    => zero
      invtrwdzdn__  => zero
      invtrwhzdn__  => zero
      invtrwdzhn__  => zero
    case (-2) ! Binomial robit workaround
      transfw__     => flink_wallace
      invtrw__      => invlink_wallace
      invtrwdz__    => invlinkdz_wallace
      invtrwhz__    => invlinkhz_wallace
      invtrw3z__    => invlink3z_wallace
      invtrwdn__    => invlinkdn_wallace
      invtrwhn__    => invlinkhn_wallace
      invtrwdzdn__  => invlinkdzdn_wallace
      invtrwhzdn__  => invlinkhzdn_wallace
      invtrwdzhn__  => invlinkdzhn_wallace
    case (-7) ! Poisson boxcox workaround
      transfw__     => flink_modbc
      invtrw__      => invlink_modbc
      invtrwdz__    => invlinkdz_modbc
      invtrwhz__    => invlinkhz_modbc
      invtrw3z__    => invlink3z_modbc
      invtrwdn__    => invlinkdn_modbc
      invtrwhn__    => invlinkhn_modbc
      invtrwdzdn__  => invlinkdzdn_modbc
      invtrwhzdn__  => invlinkhzdn_modbc
      invtrwdzhn__  => invlinkdzhn_modbc
    case (-11) ! Binomial GEV workaround
      transfw__     => flink_modgev
      invtrw__      => invlink_modgev
      invtrwdz__    => invlinkdz_modgev
      invtrwhz__    => invlinkhz_modgev
      invtrw3z__    => invlink3z_modgev
      invtrwdn__    => invlinkdn_modgev
      invtrwhn__    => invlinkhn_modgev
      invtrwdzdn__  => invlinkdzdn_modgev
      invtrwhzdn__  => invlinkhzdn_modgev
      invtrwdzhn__  => invlinkdzhn_modgev
    case default
      call rexit("Model code not used.")
    end select

    MODELIS = ifam
    MODELDEF = .true.
  end subroutine create_model

  elemental double precision function logpdfy (y1, y2, par)
    double precision, intent(in) :: y1, y2, par
    !if (associated(logpdfy__))
    logpdfy = logpdfy__(y1, y2, par)
  end function logpdfy

  elemental double precision function logdffy (y1, y2, p1, p2)
    double precision, intent(in) :: y1, y2, p1, p2
    !if (associated(logdffy__))
    logdffy = logdffy__(y1, y2, p1, p2)
  end function logdffy

  elemental double precision function logpdfydlnk (y1, y2, par)
    double precision, intent(in) :: y1, y2, par
    !if (associated(logpdfydlnk__))
    logpdfydlnk = logpdfydlnk__(y1, y2, par)
  end function logpdfydlnk

  elemental double precision function logpdfyhlnk (y1, y2, par)
    double precision, intent(in) :: y1, y2, par
    !if (associated(logpdfyhlnk__))
    logpdfyhlnk = logpdfyhlnk__(y1, y2, par)
  end function logpdfyhlnk

  elemental double precision function mustart (y1, y2)
    double precision, intent(in) :: y1, y2
    !if (associated(mustart__))
    mustart = mustart__(y1, y2)
  end function mustart

  elemental double precision function fcntruemu (w)
    double precision, intent(in) :: w
    !if (associated(fcntruemu__))
    fcntruemu = fcntruemu__(w)
  end function fcntruemu

  elemental double precision function invtruemu (mu)
    double precision, intent(in) :: mu
    !if (associated(invtruemu__))
    invtruemu = invtruemu__(mu)
  end function invtruemu

  elemental double precision function fcncum (mu)
    double precision, intent(in) :: mu
    !if (associated(fcncum__))
    fcncum = fcncum__(mu)
  end function fcncum

  elemental double precision function fcncumd2 (mu)
    double precision, intent(in) :: mu
    !if (associated(fcncumd2__))
    fcncumd2 = fcncumd2__(mu)
  end function fcncumd2

  elemental double precision function fcncumd3 (mu)
    double precision, intent(in) :: mu
    !if (associated(fcncumd3__))
    fcncumd3 = fcncumd3__(mu)
  end function fcncumd3

  elemental double precision function flink (w, d)
    double precision, intent(in) :: w, d
    !if (associated(flink__))
    flink = flink__(w, d)
  end function flink

  elemental double precision function invlink (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invlink__))
    invlink = invlink__(z, d)
  end function invlink

  elemental double precision function invlinkdz (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invlinkdz__))
    invlinkdz = invlinkdz__(z, d)
  end function invlinkdz

  elemental double precision function invlinkhz (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invlinkhz__))
    invlinkhz = invlinkhz__(z, d)
  end function invlinkhz

  elemental double precision function invlink3z (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invlink3z__))
    invlink3z = invlink3z__(z, d)
  end function invlink3z

  elemental double precision function invlinkdn (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invlinkdn__))
    invlinkdn = invlinkdn__(z, d)
  end function invlinkdn

  elemental double precision function invlinkhn (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invlinkhn__))
    invlinkhn = invlinkhn__(z, d)
  end function invlinkhn

  elemental double precision function invlinkdzdn (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invlinkdzdn__))
    invlinkdzdn = invlinkdzdn__(z, d)
  end function invlinkdzdn

  elemental double precision function invlinkhzdn (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invlinkhzdn__))
    invlinkhzdn = invlinkhzdn__(z, d)
  end function invlinkhzdn

  elemental double precision function invlinkdzhn (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invlinkdzhn__))
    invlinkdzhn = invlinkdzhn__(z, d)
  end function invlinkdzhn

  elemental double precision function transfw (w,d)
    double precision, intent(in) :: w, d
    !if (associated(transfw__))
    transfw = transfw__(w, d)
  end function transfw

  elemental double precision function invtrw (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invtrw__))
    invtrw = invtrw__(z, d)
  end function invtrw

  elemental double precision function invtrwdz (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invtrwdz__))
    invtrwdz = invtrwdz__(z, d)
  end function invtrwdz

  elemental double precision function invtrwhz (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invtrwhz__))
    invtrwhz = invtrwhz__(z, d)
  end function invtrwhz

  elemental double precision function invtrw3z (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invtrw3z__))
    invtrw3z = invtrw3z__(z, d)
  end function invtrw3z

  elemental double precision function invtrwdn (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invtrwdn__))
    invtrwdn = invtrwdn__(z, d)
  end function invtrwdn

  elemental double precision function invtrwhn (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invtrwhn__))
    invtrwhn = invtrwhn__(z, d)
  end function invtrwhn

  elemental double precision function invtrwdzdn (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invtrwdzdn__))
    invtrwdzdn = invtrwdzdn__(z, d)
  end function invtrwdzdn

  elemental double precision function invtrwhzdn (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invtrwhzdn__))
    invtrwhzdn = invtrwhzdn__(z, d)
  end function invtrwhzdn

  elemental double precision function invtrwdzhn (z,d)
    double precision, intent(in) :: z, d
    !if (associated(invtrwdzhn__))
    invtrwdzhn = invtrwdzhn__(z, d)
  end function invtrwdzhn

!   TODO
!   elemental double precision function logilinkdz (z,d)
!     double precision, intent(in) :: z, d
!     if (associated(logilinkdz__)) logilinkdz = logilinkdz__(z, d)
!   end function logilinkdz

  elemental double precision function logilinkdz (z,d)
    double precision, intent(in) :: z, d
    double precision, parameter :: bigneg = -huge(1d0)
    logilinkdz = invlinkdz(z, d)
    if (logilinkdz .gt. 0d0) then
      logilinkdz = log(logilinkdz)
    else
      logilinkdz = bigneg
    end if
  end function logilinkdz

  elemental double precision function logitrwdz (z,d)
    double precision, intent(in) :: z, d
    double precision, parameter :: bigneg = -huge(1d0)
    logitrwdz = invtrwdz(z, d)
    if (logitrwdz .gt. 0d0) then
      logitrwdz = log(logitrwdz)
    else
      logitrwdz = bigneg
    end if
  end function logitrwdz

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

  subroutine logpdfzdz (fc, gr, z, Ups, ldh_Ups, xi, lmxi, ssq, n)
!! log-pdf of z and its derivative after integrating out beta. The 2*pi
!! constant is removed.
    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: lmxi
    double precision, intent(in) :: z(n), Ups(n,n), ldh_Ups, xi(n), ssq
    double precision, intent(out) :: fc, gr(n)
    double precision zmxi(n)
    if (lmxi) then
      zmxi = z - xi
    else
      zmxi = z
    end if
    call dsymv ('u',n,1d0,Ups,n,zmxi,1,0d0,gr,1) ! gr = Ups*(z-xi)
    gr = -gr/ssq ! gr = -Ups*(z-x)/ssq
    fc = -.5d0*n*log(ssq) + ldh_Ups + .5d0*dot_product(zmxi,gr)
  end subroutine logpdfzdz

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
      logpdfmu = logpdfmu - logilinkdz(z(i),nu)
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

  function condymu (n, y1, y2, mu, tsq)
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

  subroutine logcondyzdz (fc, gr, nu, y1, y2, z, n, tsq)
!! log-pdf of y|z and its derivative w.r.t. z.
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz
    fc = 0d0
    do i = 1, n
      par = invlink(z(i),nu)
      pardz = invlinkdz(z(i),nu)
      fc = fc + logpdfy(y1(i),y2(i),par)
      gr(i) = logpdfydlnk(y1(i),y2(i),par) * pardz
    end do
    fc = fc/tsq
    gr = gr/tsq
  end subroutine logcondyzdz

  subroutine logcondyzhs (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision par, pardz
    do i = 1, n
      par = invlink(z(i),nu)
      pardz = invlinkdz(z(i),nu)
      hs(i) = logpdfyhlnk(y1(i),y2(i),par) * pardz*pardz
    end do
    hs = hs/tsq
  end subroutine logcondyzhs

  pure double precision function self (x)
    double precision, intent(in) :: x
    self = x
  end function self

  pure double precision function slog (x)
    double precision, intent(in) :: x
    slog = log(x)
  end function slog

  pure double precision function sexp (x)
    double precision, intent(in) :: x
    sexp = exp(x)
  end function sexp

end module modelfcns

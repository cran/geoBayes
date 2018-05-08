!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Derivatives of pdfs w.r.t. z.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module pdfdz
  
contains

  subroutine logcondyzdz_gt (fc, gr, nu, y1, y2, z, n, tsq)
!! TODO respdf is not needed here. Can be computed outside the fcn. XXX
!! log-pdf of y|z and its derivative w.r.t. z.
!! Gaussian transformed
!! tsq is tsqdf*tsqsc on input
!! respdf is n + tsqdf on input
    use linkfcns, only: invlink_ga
    use linkdz, only: invlinkdz_ga
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz
    fc = tsq
    do i = 1, n
      par = invlink_ga(z(i),nu)
      pardz = invlinkdz_ga(z(i),nu)
      par = y1(i) - par
      gr(i) = y2(i)*par
      fc = fc + gr(i)*par
      gr(i) = gr(i)*pardz
    end do
    gr = gr/fc ! XXX Missing respdf factor. Must be computed outside the fcn
    fc = -.5d0*log(fc) ! XXX Also missing respdf = n + tsqdf
  end subroutine logcondyzdz_gt

  subroutine logcondyzhs_gt (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
!! Gaussian transformed
!! tsq is tsqdf*tsqsc on input
    use linkfcns, only: invlink_ga
    use linkdz, only: invlinkdz_ga
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision par, pardz, fc
    fc = tsq
    do i = 1, n
      par = invlink_ga(z(i),nu)
      pardz = invlinkdz_ga(z(i),nu)
      hs(i) = y2(i)*par*pardz*pardz
      par = y1(i) - par
      fc = fc + y2(i)*par*par
    end do
    hs = -hs/fc ! XXX Missing respdf factor. Must be computed outside the fcn
  end subroutine logcondyzhs_gt

  subroutine logcondyzdz_ga (fc, gr, nu, y1, y2, z, n, tsq)
!! log-pdf of y|z and its derivative w.r.t. z.
!! Gaussian 
    use linkfcns, only: invlink_ga
    use linkdz, only: invlinkdz_ga
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz, tmp
    fc = 0d0
    do i = 1, n
      par = invlink_ga(z(i),nu)
      pardz = invlinkdz_ga(z(i),nu)
      tmp = y2(i)*par
      fc = fc + y1(i)*par - .5d0*tmp*par
      gr(i) = (y1(i) - tmp)*pardz
    end do
    fc = fc/tsq
    gr = gr/tsq
  end subroutine logcondyzdz_ga

  subroutine logcondyzhs_ga (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
!! Gaussian 
    use linkdz, only: invlinkdz_ga
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision pardz
    do i = 1, n
      pardz = invlinkdz_ga(z(i),nu)
      hs(i) = y2(i)*pardz*pardz
    end do
    hs = -hs/tsq
  end subroutine logcondyzhs_ga

  subroutine logcondyzdz_bi (fc, gr, nu, y1, y2, z, n, tsq)
!! log-pdf of y|z and its derivative w.r.t. z.
!! Binomial robit
    use interfaces, only: fexpm1, flog1mexp
    use linkfcns, only: invlink_bi
    use linkdz, only: invlinkdz_bi
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz
    fc = 0d0
    do i = 1, n
      par = invlink_bi(z(i),nu)
      pardz = invlinkdz_bi(z(i),nu)
      fc = fc + y1(i)*par + y2(i)*flog1mexp(par)
      gr(i) = (y1(i) - y2(i)/fexpm1(-par))*pardz
    end do
    fc = fc/tsq
    gr = gr/tsq
  end subroutine logcondyzdz_bi

  subroutine logcondyzhs_bi (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
!! Binomial robit
    use interfaces, only: fexpm1, flog1mexp
    use linkfcns, only: invlink_bi
    use linkdz, only: invlinkdz_bi
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision par, pardz, tmp
    do i = 1, n
      par = invlink_bi(z(i),nu)
      pardz = invlinkdz_bi(z(i),nu)
      tmp = 1d0/fexpm1(-par)
      hs(i) = y2(i)*tmp*(1d0+tmp)*pardz*pardz
    end do
    hs = -hs/tsq
  end subroutine logcondyzhs_bi

  subroutine logcondyzdz_po (fc, gr, nu, y1, y2, z, n, tsq)
!! log-pdf of y|z and its derivative w.r.t. z.
!! Poisson
    use linkfcns, only: invlink_po
    use linkdz, only: invlinkdz_po
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz, tmp
    fc = 0d0
    do i = 1, n
      par = invlink_po(z(i),nu)
      pardz = invlinkdz_po(z(i),nu)
      tmp = y2(i)*exp(par)
      fc = fc + y1(i)*par - tmp
      gr(i) = (y1(i) - tmp)*pardz
    end do
    fc = fc/tsq
    gr = gr/tsq
  end subroutine logcondyzdz_po

  subroutine logcondyzhs_po (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
!! Poisson
    use linkfcns, only: invlink_po
    use linkdz, only: invlinkdz_po
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision par, pardz, tmp
    do i = 1, n
      par = invlink_po(z(i),nu)
      pardz = invlinkdz_po(z(i),nu)
      hs(i) = y2(i)*exp(par)*pardz*pardz
    end do
    hs = -hs/tsq
  end subroutine logcondyzhs_po

  subroutine logcondyzdz_gm (fc, gr, nu, y1, y2, z, n, tsq)
!! log-pdf of y|z and its derivative w.r.t. z.
!! Gamma
    use linkfcns, only: invlink_gm
    use linkdz, only: invlinkdz_gm
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz, tmp
    fc = 0d0
    do i = 1, n
      par = invlink_gm(z(i),nu)
      pardz = invlinkdz_gm(z(i),nu)
      tmp = y1(i)*exp(-par)
      fc = fc - tmp - y2(i)*par
      gr(i) = (tmp - y2(i))*pardz
    end do
    fc = fc/tsq
    gr = gr/tsq
  end subroutine logcondyzdz_gm

  subroutine logcondyzhs_gm (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
!! Gamma
    use linkfcns, only: invlink_gm
    use linkdz, only: invlinkdz_gm
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision par, pardz
    do i = 1, n
      par = invlink_gm(z(i),nu)
      pardz = invlinkdz_gm(z(i),nu)
      hs(i) = y1(i)*exp(-par)*pardz*pardz
    end do
    hs = -hs/tsq
  end subroutine logcondyzhs_gm

  subroutine logcondyzdz_ba (fc, gr, nu, y1, y2, z, n, tsq)
!! log-pdf of y|z and its derivative w.r.t. z.
!! Binomial asymmetric
    use interfaces, only: fexpm1, flog1mexp
    use linkfcns, only: invlink_ba
    use linkdz, only: invlinkdz_ba
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz
    fc = 0d0
    do i = 1, n
      par = invlink_ba(z(i),nu)
      pardz = invlinkdz_ba(z(i),nu)
      fc = fc + y1(i)*par + y2(i)*flog1mexp(par)
      gr(i) = (y1(i) - y2(i)*fexpm1(-par))*pardz
    end do
    fc = fc/tsq
    gr = gr/tsq
  end subroutine logcondyzdz_ba

  subroutine logcondyzhs_ba (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
!! Binomial asymmetric
    use interfaces, only: fexpm1, flog1mexp
    use linkfcns, only: invlink_ba
    use linkdz, only: invlinkdz_ba
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision par, pardz, tmp
    do i = 1, n
      par = invlink_ba(z(i),nu)
      pardz = invlinkdz_ba(z(i),nu)
      tmp = 1d0/fexpm1(-par)
      hs(i) = y2(i)*tmp*(1d0+tmp)*pardz*pardz
    end do
    hs = -hs/tsq
  end subroutine logcondyzhs_ba

  subroutine logcondyzdz_bd (fc, gr, nu, y1, y2, z, n, tsq)
!! log-pdf of y|z and its derivative w.r.t. z.
!! Binomial asymmetric decreasing
    use interfaces, only: fexpm1, flog1mexp
    use linkfcns, only: invlink_bd
    use linkdz, only: invlinkdz_bd
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz
    fc = 0d0
    do i = 1, n
      par = invlink_bd(z(i),nu)
      pardz = invlinkdz_bd(z(i),nu)
      fc = fc + y2(i)*par + y1(i)*flog1mexp(par)
      gr(i) = (y2(i) - y1(i)*fexpm1(-par))*pardz
    end do
    fc = fc/tsq
    gr = gr/tsq
  end subroutine logcondyzdz_bd

  subroutine logcondyzhs_bd (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
!! Binomial asymmetric decreasing
    use interfaces, only: fexpm1, flog1mexp
    use linkfcns, only: invlink_bd
    use linkdz, only: invlinkdz_bd
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision par, pardz, tmp
    fc = 0d0
    do i = 1, n
      par = invlink_bd(z(i),nu)
      pardz = invlinkdz_bd(z(i),nu)
      tmp = 1d0/fexpm1(-par)
      hs(i) = y1(i)*tmp*(1d0+tmp)*pardz*pardz
    end do
    hs = -hs/tsq
  end subroutine logcondyzhs_bd

  subroutine logcondyzdz_bw (fc, gr, nu, y1, y2, z, n, tsq)
!! log-pdf of y|z and its derivative w.r.t. z.
!! Binomial Wallace
    use interfaces, only: fexpm1, flog1mexp
    use linkfcns, only: invlink_bw
    use linkdz, only: invlinkdz_bw
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: fc, gr(n)
    integer i
    double precision par, pardz
    fc = 0d0
    do i = 1, n
      par = invlink_bw(z(i),nu)
      pardz = invlinkdz_bw(z(i),nu)
      fc = fc + y1(i)*par + y2(i)*flog1mexp(par)
      gr(i) = (y1(i) - y2(i)*fexpm1(-par))*pardz
    end do
    fc = fc/tsq
    gr = gr/tsq
  end subroutine logcondyzdz_bw

  subroutine logcondyzhs_bw (hs, nu, y1, y2, z, n, tsq)
!! Component used in the calculation of the Hessian
!! Binomial Wallace
    use interfaces, only: fexpm1, flog1mexp
    use linkfcns, only: invlink_bw
    use linkdz, only: invlinkdz_bw
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), z(n), nu, tsq
    double precision, intent(out) :: hs(n)
    integer i
    double precision par, pardz, tmp
    fc = 0d0
    do i = 1, n
      par = invlink_bw(z(i),nu)
      pardz = invlinkdz_bw(z(i),nu)
      tmp = 1d0/fexpm1(-par)
      hs(i) = y2(i)*tmp*(1d0+tmp)*pardz*pardz
    end do
    hs = -hs/tsq
  end subroutine logcondyzhs_bw
end module pdfdz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!! Commentary: Alternative binomial
!!! Use the transformation
!!! w = sign(z) * (1/c) * sqrt(nu*log1p(z*z/nu))
!!! z = sign(w) * sqrt(nu*(expm1((c*w)*(c*w)/nu)))
!!! where c = (8*nu+3)/(8*nu+1)
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module transfbinomial
  implicit none
contains
  function condyw_bi (n, y, l, w, nu, tsq)
    use condyz, only: condyz_bi
    use interfaces, only: fexpm1, quantnorm
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: w(n), y(n), l(n), nu, tsq
    double precision condyw_bi
    double precision z(n), cnu, ccqqn, cnusqnu, qw
    integer i
    cnu = 1d0 + 2d0/(8d0*nu+1d0)
    cnusqnu = cnu*cnu/nu
    do i = 1, n
      qw = w(i)
      if (qw .eq. 0d0) then
        z(i) = 0d0
      else
        ccqqn = cnusqnu*qw*qw
        z(i) = sign(sqrt(nu*fexpm1(ccqqn)), qw)
      end if
    end do
    condyw_bi = condyz_bi(n, y, l, z, nu, tsq)
  end function condyw_bi

  function logpdfw_bi (n, w, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, modeldfh)
    use pdfz, only: logpdfz
    use interfaces, only: flogexpm1, fexpm1, quantnorm
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: w(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, modeldfh
    double precision logpdfw_bi
    double precision lfw, logjac, z(n), cnu, ccqqn, nulogh, &
       cnusqnu, cnulog, nulgmc, tmp, qw, qwsq
    integer i
    cnu = 1d0 + 2d0/(8d0*nu+1d0)
    cnusqnu = cnu*cnu/nu
    cnulog = log(cnu)
    nulogh = .5d0*log(nu)
    nulgmc = nulogh - cnulog
    logjac = 0d0
    do i = 1, n
      qw = w(i)
      qwsq = qw*qw
      if (qw .eq. 0d0) then
        z(i) = 0d0
        tmp = nulgmc
      else
        ccqqn = cnusqnu*qwsq
        z(i) = sign(sqrt(nu*fexpm1(ccqqn)), qw)
        tmp = ccqqn - .5*flogexpm1(ccqqn)
      end if
      logjac = logjac + tmp
    end do
    logjac = logjac + n*(cnulog - nulgmc)
    lfw = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    logpdfw_bi = lfw + logjac
  end function logpdfw_bi

  function jointyw_bi (n, w, y, l, Ups, ldh_Ups, &
     nu, xi, lmxi, ssqdfsc, tsq, modeldfh)
    use condyz, only: condyz_bi
    use pdfz, only: logpdfz
    use interfaces, only: flogexpm1, fexpm1, quantnorm
    implicit none
    logical, intent(in) :: lmxi
    integer, intent(in) :: n
    double precision, intent(in) :: w(n), y(n), l(n), Ups(n, n), &
       ldh_Ups, nu, xi(n), ssqdfsc, tsq, modeldfh
    double precision jointyw_bi
    double precision lfw, lfy, logjac, z(n), cnu, ccqqn, nulogh, &
       cnusqnu, cnulog, nulgmc, tmp, qw, qwsq
    integer i
    cnu = 1d0 + 2d0/(8d0*nu+1d0)
    cnusqnu = cnu*cnu/nu
    cnulog = log(cnu)
    nulogh = .5d0*log(nu)
    nulgmc = nulogh - cnulog
    logjac = 0d0
    do i = 1, n
      qw = w(i) ! quantnorm(w(i))
      qwsq = qw*qw
      if (qw .eq. 0d0) then
        z(i) = 0d0
        tmp = nulgmc
      else
        ccqqn = cnusqnu*qwsq
        z(i) = sign(sqrt(nu*fexpm1(ccqqn)), qw)
        tmp = ccqqn - .5*flogexpm1(ccqqn) !+ log(abs(qw))
      end if
      logjac = logjac + tmp !+ .5d0*qwsq
    end do
    logjac = logjac + n*(cnulog - nulgmc)
    lfw = logpdfz(n, z, Ups, ldh_Ups, xi, lmxi, ssqdfsc, modeldfh)
    lfy = condyz_bi(n, y, l, z, nu, tsq)
    jointyw_bi = lfw + lfy + logjac
  end function jointyw_bi
end module transfbinomial

module calcbd_fcns

  private :: condymu_sp, weigh_llik_st, weigh_llik_cv, weigh_llik_deriv_st,&
     weigh_llik_deriv_cv

contains

  pure function condymuf (ifam, n, y1, y2, mu, tsqval, respdfh)
    use condymu, only: condymu_gt
    implicit none
    integer, intent(in) :: n, ifam
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

  pure function condymu_sp (n, y1, y2, mu, tsqval, respdfh)
    use modelfcns, only: condymu_mf => condymu
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: y1(n), y2(n), mu(n), tsqval, &
       respdfh
    double precision condymu_sp
    condymu_sp = condymu_mf(n,y1,y2,mu,tsqval)
  end function condymu_sp

  pure double precision function weigh_llik (icv, llik, lw, Qcv, n)
    implicit none
    integer, intent(in) :: n, icv
    double precision, intent(in) :: llik(n), lw(n), Qcv(n)
    select case (icv)
    case (0) ! No control variates
      weigh_llik = weigh_llik_st(llik,lw,Qcv,n)
    case (1) ! With control variates
      weigh_llik = weigh_llik_cv(llik,lw,Qcv,n)
    end select
  end function weigh_llik

  pure double precision function weigh_llik_st (llik, lw, Qcv, n) result (out)
    ! out = log(sum(exp(llik - lw)))
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: llik(n), lw(n), Qcv(n)
    double precision llikw(n)
    llikw = llik - lw
    out = maxval(llikw)
    llikw = exp(llikw - out)
    out = out + log(sum(llikw))
  end function weigh_llik_st

  pure double precision function weigh_llik_cv (llik, lw, Qcv, n) result (out)
    ! out = log(Qcv'*(n*exp(llik - lw)))
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: llik(n), lw(n), Qcv(n)
    double precision llikw(n), ql
    llikw = llik - lw
    out = maxval(llikw)
    llikw = exp(llikw - out)
    ql = dot_product(Qcv,llikw)
    if (ql .gt. 0d0) then
      out = out + log(n*ql)
    else
      out = -huge(1d0)
    end if
  end function weigh_llik_cv

  pure double precision function weigh_llik_deriv (icv, dllik, llik, lw, Qcv, n)
    implicit none
    integer, intent(in) :: n, icv
    double precision, intent(in) :: dllik(n), llik(n), lw(n), Qcv(n)
    select case (icv)
    case (0) ! No control variates
      weigh_llik_deriv = weigh_llik_deriv_st(dllik,llik,lw,Qcv,n)
    case (1) ! With control variates
      weigh_llik_deriv = weigh_llik_deriv_cv(dllik,llik,lw,Qcv,n)
    end select
  end function weigh_llik_deriv

  pure double precision function weigh_llik_deriv_st (dllik, llik, lw, Qcv, n) &
     result (out)
    ! out = log(sum(exp(llik - lw)))
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: dllik(n), llik(n), lw(n), Qcv(n)
    double precision llikw(n)
    llikw = llik - lw
    out = maxval(llikw)
    llikw = exp(llikw - out)
    llikw = llikw/sum(llikw)
    out = dot_product(dllik,llikw)
  end function weigh_llik_deriv_st

  pure double precision function weigh_llik_deriv_cv (dllik, llik, lw, Qcv, n) &
     result (out)
    ! out = log(Qcv'*(n*exp(llik - lw)))
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: dllik(n), llik(n), lw(n), Qcv(n)
    double precision llikw(n)
    llikw = llik - lw
    out = maxval(llikw)
    llikw = exp(llikw - out)
    llikw = Qcv*llikw
    llikw = llikw/sum(llikw)
    out = dot_product(dllik,llikw)
  end function weigh_llik_deriv_cv

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

  double precision function cor_dcov (n, dm, phi, nsq, kappa, id) result (DT)
    ! id = which derivative do you want? phi, nsq, kappa
    ! Return dT where T = nsq*I + R(phi,kappa)
    use covfun, only: spcor_dh, spcor_dk, fill_symmetric_matrix, upper_tri
    implicit none
    integer, intent(in) :: n, id
    double precision, intent(in) :: dm(n,n), phi, nsq, kappa
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
end module calcbd_fcns

module bmargin
contains
  subroutine revlogistic (eta,llik,kg,Ntot,Nout)

    ! use lbfgsqck
    use lbfgsbmod
    use flogsumexp, only: logrsumexp

    implicit none

    integer, intent(in) :: kg, Ntot, Nout(kg)
    double precision, intent(in) :: llik(Ntot,kg)
    double precision, intent(inout) :: eta(kg)
    double precision lgdenom(Ntot), logp(Ntot,kg), lliketa(Ntot,kg)
    integer, parameter :: maxit=500, iprint(2) = (/-1,0/)
    ! double precision, parameter :: epstol = 5d-5
    double precision, parameter :: pgtol = 0d-5, factr = 1d7
    double precision neglogfun, neglogfun_grad(kg)
    integer ia, ie, i
    integer iflag
    logical lidx(Ntot,kg)
    double precision etalo(kg), etaup(kg)
    integer etabd(kg)


    lidx = .false.
    ia = 1
    ie = 0
    do i = 1, kg
      ie = ie + Nout(i)
      lidx(ia:ie,i) = .true. 
      ia = ia + Nout(i)
    end do

    etalo = eta
    etaup = eta
    etabd = 0
    ! etabd(1) = 2
    iflag = 0
    do i = 1, maxit
      lliketa = llik + spread(eta,1,Ntot)
      lgdenom = logrsumexp(lliketa,Ntot,kg)
      logp = spread(lgdenom,2,kg)
      logp = lliketa - logp
      neglogfun = -sum(logp,lidx)
      neglogfun_grad = -dble(Nout) + sum(exp(logp),1)
      !       call lbfgs (kg-1,eta(2:kg),neglogfun,neglogfun_grad(2:kg),&
      !          iprint,epstol,iflag)
      call lbfgsb (kg,eta,etalo,etaup,etabd,neglogfun,neglogfun_grad,&
         iprint(1),factr,pgtol,iflag)
      if (iflag == 0) then
        !       print '(" Convergence in ",I0," iterations")', i
        exit
      else if (iflag < 0) then
        call rwarn ('The reverse logistic agorithm didn''t converge')
        exit
      end if
    end do
    if (iflag > 0) then
      call rwarn ('The reverse logistic algorithm needs more iterations')
    end if
    eta = eta - eta(1) + log(dble(Nout))
  end subroutine revlogistic


  subroutine mengwong (eta,llik,kg,Ntot,Nout)

    use flogsumexp, only: logrsumexp

    implicit none

    integer, intent(in) :: kg, Ntot, Nout(kg)
    double precision, intent(in) :: llik(Ntot,kg)
    double precision, intent(inout) :: eta(kg)
    integer, parameter :: maxit=500
    double precision, parameter :: pgtol = 5d-5
    integer ia, ie, i, j, it, ipiv(kg-1)
    double precision llik_mix_all(Ntot), a(kg), Bmat(2:kg,2:kg), bvec(2:kg), &
       r(kg), likw(Ntot,kg)
    logical convergence
    double precision bnorm, work(4*kg-4), rcond
    integer iwork(kg-1)
    double precision, external :: dlange

    a = dble(Nout)
    r = exp(eta)/a
    a = a/dble(Ntot)

    iterations: do it = 1, maxit
      call rchkusr ! Check if user has requested interrupt
      ia = 1
      ie = 0
      likw = spread(log(a*r),1,Ntot)
      likw = likw + llik
      llik_mix_all = -logrsumexp(likw,Ntot,kg)
      likw = spread(llik_mix_all,2,kg)
      likw = llik + likw - log(dble(Ntot))
      likw = exp(likw)
      ia = 1 + Nout(1)
      ie = Nout(1)
      do i = 2, kg
        ie = ie + Nout(i)
        forall (j=2:kg, j .ne. i)
          Bmat(i,j) = -a(j)*sum(likw(ia:ie,j))
        end forall
        bvec(i) = a(1)*sum(likw(ia:ie,1))
        ia = ia + Nout(i)
      end do
      do i = 2, kg
        ia = 1 
        ie = 0
        Bmat(i,i) = 0d0
        do j = 1, kg
          ie = ie + Nout(j)
          if (j .ne. i) then
            Bmat(i,i) = Bmat(i,i) + a(i)*sum(likw(ia:ie,i))
          end if
          ia = ia + Nout(j)
        end do
      end do
      j = kg - 1
      bnorm = dlange('1',j,j,Bmat,j,work) ! Matrix norm used to check singularity
      call dgesv (j,1,Bmat,j,ipiv,bvec,j,i)
      if (i .ne. 0) then
        call rwarn ('mengwong - Cannot find solution to linear system. &
           &Previous solution restored.')
        exit iterations
      end if
      if (any(bvec .le. 0d0)) then
        call rwarn ('mengwong - Non-positive solution to the linear system. &
           &Previous solution restored.')
        exit iterations
      end if
      call dgecon ('1',j-1,Bmat,j-1,bnorm,rcond,work,iwork,i)
      if (rcond .lt. 1d-17) then
        call rwarn ('mengwong - System is computationally singular. &
           &Previous solution restored.')
        exit iterations
      end if
      convergence = sum(abs(bvec - r(2:kg))) < pgtol
      r(2:kg) = bvec
      if (convergence) exit iterations
    end do iterations
    if (.not. convergence) then
      call rwarn ('mengwong - Algorithm did''t converge.')
    end if
    eta = log(dble(Nout)) + log(r)
  end subroutine mengwong
end module bmargin

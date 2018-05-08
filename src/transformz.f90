subroutine transformz (sample, nu, n, ifam)
  use modelfcns
  integer, intent(in) :: n, ifam
  double precision, intent(in) :: nu
  double precision, intent(inout) :: sample(n)
  integer i
  call create_model (ifam)
  do i = 1, n
    sample(i) = invtrw(sample(i),nu)
  end do
end subroutine transformz

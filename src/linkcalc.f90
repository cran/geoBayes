subroutine flinkfcn (z, n, mu, linkp, ifam)
  use modelfcns, only: create_model, invtruemu, flink
  implicit none
  integer, intent(in) :: n, ifam
  double precision, intent(in) :: mu(n), linkp
  double precision, intent(out) :: z(n)
  call create_model(ifam)
  z = invtruemu(mu)
  z = flink(z,linkp)
end subroutine flinkfcn

subroutine flinkinv (mu, n, z, linkp, ifam)
  use modelfcns, only: create_model, fcntruemu, invlink
  implicit none
  integer, intent(in) :: n, ifam
  double precision, intent(in) :: z(n), linkp
  double precision, intent(out) :: mu(n)
  call create_model(ifam)
  mu = invlink(z,linkp)
  mu = fcntruemu(mu)
end subroutine flinkinv

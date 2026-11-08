module vmath_mod
  ! exp, expm1 and log on arrays: plain loops (compilers with a vector math
  ! library vectorize them, e.g. Intel's SVML); src/vmath_accelerate.f90 on
  ! macOS.
  use kinds_mod, only: dp
  use eos_math, only: expm1
  implicit none
  private
  public :: vexp, vexpm1, vlog
contains
  subroutine vexp(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: y(n)
    integer :: i
    do i = 1, n
      y(i) = exp(x(i))
    end do
  end subroutine vexp
  subroutine vexpm1(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: y(n)
    integer :: i
    do i = 1, n
      y(i) = expm1(x(i))
    end do
  end subroutine vexpm1
  subroutine vlog(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: y(n)
    integer :: i
    do i = 1, n
      y(i) = log(x(i))
    end do
  end subroutine vlog
end module vmath_mod

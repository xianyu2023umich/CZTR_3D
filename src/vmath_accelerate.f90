module vmath_mod
  ! exp, expm1 and log on arrays (macOS: Apple's Accelerate framework, vvexp
  ! etc.; src/vmath_plain.f90 elsewhere). Results agree with the scalar
  ! library functions to the last digit or so, not bit for bit.
  use kinds_mod, only: dp
  implicit none
  private
  public :: vexp, vexpm1, vlog
  interface
    subroutine vvexp(y, x, n) bind(C, name='vvexp')
      import :: dp
      real(dp), intent(out) :: y(*)
      real(dp), intent(in) :: x(*)
      integer, intent(in) :: n
    end subroutine vvexp
    subroutine vvexpm1(y, x, n) bind(C, name='vvexpm1')
      import :: dp
      real(dp), intent(out) :: y(*)
      real(dp), intent(in) :: x(*)
      integer, intent(in) :: n
    end subroutine vvexpm1
    subroutine vvlog(y, x, n) bind(C, name='vvlog')
      import :: dp
      real(dp), intent(out) :: y(*)
      real(dp), intent(in) :: x(*)
      integer, intent(in) :: n
    end subroutine vvlog
  end interface
contains
  subroutine vexp(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: y(n)
    call vvexp(y, x, n)
  end subroutine vexp
  subroutine vexpm1(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: y(n)
    call vvexpm1(y, x, n)
  end subroutine vexpm1
  subroutine vlog(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: y(n)
    call vvlog(y, x, n)
  end subroutine vlog
end module vmath_mod

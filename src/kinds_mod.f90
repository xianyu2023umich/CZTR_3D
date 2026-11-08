module kinds_mod
  ! Single place to control working precision.
  use iso_fortran_env, only: real64
  implicit none
  private
  integer, parameter, public :: dp = real64
end module kinds_mod

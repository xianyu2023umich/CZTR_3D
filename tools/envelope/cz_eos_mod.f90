module cz_eos_mod
  ! Single EOS interface for the envelope code. Every EOS use (surface
  ! point, structure equations, mixing-length inputs, entropy target) goes
  ! through this module, and one run-time switch selects the backend:
  !
  !   EOS_FREEEOS  the FreeEOS table (eos_table_mod + thermo_mod), the
  !                original behavior and the default
  !   EOS_HHE      the eos_hhe library physics core (exact, no table),
  !                gas-only
  !
  ! Programs pick the backend from the command line with cz_eos_from_args:
  !   ./main.x                  FreeEOS (default)
  !   ./main.x --eos=hhe        eos_hhe, composition X=0.75 Z=0.02 (as the table)
  !   ./main.x --eos=hhe --X=0.70 --Z=0.02
  !
  ! Entropy zero points differ between backends, so the deep-CZ target is
  ! backend-specific: FreeEOS keeps its tabulated value logS_deep_target;
  ! eos_hhe evaluates its own entropy at the same physical point
  ! (logQ_target, logT_target in const_mod).
  use kinds_mod, only: dp
  use const_mod, only: logS_deep_target, logQ_target, logT_target
  use eos_table_mod, only: freeeos_table_t => eos_table_t, freeeos_state_t => eos_state_t, &
      eos_load, eos_get_rhoT, logrho_of_QT
  use thermo_mod, only: rho_of_PT
  use eos_types, only: hhe_state_t => eos_state, eos_composition
  use eos_thermo, only: hhe_from_rho_T => eos_get_from_rho_T, hhe_from_p_T => eos_get_from_p_T
  implicit none
  private

  integer, parameter, public :: EOS_FREEEOS = 1, EOS_HHE = 2

  type, public :: cz_eos_t
    integer :: backend = EOS_FREEEOS
    character(len=256) :: table_path = ''
    type(freeeos_table_t) :: table          ! FreeEOS backend only
    type(eos_composition) :: comp           ! eos_hhe backend only
    logical :: coulomb = .true.             ! eos_hhe backend only
  end type cz_eos_t

  ! What the envelope code needs from the EOS, identical for both backends.
  ! logS is log10 of the specific entropy in the backend's own zero point.
  type, public :: cz_eos_state_t
    real(dp) :: rho = 0, T = 0, P = 0, logS = 0
    real(dp) :: chiRho = 0, chiT = 0, Cp = 0, grad_ad = 0
  end type cz_eos_state_t

  public :: cz_eos_from_args, cz_eos_init, cz_eos_rhoT, cz_eos_PT, cz_eos_logS_target
  public :: cz_eos_describe, cz_eos_tag

contains

  ! Parse --eos=freeeos|hhe, --X=..., --Z=... (Y = 1 - X - Z) from the command
  ! line and initialize. default_table is the FreeEOS file to use when the
  ! FreeEOS backend is selected.
  subroutine cz_eos_from_args(eos, default_table)
    type(cz_eos_t), intent(out) :: eos
    character(*), intent(in) :: default_table
    character(len=256) :: arg
    integer :: i, backend, ios
    real(dp) :: X, Z

    backend = EOS_FREEEOS
    X = 0.75_dp
    Z = 0.02_dp
    do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      if (arg(1:6) == '--eos=') then
        select case (trim(arg(7:)))
        case ('freeeos'); backend = EOS_FREEEOS
        case ('hhe');     backend = EOS_HHE
        case default
          print *, 'cz_eos: unknown backend "', trim(arg(7:)), '" (use --eos=freeeos or --eos=hhe)'
          stop 1
        end select
      else if (arg(1:4) == '--X=') then
        read(arg(5:), *, iostat=ios) X
        if (ios /= 0) stop 'cz_eos: cannot read --X='
      else if (arg(1:4) == '--Z=') then
        read(arg(5:), *, iostat=ios) Z
        if (ios /= 0) stop 'cz_eos: cannot read --Z='
      end if
    end do
    call cz_eos_init(eos, backend, default_table, eos_composition(X=X, Y=1.0_dp - X - Z, Z=Z))
  end subroutine cz_eos_from_args

  subroutine cz_eos_init(eos, backend, table_path, comp, coulomb)
    type(cz_eos_t), intent(out) :: eos
    integer, intent(in) :: backend
    character(*), intent(in) :: table_path
    type(eos_composition), intent(in), optional :: comp
    logical, intent(in), optional :: coulomb

    eos%backend = backend
    eos%table_path = table_path
    if (present(comp)) eos%comp = comp
    if (present(coulomb)) eos%coulomb = coulomb
    if (backend == EOS_FREEEOS) call eos_load(eos%table, table_path)
    call cz_eos_describe(eos)
  end subroutine cz_eos_init

  subroutine cz_eos_describe(eos)
    type(cz_eos_t), intent(in) :: eos
    if (eos%backend == EOS_FREEEOS) then
      print '(A,A)', 'EOS backend: FreeEOS table ', trim(eos%table_path)
    else
      print '(A,3(A,F7.4),A,L1)', 'EOS backend: eos_hhe physics core (gas only),', &
          ' X=', eos%comp%X, ' Y=', eos%comp%Y, ' Z=', eos%comp%Z, '  Coulomb=', eos%coulomb
    end if
  end subroutine cz_eos_describe

  ! Short tag for output file names ('' for FreeEOS keeps the original names).
  function cz_eos_tag(eos) result(tag)
    type(cz_eos_t), intent(in) :: eos
    character(len=8) :: tag
    tag = ''
    if (eos%backend == EOS_HHE) tag = '_hhe'
  end function cz_eos_tag

  function cz_eos_rhoT(eos, rho, T) result(st)
    type(cz_eos_t), intent(in) :: eos
    real(dp), intent(in) :: rho, T
    type(cz_eos_state_t) :: st
    type(freeeos_state_t) :: ft
    type(hhe_state_t) :: hs

    if (eos%backend == EOS_FREEEOS) then
      ft = eos_get_rhoT(eos%table, rho, T)
      st = from_freeeos(ft, rho, T)
    else
      call hhe_from_rho_T(rho, T, hs, eos%comp, eos%coulomb)
      st = from_hhe(hs)
    end if
  end function cz_eos_rhoT

  ! State at given pressure and temperature; rho_guess is a warm start.
  function cz_eos_PT(eos, P, T, rho_guess) result(st)
    type(cz_eos_t), intent(in) :: eos
    real(dp), intent(in) :: P, T, rho_guess
    type(cz_eos_state_t) :: st
    type(freeeos_state_t) :: ft
    type(hhe_state_t) :: hs
    real(dp) :: rho
    integer :: ierr

    if (eos%backend == EOS_FREEEOS) then
      rho = rho_of_PT(eos%table, P, T, rho_guess, ft)
      st = from_freeeos(ft, rho, T)
    else
      call hhe_from_p_T(P, T, hs, eos%comp, eos%coulomb, rho_guess=rho_guess, ierr=ierr)
      if (ierr /= 0) print '(A,ES12.4,A,ES12.4)', 'cz_eos: eos_hhe (P,T) inversion did not converge at P=', &
          P, ' T=', T
      st = from_hhe(hs)
    end if
  end function cz_eos_PT

  ! Deep-CZ entropy target in the backend's own entropy convention.
  function cz_eos_logS_target(eos) result(logS)
    type(cz_eos_t), intent(in) :: eos
    real(dp) :: logS
    type(cz_eos_state_t) :: st

    if (eos%backend == EOS_FREEEOS) then
      logS = logS_deep_target
    else
      st = cz_eos_rhoT(eos, 10.0_dp**logrho_of_QT(logQ_target, logT_target), 10.0_dp**logT_target)
      logS = st%logS
    end if
  end function cz_eos_logS_target

  function from_freeeos(ft, rho, T) result(st)
    type(freeeos_state_t), intent(in) :: ft
    real(dp), intent(in) :: rho, T
    type(cz_eos_state_t) :: st
    ! rho is the caller's value, exactly as the original code used it.
    st%rho = rho
    st%T = T
    st%P = 10.0_dp**ft%logPgas
    st%logS = ft%logS
    st%chiRho = ft%chiRho
    st%chiT = ft%chiT
    st%Cp = ft%Cp
    st%grad_ad = ft%grad_ad
  end function from_freeeos

  function from_hhe(hs) result(st)
    type(hhe_state_t), intent(in) :: hs
    type(cz_eos_state_t) :: st
    st%rho = hs%rho
    st%T = hs%T
    st%P = hs%p
    st%logS = log10(hs%s)
    st%chiRho = hs%chi_rho
    st%chiT = hs%chi_T
    st%Cp = hs%cp
    st%grad_ad = hs%grad_ad
  end function from_hhe

end module cz_eos_mod

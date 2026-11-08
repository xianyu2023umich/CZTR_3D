module eos_table_mod
  ! Reader + bilinear interpolator for the FreeEOS table.
  !
  ! File format (data/FreeEOS_v51_x0.75_z0.02.dat):
  !   line 1: "EOS"
  !   line 2: number of header lines before the column-name row (informational)
  !   line 3: "<nQ> <nT> <ncol>"
  !   line 4: column names
  !   then nQ*nT rows of ncol reals, logT varying fastest, logQ slowest.
  !
  ! IMPORTANT: this table's "Q" is NOT the standard OPAL R = rho/T6^3.
  ! Empirically (checked against the table's own logRho column at three
  ! independent grid points): logRho = logQ + 2*logT - 12, i.e. Q = rho/T6^2.
  ! Do not reuse logQ as if it were the opacity table's logR -- convert
  ! through logRho first (see eos_get_rhoT / opacity_table_mod).
  use kinds_mod, only: dp
  implicit none
  private

  integer, parameter :: NCOL = 23
  ! column indices into eos_table_t%data(:, iT, iQ)
  integer, parameter, public :: COL_logPgas = 3
  integer, parameter, public :: COL_logE    = 4
  integer, parameter, public :: COL_logS    = 5
  integer, parameter, public :: COL_chiRho  = 6
  integer, parameter, public :: COL_chiT    = 7
  integer, parameter, public :: COL_Cp      = 8
  integer, parameter, public :: COL_Cv      = 9
  integer, parameter, public :: COL_mu      = 13
  integer, parameter, public :: COL_log10ne = 14
  integer, parameter, public :: COL_gamma1  = 15
  integer, parameter, public :: COL_gamma3  = 16
  integer, parameter, public :: COL_grad_ad = 17
  integer, parameter, public :: COL_logRho  = 20

  type, public :: eos_table_t
    integer  :: nQ = 0, nT = 0
    real(dp) :: Qmin = 0, dQ = 0, Tmin = 0, dT = 0
    real(dp), allocatable :: dat(:,:,:)   ! (NCOL, nT, nQ)
  end type eos_table_t

  type, public :: eos_state_t
    real(dp) :: logPgas, logE, logS, chiRho, chiT, Cp, Cv
    real(dp) :: mu, log10ne, gamma1, gamma3, grad_ad, logRho
  end type eos_state_t

  public :: eos_load, eos_get, eos_get_rhoT, logrho_of_QT, logQ_of_rhoT

contains

  subroutine eos_load(tab, path)
    type(eos_table_t), intent(out) :: tab
    character(*), intent(in) :: path
    integer :: unit, nQ, nT, ncol_check, ios
    character(len=256) :: line
    real(dp), allocatable :: logQ_axis(:), logT_axis(:)

    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'eos_load: cannot open ', trim(path)
      stop 1
    end if

    read(unit,'(A)') line   ! "EOS"
    read(unit,'(A)') line   ! number of header lines (unused)
    read(unit,*) nQ, nT, ncol_check
    read(unit,'(A)') line   ! column-name row

    if (ncol_check /= NCOL) then
      write(*,*) 'eos_load: unexpected column count', ncol_check, 'expected', NCOL
      stop 1
    end if

    tab%nQ = nQ
    tab%nT = nT
    allocate(tab%dat(NCOL, nT, nQ))
    allocate(logQ_axis(nQ), logT_axis(nT))

    ! bulk read: file order is logT fastest, logQ slowest, matching (NCOL,nT,nQ)
    read(unit,*) tab%dat
    close(unit)

    logQ_axis = tab%dat(1, 1, :)
    logT_axis = tab%dat(2, :, 1)
    tab%Qmin = logQ_axis(1)
    tab%dQ   = logQ_axis(2) - logQ_axis(1)
    tab%Tmin = logT_axis(1)
    tab%dT   = logT_axis(2) - logT_axis(1)

    write(*,'(A,I0,A,I0,A)') 'eos_load: loaded ', nQ, ' x ', nT, ' grid'
    write(*,'(A,F10.4,A,F10.4,A,F10.6)') '  logQ in [', logQ_axis(1), ', ', logQ_axis(nQ), &
         '],  step ', tab%dQ
    write(*,'(A,F10.4,A,F10.4,A,F10.6)') '  logT in [', logT_axis(1), ', ', logT_axis(nT), &
         '],  step ', tab%dT
  end subroutine eos_load

  ! ρ/T6^2 convention used by this specific table (see module header note)
  pure function logrho_of_QT(logQ, logT) result(logRho)
    real(dp), intent(in) :: logQ, logT
    real(dp) :: logRho
    logRho = logQ + 2.0_dp*logT - 12.0_dp
  end function logrho_of_QT

  pure function logQ_of_rhoT(logRho, logT) result(logQ)
    real(dp), intent(in) :: logRho, logT
    real(dp) :: logQ
    logQ = logRho - 2.0_dp*logT + 12.0_dp
  end function logQ_of_rhoT

  subroutine locate(x, xmin, dx, n, i0, frac)
    real(dp), intent(in)  :: x, xmin, dx
    integer,  intent(in)  :: n
    integer,  intent(out) :: i0
    real(dp), intent(out) :: frac
    real(dp) :: pos
    pos = (x - xmin)/dx + 1.0_dp
    i0 = floor(pos)
    if (i0 < 1)   i0 = 1
    if (i0 > n-1) i0 = n-1
    frac = pos - real(i0, dp)
    frac = max(0.0_dp, min(1.0_dp, frac))
  end subroutine locate

  function eos_get(tab, logQ, logT) result(st)
    type(eos_table_t), intent(in) :: tab
    real(dp), intent(in) :: logQ, logT
    type(eos_state_t) :: st
    integer :: iQ, iT
    real(dp) :: fQ, fT, row(NCOL)

    call locate(logQ, tab%Qmin, tab%dQ, tab%nQ, iQ, fQ)
    call locate(logT, tab%Tmin, tab%dT, tab%nT, iT, fT)

    row = (1-fQ)*(1-fT)*tab%dat(:, iT,   iQ)   &
        + (1-fQ)*fT    *tab%dat(:, iT+1, iQ)   &
        + fQ    *(1-fT)*tab%dat(:, iT,   iQ+1) &
        + fQ    *fT    *tab%dat(:, iT+1, iQ+1)

    st%logPgas = row(COL_logPgas)
    st%logE    = row(COL_logE)
    st%logS    = row(COL_logS)
    st%chiRho  = row(COL_chiRho)
    st%chiT    = row(COL_chiT)
    st%Cp      = row(COL_Cp)
    st%Cv      = row(COL_Cv)
    st%mu      = row(COL_mu)
    st%log10ne = row(COL_log10ne)
    st%gamma1  = row(COL_gamma1)
    st%gamma3  = row(COL_gamma3)
    st%grad_ad = row(COL_grad_ad)
    st%logRho  = row(COL_logRho)
  end function eos_get

  function eos_get_rhoT(tab, rho, T) result(st)
    type(eos_table_t), intent(in) :: tab
    real(dp), intent(in) :: rho, T
    type(eos_state_t) :: st
    real(dp) :: logQ, logT
    logT = log10(T)
    logQ = logQ_of_rhoT(log10(rho), logT)
    st = eos_get(tab, logQ, logT)
  end function eos_get_rhoT

end module eos_table_mod

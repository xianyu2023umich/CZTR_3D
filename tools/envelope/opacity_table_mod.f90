module opacity_table_mod
  ! Reader + bilinear interpolator for the merged (OPAL + Ferguson05) Rosseland
  ! mean opacity table, data/opacity_merged_x70_z02.dat.
  !
  ! File format:
  !   line 1: "OPACITY"
  !   line 2: number of header lines (informational)
  !   line 3: "<nT> <nR> 3"
  !   line 4: column names "logT logR logkappa"
  !   then nT*nR rows of (logT, logR, logkappa), logR varying fastest.
  !
  ! logR == log10(rho / T6^3), T6 = T/1e6  (standard OPAL convention).
  ! Cells outside either source table's valid range are flagged -99 (see
  ! scripts/merge_opacity.py); opac_get_logkappa reports those via out_of_table.
  use kinds_mod, only: dp
  implicit none
  private

  real(dp), parameter :: MISSING_SENTINEL = -50.0_dp   ! anything <= this is "no data"

  type, public :: opac_table_t
    integer  :: nT = 0, nR = 0
    real(dp) :: Tmin = 0, dT = 0, Rmin = 0, dR = 0
    real(dp), allocatable :: logk(:,:)   ! (nR, nT)
  end type opac_table_t

  public :: opac_load, opac_get_logkappa, opac_get_kappa_rhoT

contains

  subroutine opac_load(tab, path)
    type(opac_table_t), intent(out) :: tab
    character(*), intent(in) :: path
    integer :: unit, nT, nR, ncol_check, ios
    character(len=256) :: line
    real(dp), allocatable :: buf(:,:,:)   ! (3, nR, nT): logT, logR, logkappa
    real(dp), allocatable :: Taxis(:), Raxis(:)

    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'opac_load: cannot open ', trim(path)
      stop 1
    end if

    read(unit,'(A)') line   ! "OPACITY"
    read(unit,'(A)') line   ! header-line count (unused)
    read(unit,*) nT, nR, ncol_check
    read(unit,'(A)') line   ! column-name row

    if (ncol_check /= 3) then
      write(*,*) 'opac_load: unexpected column count', ncol_check
      stop 1
    end if

    tab%nT = nT
    tab%nR = nR
    allocate(buf(3, nR, nT))
    allocate(tab%logk(nR, nT))
    allocate(Taxis(nT), Raxis(nR))

    read(unit,*) buf   ! R fastest, T slowest -- matches (3,nR,nT)
    close(unit)

    tab%logk = buf(3, :, :)
    Taxis = buf(1, 1, :)
    Raxis = buf(2, :, 1)
    tab%Tmin = Taxis(1);  tab%dT = Taxis(2) - Taxis(1)
    tab%Rmin = Raxis(1);  tab%dR = Raxis(2) - Raxis(1)

    write(*,'(A,I0,A,I0,A)') 'opac_load: loaded ', nT, ' x ', nR, ' grid'
    write(*,'(A,F10.4,A,F10.4,A,F10.6)') '  logT in [', Taxis(1), ', ', Taxis(nT), &
         '],  step ', tab%dT
    write(*,'(A,F10.4,A,F10.4,A,F10.6)') '  logR in [', Raxis(1), ', ', Raxis(nR), &
         '],  step ', tab%dR
  end subroutine opac_load

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

  subroutine opac_get_logkappa(tab, logT, logR, logk, out_of_table)
    type(opac_table_t), intent(in)  :: tab
    real(dp), intent(in)  :: logT, logR
    real(dp), intent(out) :: logk
    logical,  intent(out) :: out_of_table
    integer :: iT, iR
    real(dp) :: fT, fR
    real(dp) :: c00, c10, c01, c11

    call locate(logT, tab%Tmin, tab%dT, tab%nT, iT, fT)
    call locate(logR, tab%Rmin, tab%dR, tab%nR, iR, fR)

    c00 = tab%logk(iR,   iT)
    c10 = tab%logk(iR+1, iT)
    c01 = tab%logk(iR,   iT+1)
    c11 = tab%logk(iR+1, iT+1)

    out_of_table = (c00 <= MISSING_SENTINEL) .or. (c10 <= MISSING_SENTINEL) .or. &
                   (c01 <= MISSING_SENTINEL) .or. (c11 <= MISSING_SENTINEL)

    logk = (1-fR)*(1-fT)*c00 + fR*(1-fT)*c10 + (1-fR)*fT*c01 + fR*fT*c11
  end subroutine opac_get_logkappa

  ! Convenience wrapper taking physical rho [g/cm^3], T [K]; returns kappa in cm^2/g.
  function opac_get_kappa_rhoT(tab, rho, T, out_of_table) result(kappa)
    type(opac_table_t), intent(in) :: tab
    real(dp), intent(in) :: rho, T
    logical, intent(out) :: out_of_table
    real(dp) :: kappa, logT, logR, logk
    logT = log10(T)
    logR = log10(rho) - 3.0_dp*(logT - 6.0_dp)
    call opac_get_logkappa(tab, logT, logR, logk, out_of_table)
    kappa = 10.0_dp**logk
  end function opac_get_kappa_rhoT

end module opacity_table_mod

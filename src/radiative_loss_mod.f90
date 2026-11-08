module radiative_loss_mod
  ! Chromospheric (Carlsson & Leenaarts 2012) + coronal (CHIANTI, optically
  ! thin) radiative loss, smoothly blended, plus the companion hydrogen
  ! conductivity correction -- all three lookup-table sets taken directly
  ! from the HYDRAD code (rice-solar-physics/HYDRAD), which implements
  ! exactly this recipe (see OpticallyThickIon.cpp, traced by hand and
  ! ported here; not re-derived from the paper's figures).
  !
  ! Combination formula (matches HYDRAD's GetVolumetricLossRate exactly,
  ! itself Carlsson & Leenaarts Eq. 1):
  !   Loss_X = Q_X(T) * Abund_X * IonFrac_X(T) * EscProb_X(depth) * N_H * n_e*rho
  ! where N_H is "H nuclei per gram of stellar material" (a constant fixed by
  ! the abundance mixture), Abund_X*IonFrac_X are pre-multiplied into a single
  ! "combined emissivity" table at load time (same as HYDRAD's ppEmiss), and
  ! depth is log10(column mass) for Ca II/Mg II or log10(4e-14*N_HI_column)
  ! for H I (the depth variable HYDRAD actually feeds its H table -- an
  ! approximate Lyman-alpha-center optical depth built from the neutral-H
  ! column, not plain column mass).
  !
  ! Deliberate departure from HYDRAD: instead of a hard on/off switch at
  ! T=2.4e4 K between this chromospheric treatment and the coronal Lambda(T)
  ! treatment, we blend smoothly in log(T) (see blend_weight) -- the same
  ! lesson the pre-existing n_thick cutoff in this codebase already learned
  ! the hard way (see tr_equilibrium_20260911.f90's header comments).
  use kinds_mod, only: dp
  implicit none
  private

  public :: species_table_t, load_species_table, species_emiss_ionfrac, species_escape_prob
  public :: load_abundance_mixture, N_H_per_gram_global
  public :: coronal_table_t, load_coronal_table, coronal_lambda
  public :: hconduct_table_t, load_hconduct_table, extra_h_kappa0
  public :: blend_weight

  real(dp) :: N_H_per_gram_global = 0.0_dp   ! H nuclei per gram, set once by load_abundance_mixture

  type :: species_table_t
    integer :: n_em = 0, n_esc = 0
    real(dp), allocatable :: logT_em(:), Qcombined_em(:)   ! Q(T)*Abund*IonFrac(T), precombined
    real(dp), allocatable :: x_esc(:), prob_esc(:)
    ! kept separately so callers needing just the ionisation fraction (e.g.
    ! to build the neutral-H column for the H I depth variable) can get it
    integer :: n_bal = 0
    real(dp), allocatable :: logT_bal(:), frac_bal(:)
  end type species_table_t

  type :: coronal_table_t
    integer :: n = 0
    real(dp) :: logTmin = 0, dlogT = 0
    real(dp), allocatable :: loglambda(:)
  end type coronal_table_t

  type :: hconduct_table_t
    integer :: n = 0
    real(dp), allocatable :: logT(:), log10kappa(:)
  end type hconduct_table_t

contains

  !--------------------------------------------------------------
  ! Reads an element/isotope abundance list (HYDRAD .ab format: "Z  log10Ab"
  ! pairs relative to H=12.00, terminated by a "-1" line) plus the matching
  ! atomic-mass file, and computes (a) Abund_Z, the abundance of element Z
  ! relative to H by number, for each of the three species we need, and
  ! (b) N_H_per_gram_global, H nuclei per gram of material for this mixture
  ! -- algorithm matches HYDRAD's COpticallyThickIon::GetAbundData exactly.
  subroutine load_abundance_mixture(path_ab, path_amu, Abund_H, Abund_Ca, Abund_Mg)
    character(*), intent(in) :: path_ab, path_amu
    real(dp), intent(out) :: Abund_H, Abund_Ca, Abund_Mg

    integer, parameter :: MAXEL = 128
    integer :: iZ(MAXEL), n, i, iu, nmasses, iZm, j, ios
    real(dp) :: fAb(MAXEL), fm_g(MAXEL), fAbH, fSumAbund, logAb, amu_val
    real(dp), parameter :: AMU = 1.660538782e-24_dp
    character(len=200) :: line

    ! Read each "Z  log10Ab" line as a whole record and parse both values
    ! from it together -- list-directed READ discards the rest of a record
    ! once a statement's item list is satisfied, so issuing separate reads
    ! for Z and log10Ab in different statements desynchronises immediately
    ! (each subsequent read jumps to the next line rather than continuing
    ! the current one). Terminated by a lone "-1" line.
    open(newunit=iu, file=path_ab, status='old', action='read')
    n = 0
    do
      read(iu,'(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (line(1:2) == '-1') exit
      n = n + 1
      read(line,*) iZ(n), logAb
      if (n == 1) then
        fAbH = 10.0_dp**logAb
        fAb(1) = 1.0_dp
      else
        fAb(n) = 10.0_dp**logAb/fAbH
      end if
    end do
    close(iu)

    Abund_H  = fAb(1)
    Abund_Ca = 0.0_dp
    Abund_Mg = 0.0_dp
    do i = 1, n
      if (iZ(i) == 20) Abund_Ca = fAb(i)
      if (iZ(i) == 12) Abund_Mg = fAb(i)
    end do

    fSumAbund = sum(fAb(1:n))
    fAb(1:n) = fAb(1:n)/fSumAbund

    open(newunit=iu, file=path_amu, status='old', action='read')
    read(iu,*) nmasses
    j = 1
    do i = 1, nmasses
      read(iu,*) iZm, amu_val
      if (j <= n) then
        if (iZm == iZ(j)) then
          fm_g(j) = amu_val*AMU
          j = j + 1
        end if
      end if
    end do
    close(iu)

    fSumAbund = sum(fAb(1:n)*fm_g(1:n))
    N_H_per_gram_global = fAb(1)/fSumAbund
  end subroutine load_abundance_mixture

  !--------------------------------------------------------------
  subroutine load_species_table(species_prefix, base_dir, Abund, tab)
    character(*), intent(in) :: species_prefix, base_dir
    real(dp), intent(in) :: Abund
    type(species_table_t), intent(out) :: tab

    integer :: iu, i
    real(dp) :: logT_raw, log10Q_raw

    ! ionisation balance table (log10 T, fraction)
    open(newunit=iu, file=trim(base_dir)//'/balances/'//trim(species_prefix)//'.bal', &
         status='old', action='read')
    read(iu,*) tab%n_bal
    allocate(tab%logT_bal(tab%n_bal), tab%frac_bal(tab%n_bal))
    do i = 1, tab%n_bal
      read(iu,*) tab%logT_bal(i), tab%frac_bal(i)
    end do
    close(iu)

    ! emissivity table (log10 T, log10 Q); precombine with Abund*IonFrac(T),
    ! matching HYDRAD's ppEmiss precomputation exactly
    open(newunit=iu, file=trim(base_dir)//'/emissivities/'//trim(species_prefix)//'.em', &
         status='old', action='read')
    read(iu,*) tab%n_em
    allocate(tab%logT_em(tab%n_em), tab%Qcombined_em(tab%n_em))
    do i = 1, tab%n_em
      read(iu,*) logT_raw, log10Q_raw
      tab%logT_em(i) = logT_raw
      tab%Qcombined_em(i) = 10.0_dp**log10Q_raw * Abund * &
           interp_linear(tab%logT_bal, tab%frac_bal, tab%n_bal, logT_raw, clamp=.true.)
    end do
    close(iu)

    ! escape probability table (depth variable, probability)
    open(newunit=iu, file=trim(base_dir)//'/escape_probabilities/'//trim(species_prefix)//'.esc', &
         status='old', action='read')
    read(iu,*) tab%n_esc
    allocate(tab%x_esc(tab%n_esc), tab%prob_esc(tab%n_esc))
    do i = 1, tab%n_esc
      read(iu,*) tab%x_esc(i), tab%prob_esc(i)
    end do
    close(iu)
  end subroutine load_species_table

  !--------------------------------------------------------------
  ! Returns Q(T)*Abund*IonFrac(T) (the combined emissivity factor) and,
  ! separately, the raw ionisation fraction (needed by the caller to build
  ! the neutral-H column density for the H I depth variable).
  subroutine species_emiss_ionfrac(tab, log10T, Qcombined, ionfrac)
    type(species_table_t), intent(in) :: tab
    real(dp), intent(in) :: log10T
    real(dp), intent(out) :: Qcombined, ionfrac
    Qcombined = interp_linear(tab%logT_em, tab%Qcombined_em, tab%n_em, log10T, clamp=.true.)
    ionfrac   = interp_linear(tab%logT_bal, tab%frac_bal, tab%n_bal, log10T, clamp=.true.)
    ionfrac   = max(0.0_dp, min(1.0_dp, ionfrac))
  end subroutine species_emiss_ionfrac

  function species_escape_prob(tab, log10depth) result(prob)
    type(species_table_t), intent(in) :: tab
    real(dp), intent(in) :: log10depth
    real(dp) :: prob
    ! HYDRAD extrapolates the escape probability linearly beyond the
    ! table's own range, then clamps the result to [0,1] -- reproduced here.
    prob = interp_linear(tab%x_esc, tab%prob_esc, tab%n_esc, log10depth, clamp=.false.)
    prob = max(0.0_dp, min(1.0_dp, prob))
  end function species_escape_prob

  !--------------------------------------------------------------
  ! CHIANTI coronal (optically thin) Lambda(T) table, same file/format as
  ! the original tr_equilibrium code's load_radiative_table.
  subroutine load_coronal_table(path, tab)
    character(*), intent(in) :: path
    type(coronal_table_t), intent(out) :: tab
    integer :: iu, i
    character(len=200) :: dummy
    real(dp) :: logT, lambda_raw
    real(dp), allocatable :: logT_arr(:)

    open(newunit=iu, file=path, status='old', action='read')
    do i = 1, 5
      read(iu,'(A)') dummy
    end do
    tab%n = 600
    allocate(tab%loglambda(tab%n), logT_arr(tab%n))
    do i = 1, tab%n
      read(iu,*) logT, lambda_raw
      logT_arr(i) = logT
      tab%loglambda(i) = log10(lambda_raw*1.0e-22_dp)
    end do
    close(iu)
    tab%logTmin = logT_arr(1)
    tab%dlogT = logT_arr(2) - logT_arr(1)
  end subroutine load_coronal_table

  function coronal_lambda(tab, log10T) result(lam)
    type(coronal_table_t), intent(in) :: tab
    real(dp), intent(in) :: log10T
    real(dp) :: lam, pos, frac
    integer :: idx
    pos = (log10T - tab%logTmin)/tab%dlogT + 1.0_dp
    idx = max(1, min(int(pos), tab%n - 1))
    frac = max(0.0_dp, min(1.0_dp, pos - real(idx, dp)))
    lam = 10.0_dp**(tab%loglambda(idx) + frac*(tab%loglambda(idx+1) - tab%loglambda(idx)))
  end function coronal_lambda

  !--------------------------------------------------------------
  subroutine load_hconduct_table(path, tab)
    character(*), intent(in) :: path
    type(hconduct_table_t), intent(out) :: tab
    integer :: iu, i
    open(newunit=iu, file=path, status='old', action='read')
    read(iu,*) tab%n
    allocate(tab%logT(tab%n), tab%log10kappa(tab%n))
    do i = 1, tab%n
      read(iu,*) tab%logT(i), tab%log10kappa(i)
    end do
    close(iu)
  end subroutine load_hconduct_table

  ! Extra (neutral-hydrogen) conductivity coefficient, to be ADDED to the
  ! standard Spitzer kappa_0 before the shared T^2.5 scaling is applied
  ! (matches HYDRAD: Kappa_B = (Kappa_Spitzer + Getkappa_0(T)) * T^2.5).
  ! Returns 0 above the table's own range -- extrapolating a neutral-H
  ! conduction enhancement where hydrogen is no longer significantly
  ! neutral would not be physical.
  function extra_h_kappa0(tab, log10T) result(kap0)
    type(hconduct_table_t), intent(in) :: tab
    real(dp), intent(in) :: log10T
    real(dp) :: kap0
    if (log10T > tab%logT(tab%n)) then
      kap0 = 0.0_dp
      return
    end if
    kap0 = 10.0_dp**interp_linear(tab%logT, tab%log10kappa, tab%n, log10T, clamp=.true.)
  end function extra_h_kappa0

  !--------------------------------------------------------------
  ! Smooth replacement for HYDRAD's hard switch at T=2.4e4 K: weight -> 1
  ! (fully chromospheric/Carlsson-Leenaarts) for T well below T_switch,
  ! -> 0 (fully coronal/optically-thin) well above it. Width in K.
  function blend_weight(T, T_switch, width) result(w)
    real(dp), intent(in) :: T, T_switch, width
    real(dp) :: w
    w = 0.5_dp*(1.0_dp - tanh((T - T_switch)/width))
  end function blend_weight

  !--------------------------------------------------------------
  ! Plain linear interpolation on a (not necessarily uniform) grid via
  ! bisection-free linear search -- these tables are all a few tens of
  ! points, called a handful of times per cell per step, so simplicity
  ! wins over building a fancier search for a negligible cost difference.
  ! clamp=.true. clamps to the edge VALUE outside the table range (matches
  ! HYDRAD's GetEmiss); clamp=.false. linearly extrapolates using the two
  ! nearest points instead (matches HYDRAD's GetEscProb).
  function interp_linear(x, y, n, xq, clamp) result(yq)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), y(n), xq
    logical, intent(in) :: clamp
    real(dp) :: yq, frac
    integer :: i

    if (clamp) then
      if (xq <= x(1)) then
        yq = y(1); return
      else if (xq >= x(n)) then
        yq = y(n); return
      end if
    end if

    i = 1
    do while (i < n-1 .and. x(i+1) < xq)
      i = i + 1
    end do

    frac = (xq - x(i))/(x(i+1) - x(i))
    yq = y(i) + frac*(y(i+1) - y(i))
  end function interp_linear

end module radiative_loss_mod

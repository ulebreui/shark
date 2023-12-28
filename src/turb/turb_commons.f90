module turb_commons
  use precision
  real(dp), dimension(:,:), allocatable :: f_turb
  real(dp), dimension(:,:), allocatable :: f_turb_old
  integer, parameter :: n_modes = 6

  real(dp), dimension(1:n_modes) :: f_k = 1.0d0
  real(dp), dimension(1:n_modes) :: phase_turb = 1.0d0

  real(dp), dimension(1:n_modes) :: f_k_old = 1.0d0
  real(dp) :: turb_turnover = 1.0d0
  real(dp) :: turb_correl   = 1.0d0
  real(dp) :: v_rms   = 1.0d5 ! cm /s


  integer, parameter :: turb_seed = 86456

  !Turbulent commons from RAMSES
  real(dp), allocatable     :: fturb(:,:)           ! Turbulent force
  complex(cdp), allocatable :: turb_last(:,:,:,:)   ! Turbulent spectrum at time = t
  complex(cdp), allocatable :: turb_next(:,:,:,:)   ! Turbulent spectrum at time = t + dt
  real(dp), allocatable     :: afield_last(:,:,:,:) ! Forcing field at time = t
  real(dp), allocatable     :: afield_next(:,:,:,:) ! Forcing field at time = t + dt
  real(dp), allocatable     :: afield_now(:,:,:,:)  ! Forcing field now
  real(dp), allocatable     :: power_spec(:,:,:)    ! Power spectrum of turbulence

  real(dp)    :: sol_frac         ! Solenoidal fraction
  real(dp)    :: turb_last_time   ! Time of old turbulent field
  real(dp)    :: turb_next_time   ! Time of next turbulent field
  real(dp)    :: turb_dt          ! Turbulent velocity evolution timestep
  real(dp)    :: turb_decay_frac  ! Decay fraction per dt
  real(dp)    :: turb_space(1:3)  ! Grid spacing
  real(dp)    :: turb_norm        ! Normalizing constant from combination
                                       ! of Ornstein-Uhlenbeck process, initial
                                       ! power spectrum and projection

  character(len=256) :: turb_file_last    ! filename for 'last' field
  character(len=256) :: turb_file_next    ! filename for 'last' field
  character(len=256) :: turb_file_header  ! filename for header file

  integer(ILP)       :: kiss64_state(1:4) ! State variables for KISS64 PRNG (Marsaglia)

end module turb_commons

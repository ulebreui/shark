module setup_parameters
  use precision
  
  !Cloud & Gas properties 
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0
  real(dp),dimension(1:NDUST) :: dust2gas_ratio = 1.0d-2 !TODO to check/change when considering a distribution
  real(dp),dimension(1:NDUST) :: St_0 = 1.0d0
  real(dp) :: kx_wave_B   = 1.0d0
  real(dp) :: ky_wave  = 1.0d0
  real(dp) :: kz_wave   = 1.0d0

  real(dp) :: kx_wave_d   = 1.0d0
  real(dp) :: ky_wave_d  = 1.0d0
  real(dp) :: kz_wave_d   = 1.0d0

  real(dp) :: beta_0  = 1.0d0
  real(dp) :: delta_B   = 1.0d0
  real(dp) :: vy_0  = 0d0
  real(dp) :: vz_0  = 0d0


  real(dp) :: vdx_0   = 0.0d0
  real(dp) :: vdy_0   = 0.0d0
  real(dp) :: vdz_0   = 0.0d0
  real(dp) :: cs_0   = 1.0d-1


  real(dp) :: rho_0 = 1.0d0
  real(dp) :: delta_rho   = 0.0d0
  real(dp) :: delta_rho_d   = 0.0d0

  logical :: decaying_turb_compressive = .false.   ! Initialize compressive modes for decaying turbulence
  logical :: decaying_turb_solenoidal = .false.   ! Initialize solenoidal modes for decaying turbulence
  character (len=300) :: decay_turb_random_path = "test"



  integer :: nb_turb_modes   = 1

  integer , dimension(:)  , allocatable    :: k_turb
  real(dp), dimension(:)  , allocatable    :: vx_turb
  real(dp), dimension(:)  , allocatable    :: vy_turb
  real(dp), dimension(:)  , allocatable    :: vz_turb
  real(dp), dimension(:)  , allocatable    :: phix_turb
  real(dp), dimension(:)  , allocatable    :: phiy_turb
  real(dp), dimension(:)  , allocatable    :: phiz_turb









end module setup_parameters

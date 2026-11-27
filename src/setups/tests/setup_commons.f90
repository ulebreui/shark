module setup_parameters
  use precision
  implicit none

  ! Standard 1D validation tests (Cartesian, GEOM=0, NY=1), selected with test_case:
  !   'sod'        Sod shock tube (adiabatic, NDUST=0)
  !   'soundwave'  isothermal linear sound wave (NDUST=0)
  !   'dustybox'   uniform drag relaxation (Laibe & Price 2011)
  !   'dustywave'  linear gas+dust acoustic wave (Laibe & Price 2011)
  !   'dustyshock' stationary isothermal dusty shock (Loren-Aguilar & Bate 2014), constant drag coefficient
  !   'smoluchowski' coagulation with constant/additive kernel (kernel_type in DUST_PARAMS), uniform medium

  real(dp), parameter :: mu_gas = 2.31d0      ! needed by the dust growth module

  character(len=20) :: test_case = 'sod'
  real(dp) :: box_l   = 1.0d0                 ! domain length
  real(dp) :: box_l_y = 1.0d0                 ! (unused in 1D)
  real(dp) :: cs0     = 1.0d0                 ! isothermal sound speed
  real(dp) :: amp     = 1.0d-4                ! amplitude of the linear waves
  real(dp) :: nwave   = 1.0d0                 ! number of wavelengths in the box

  ! Sod
  real(dp) :: sod_rho_l = 1.0d0, sod_p_l = 1.0d0, sod_rho_r = 0.125d0, sod_p_r = 0.1d0
  real(dp) :: x_disc = 0.5d0                  ! position of the initial discontinuity (sod, dustyshock)

  ! Dust (dustybox, dustywave): dust-to-gas ratio, stopping time, initial velocity of each species
  real(dp), dimension(1:NDUST) :: eps_dust = 1.0d0
  real(dp), dimension(1:NDUST) :: ts_dust  = 0.1d0
  real(dp), dimension(1:NDUST) :: vd0      = 0.0d0
  real(dp) :: vg0 = 1.0d0                     ! dustybox: initial gas velocity

  ! Dusty shock: upstream velocity, dust-to-gas ratio, drag coefficient K (t_s = rho_d/K)
  real(dp) :: v_shock   = 2.0d0
  real(dp) :: eps_shock = 1.0d0
  real(dp) :: K_drag    = 1.0d0

  real(dp), dimension(:,:,:), allocatable :: uprim_condinit   ! initial state (fixed boundaries)

end module setup_parameters

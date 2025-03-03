module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l       = 6.0d0 ! in R0
  real(dp) :: r0          = 1540  !in AU
  real(dp) :: rout        = 10000 !in AU
  real(dp) :: r_cut        = 100 !in AU

  real(dp) :: smooth_r    = 6.0  ! softening lenght in AU
  real(dp) :: rsink       = 10.0  ! Sink lenght in AU

  real(dp) :: l_soft    =  4.0  ! softening lenght in dx

  real(dp) :: sigma_0     = 0.37 ! g/cm^2
  real(dp) :: sigma_crit  = 500 ! g/cm^2
  real(dp) :: rho_crit    = 1d-13 ! g/cm^3

  real(dp) :: sigma_sink  = 1500 ! g/cm^2
  real(dp) :: rho_pert    = 0.0d0! g/cm^2

  real(dp) :: Mstar    = 0.0d0           ! in solar mass
  real(dp) :: Mtot     = 0.0d0          ! in solar mass
  real(dp) :: M_cloud  = 1.0d0           ! in solar mass
  real(dp) :: Vol_tot  = 0.0d0          ! in code units

  real(dp) :: Omega_0     = 1e-14   ! in s-1 corresponds to 1.45 km/s/pc
  real(dp) :: T0          = 10.d0    ! in K

  real(dp) :: alpha_cloud = 0.25
  real(dp) :: beta_cloud  = 0.01

  real(dp),dimension(:,:,:),allocatable :: uprim_condinit 

end module setup_parameters
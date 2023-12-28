module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0
  real(dp) :: box_l_y = 1.0d0

  real(dp) :: M_cloud     = 1.0d0
  real(dp) :: alpha_cloud = 0.25d0
  real(dp) :: T0_cloud    = 10.0d0
  real(dp) :: rho_max_sim = 1d-10
  logical  :: more_outputs =.false.
  real(dp) :: order_mag
  logical  :: single_size = .false.



end module setup_parameters

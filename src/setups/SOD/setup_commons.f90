module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0
  real(dp) :: box_l_y = 1.0d0

  real(dp) :: rho_l = 1.0d0
  real(dp) :: rho_r = 0.125d0
  real(dp) :: P_l = 1.0d0
  real(dp) :: P_r = 0.1d0
  real(dp) :: v_l = 0.d0
  real(dp) :: v_r = 0.0d0
  integer  :: direction_shock=0
end module setup_parameters

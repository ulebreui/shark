module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0
  real(dp) :: box_l_y = 1.0d0

  real(dp) :: rho_dense = 1.0d0
  real(dp) :: rho_diffuse = 0.8d0
  real(dp) :: P_dense = 1.0d0
  real(dp) :: P_diffuse = 1.0d0

  real(dp) :: T_dense   = 10.0d0 !in Kelvins
  real(dp) :: T_diffuse = 10.0d0
  real(dp) :: v_dense   = 1.0d0
  real(dp) :: v_diffuse = -1.0d0
  real(dp) :: Mach = 0.5d0
  real(dp) :: vy0 = 0.1d0
  real(dp) :: kx  = 10.d0
end module setup_parameters

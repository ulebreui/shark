module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0
  real(dp) :: box_l_y = 1.0d0

  real(dp) :: rho_init    = 1.0d0
  real(dp) :: Omega_shear = 1.0d0
  real(dp) :: Q_shear   = 1.5d0
  real(dp) :: HoverR    = 0.1d0
  real(dp) :: eta_stream= 0.005d0
  real(dp) :: rad0= 1.0d0

end module setup_parameters

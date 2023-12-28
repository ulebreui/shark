module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0
  real(dp) :: box_l_y = 1.0d0

  real(dp) :: rho_ext= 1.0d0
  real(dp) :: P_cen = 1.0
  real(dp) :: P_ext = 1d-5
  real(dp) :: r_expl=0.1
end module setup_parameters

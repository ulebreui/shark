module setup_parameters
  use precision
  
  !Cloud & Gas properties 
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0
  real(dp),dimension(1:NDUST) :: dust2gas_ratio = 1.0d0 !TODO to check/change when considering a distribution/coagulation
  real(dp),dimension(1:NDUST) :: K = 1.0d0
  real(dp) :: kx_wave   = 1.0d0




  real(dp) :: vx_0   = 1.0d-2
  real(dp) :: cs_0   = 1.0d0


  real(dp) :: rho_0 = 1.0d0
  real(dp) :: delta_rho   = 1.0d-2
  real(dp) :: delta_rho_d   = 1.0d-2


end module setup_parameters

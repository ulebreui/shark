module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 4.0d0 ! in units of outer disk
  real(dp) :: box_l_y = 4.0d0 ! in units of outer disk

  real(dp) :: sigma_R0    = 1700.0d0 ! in g / cm ^2
  real(dp) :: R0_disk     = 1.0d0    ! in AU
  real(dp) :: Mstar       = 1.0d0    ! in solar mass
  real(dp) :: Mplanet     = 0.0d0    ! in solar mass
  real(dp) :: Rplanet     = 10.0d0    ! in au
  real(dp) :: buffer_size  = 0.1d0   ! in r_out

  real(dp) :: R_out_disk  = 25.0d0  ! outer radii
  real(dp) :: HoverR      = 0.05    ! H/R
  real(dp) :: smooth_r    = 1d0   ! softening lenght in dx





end module setup_parameters

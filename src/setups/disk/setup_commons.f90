module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 400.0d0 ! in au
  real(dp) :: box_l_y = 400.0d0 ! in au

  real(dp) :: sigma_R0    = 1700.0d0 ! in g / cm ^2
  real(dp) :: R0_disk     = 1.0d0    ! in AU
  real(dp) :: Mstar       = 1.0d0    ! in solar mass
  real(dp) :: Mplanet     = 0.0d0    ! in solar mass
  real(dp) :: Rplanet     = 20.0d0    ! in au

  real(dp) :: frac_R_disk  = 0.1d0  ! outer radii in box size units
  real(dp) :: buffer_in    = 0.1d-2 ! inner buffer
  real(dp) :: buffer_out   = 0.25d0  ! outer buffer radii in box size units

  real(dp) :: HoverR      = 0.05   ! H/R
  real(dp) :: smooth_r    = 0.5  ! softening lenght
  real(dp) :: alpha_visc  = 1d-2
  logical  :: accretion   = .false.
  real(dp) :: M_acc       = 1d-6 ! Solar mass per yer



end module setup_parameters
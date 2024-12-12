module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 400.0d0 ! in au
  real(dp) :: box_l_y = 400.0d0 ! in au
  real(dp) :: r_relax = 3.0d0
  real(dp) :: n_rel   = 8.0d0
  real(dp) :: sigma_R0    = 1700.0d0 ! in g / cm ^2
  real(dp) :: R0_disk     = 1.0d0    ! in AU
  real(dp) :: R_out_disk  = 40.d0
  real(dp) :: Mstar       = 1.0d0    ! in solar mass
  real(dp) :: Mplanet     = 0.0d0    ! in solar mass
  real(dp) :: Rplanet     = 20.0d0    ! in au
  real(dp) :: a_smooth_pl = 0.1
  real(dp) :: decrease_density = 1d-2
  real(dp) :: HoverR      = 0.1   ! H/R
  real(dp) :: smooth_r    = 1.0  ! softening lenght
  real(dp) :: alpha_visc  = 1d-2
  logical  :: accretion   = .false.
  logical  :: test_planet_disk = .false.
  real(dp) :: M_acc       = 1d-7! Solar mass per yer
  real(dp) :: phi_mom     = 2.0d0
  real(dp) :: L_out       = 1.0d0
  real(dp),dimension(:,:,:),allocatable :: uprim_condinit 

end module setup_parameters
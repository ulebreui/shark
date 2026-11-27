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
  real(dp) :: a_smooth_pl = 0.6
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
  
  logical  :: damping       = .true.   ! de Val-Borro (2006) wave-damping zones
  real(dp) :: damp_in_frac  = 1.25d0   ! inner zone: r < damp_in_frac  * r_in
  real(dp) :: damp_out_frac = 0.84d0   ! outer zone: r > damp_out_frac * r_out
  real(dp) :: damp_tau      = 0.1d0    ! damping time, in orbital periods at the boundary
  real(dp), dimension(:), allocatable :: damp_rate   ! R(r)/tau for each ring

  real(dp) :: bump_amp = 0.0d0   ! test: amplitude of a Gaussian dust bump (0 = no bump)
  real(dp) :: bump_r   = 10.0d0  ! its radius (au)
  real(dp) :: bump_phi = 3.0d0   ! its azimuth (rad)
  real(dp) :: bump_w   = 2.0d0   ! its width (au)
  integer, parameter :: nplanet_max = 10
  integer  :: nplanet = 0                                     ! number of planets
  real(dp), dimension(nplanet_max) :: Mplanet      = 0.0d0   ! masses (Msun)
  real(dp), dimension(nplanet_max) :: Rplanet      = 20.0d0  ! orbital radii (au)
  real(dp), dimension(nplanet_max) :: phi_planet   = 0.0d0   ! initial azimuths (rad)
  real(dp), dimension(nplanet_max) :: nramp_planet = 5.0d0   ! mass ramp-up duration (planet orbits)
  logical  :: indirect_term = .true.  ! indirect term of the planets (star-centred frame)
  integer  :: freq_torque = 50       ! torques written every freq_torque steps
end module setup_parameters
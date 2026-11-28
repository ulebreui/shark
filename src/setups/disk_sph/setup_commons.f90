module setup_parameters
  use precision
  implicit none

  ! Axisymmetric disk on a spherical (r,theta) grid (GEOM=3)
  ! radii = r, phi = theta, ivx = v_r, ivy = v_theta, ivz = v_phi

  real(dp), parameter :: mu_gas = 2.31d0      ! needed by the dust growth module

  ! Grid
  real(dp) :: r_in   = 1.0d0                  ! inner radius [au]
  real(dp) :: r_out  = 10.0d0                 ! outer radius [au]
  real(dp) :: th_min = 1.2707963267948966d0   ! pi/2 - 0.3
  real(dp) :: th_max = 1.9707963267948966d0   ! pi/2 + 0.3

  ! Star and disk (TODO: hydrostatic IC)
  real(dp) :: Mstar  = 1.0d0                  ! [Msun]
  real(dp) :: HoverR = 0.05d0                 ! H/R at r0
  real(dp) :: r0      = 1.0d0       ! reference radius [au]
  real(dp) :: rho0    = 1.0d-10     ! midplane density at r0 [g/cm^3]
  real(dp) :: p_dens  = -1.5d0      ! rho_mid ~ R^p
  real(dp) :: q_temp  = -1.0d0      ! cs^2 ~ R^q  (q=-1: constant H/R)
  real(dp), dimension(:,:,:), allocatable :: uprim_condinit   ! (ivar,ix,iy) initial state, used by the fixed boundaries
  
  ! de Val-Borro (2006) radial damping zones, relaxation towards the IC
  logical  :: damping       = .true.
  real(dp) :: damp_in_frac  = 1.25d0     ! inner zone: r < damp_in_frac *r_in
  real(dp) :: damp_out_frac = 0.84d0     ! outer zone: r > damp_out_frac*r_out
  real(dp) :: damp_tau      = 0.1d0      ! damping time in orbital periods at the boundary
  real(dp), dimension(:,:), allocatable :: damp_rate     ! 1/tau * R(r)^2, precomputed
end module setup_parameters

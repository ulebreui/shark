module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp), parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0
  real(dp) :: box_l_y = 1.0d0

  real(dp) :: M_cloud      = 1.0d0
  real(dp) :: alpha_cloud  = 0.25d0
  real(dp) :: T0_cloud     = 10.0d0
  real(dp) :: rho_max_sim  = 1d-10
  real(dp) :: s_max  = 1.0d-6

  logical  :: more_outputs = .false.
  real(dp) :: order_mag
  logical  :: single_size  = .false.
  integer  :: ntimes_part  = 1
  integer  :: ilist=1
  real(dp), dimension(:), allocatable   :: times_particles
  real(dp), dimension(:), allocatable :: rho_particles

end module setup_parameters

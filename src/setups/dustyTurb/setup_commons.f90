module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0
  real(dp) :: box_l_y = 1.0d0

  real(dp) :: rho_init    = 1.0d0
  real(dp) :: cs0 = 1.0d0
  real(dp) :: Mach = 1.0d0
  real(dp) :: tcross = 1.0d0
  real(dp) :: turnover_time=1.0d0
  real(dp), dimension(:), allocatable :: f_hat
  real(dp), dimension(:), allocatable :: phase_forcing

  real(dp), dimension(:), allocatable :: f_hat_old
  integer :: numb_modes = 20
  real(dp), dimension(1:NDUST):: Stokes_species   = 0.2d0
  real(dp), dimension(1:NDUST):: dust2gas_species = 1.0d0

  real(dp) :: Stokes_min  = 1d-3
  real(dp) :: Stokes_max  = 1.0d0
  real(dp) :: Stokes_cut  = 0.01d0
  real(dp) :: theta_dust  = 1d-11 ! Ratio between disk density and dust grain density

  logical  :: stokes_distrib = .false. ! Use a MRN like distribution

end module setup_parameters

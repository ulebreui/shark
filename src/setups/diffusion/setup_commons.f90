module setup_parameters
  use precision
  use phys_const
  

  real(dp) :: D_test_value=1.0d0
  logical :: diffusion_test = .false.

  !Cloud & Gas properties 
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0

  real(dp) :: cs_0   = 1.0d-1
  real(dp) :: rho_0 = 1.0d0


  real(dp), dimension(:), allocatable   :: D_test
  real(dp), dimension(:), allocatable   :: By_old
  real(dp), dimension(:), allocatable   :: Bz_old













end module setup_parameters

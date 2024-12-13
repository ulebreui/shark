module setup_parameters
  use precision
  use phys_const
  

  real(dp) :: D_test_value=1.0d0
  logical :: diffusion_test = .false.
  integer :: m = 2 !Characterizes the level of non-linearity
  real(dp) :: Cte = 1.0d0 !Related to initial value for BP solution
  real(dp) :: t_init = 1.0d-2 !initial time to compare with analytical sol.
  real(dp) :: x_wid= 0.05 ! Control width of initial condition.
  real(dp) :: By_init= 1.0d0 ! Height initial condition.


  !Cloud & Gas properties 
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0

  real(dp) :: cs_0   = 1.0d-1
  real(dp) :: rho_0 = 1.0d0
















end module setup_parameters

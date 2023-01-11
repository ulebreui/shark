module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp) :: M_cloud     = 1.0d0     ! in sun masses
  real(dp) :: alpha_cloud = 0.25d0    ! Virial parameter
  real(dp) :: T0_cloud    = 10.0d0    ! Isothermal temperature of the cloud at low density 
  real(dp) :: gamma       = 1.66667   ! Adiabatic index
  logical  :: ramses_baro = .true.
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: order_mag   = 1d-24     ! Low number to find the intial min density  (for outputing)
  character(len=60) :: grid_type='log'
  real(dp) :: rho_max_sim=1d-10
  logical  :: more_outputs=.false. ! To plot more outputs
end module setup_parameters

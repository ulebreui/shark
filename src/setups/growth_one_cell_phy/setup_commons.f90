module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp) :: M_cloud     = 1.0d0               ! in sun masses
  real(dp) :: alpha_cloud = 0.25d0              ! Virial parameter
  real(dp) :: T0_cloud    = 10.0d0              ! Isothermal temperature of the cloud at low density
  real(dp) :: nh_cell     = 1d5                 ! In cm^{-3}
  real(dp) :: gamma       = 1.66667             ! Adiabatic index
  logical :: isotherm     = .false.
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: order_mag   = 1d-24               ! Low number to find the intial min density  (for outputing)
  real(dp),dimension(1:100) :: list_times=0.!(0,10,100,1000,1d4,1d5,1d6)
  real(dp) :: dt_fix=0.
  integer   :: dt_strat = 0 ! 0 only 2 outputs, 1 list of times, 2 fixed frequency
  integer :: ilist=1
end module setup_parameters

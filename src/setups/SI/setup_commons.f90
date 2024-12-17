module setup_parameters
  use precision
  
  !Cloud & Gas properties
  real(dp),parameter :: mu_gas      = 2.31d0    ! Mean molecular weight
  real(dp) :: box_l   = 1.0d0
  real(dp) :: box_l_y = 1.0d0

  real(dp) :: rho_init    = 1.0d0
  real(dp) :: Omega_shear = 1.0d0
  real(dp) :: Q_shear     = 1.5d0
  real(dp) :: HoverR      = 0.1d0
  real(dp) :: eta_stream  = 0.005d0
  real(dp) :: rad0        = 1.0d0
  real(dp) :: mag_pert        = 2.0d-2 !fraction of the soundspeed to be used for initial density perturbation


  real(dp), dimension(1:NDUST):: Stokes_species   = 0.2d0
  real(dp), dimension(1:NDUST):: dust2gas_species = 1.0d0

  real(dp) :: Stokes_min  = 1d-3
  real(dp) :: Stokes_max  = 1.0d0
  real(dp) :: Stokes_cut  = 0.01d0
  real(dp) :: theta_dust  = 1d-11 ! Ratio between disk density and dust grain density

  logical  :: stokes_distrib = .false. ! Use a MRN like distribution

end module setup_parameters

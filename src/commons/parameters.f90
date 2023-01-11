module parameters
  use precision
  use phys_const
  use setup_parameters
#if NDUST>0  
  use dust_parameters
#endif  
  implicit none

  integer,parameter   ::  nghost         = NGHOST
  integer,parameter   ::  ncells         = NHBINS + 2*nghost !Number of cells + ghost cells at each sides
 
  ! Computational domain boundaries
  integer,parameter   ::  first_active   = nghost + 1
  integer,parameter   ::  last_active    = ncells-nghost
  
  ! Ghost domain boundaries
  integer,parameter   ::  inner_bound   = nghost 
  integer,parameter   ::  outer_bound   = ncells-nghost + 1

  ! Variables
  integer,parameter   ::  Ndust       = NDUST  !Number of dust species
  integer,parameter   ::  Nmhd        = MHD    !Number of mhd variables
  integer,parameter   ::  NENERGY     = ENERGY
  integer,parameter   ::  nvar        = 2 + Ndust * 2 + MHD +ENERGY!Number of variables


  logical             ::  static            = .false.
  logical             ::  zero_flux_in      = .false. ! Zero flux at physical boundaries
  integer             ::  freq_out          = 1000    ! Output frequency
  logical             ::  stop_at_first_core= .false. ! To stop the integration at the first core formation
  real(dp)            ::  NR_CLOUD = 4.0d0 ! Size of the box (in cloud radius)
  real(dp)            ::  Nff      = 1.3d0 ! Number of free-fall timescales to integrate (if stop_at_first_core)
  real(dp)            ::  rin = 1.0d0  ! Inner radius boundary
  real(dp)            ::  CFL = 0.5d0  ! CFL constant
  real(dp)            ::  t21 = 0.0d0
  real(dp)            ::  t32 = 0.0d0
  real(dp)            ::  t43 = 0.0d0
  real(dp)            ::  t54 = 0.0d0
  real(dp)            ::  t65 = 0.0d0
  real(dp)            ::  t76 = 0.0d0

  integer :: nrestart  =0 ! For restart
  integer :: restarting=0 !
 end module parameters

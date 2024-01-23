module parameters
  use precision
  use phys_const
  use setup_parameters
#if NDUST>0  
  use dust_parameters
#endif  
  implicit none

  integer,parameter   ::  nghost         = NGHOST
  integer,parameter   ::  nx             = NX
  integer,parameter   ::  ny             = NY
#if NY==1
  integer,parameter   ::  ncells        = NX + 2*nghost !Number of cells + ghost cells at each sides
  integer,parameter   ::  ndim          = 1
  integer,parameter   ::  ncells_active = NX
  integer,parameter   ::  nx_max        = NX + 2*nghost
  integer,parameter   ::  ny_max        = 1
  integer,parameter   ::  ny_max_minus_1= 1
#else
  integer,parameter   ::  ncells        = (NX+2*nghost)*(NY+2*nghost) !Number of cells + ghost cells at each sides
  integer,parameter   ::  ncells_active = NX*NY 
  integer,parameter   ::  ndim          = 2
  integer,parameter   ::  nx_max        = NX + 2*nghost
  integer,parameter   ::  ny_max        = NY + 2*nghost
  integer,parameter   ::  ny_max_minus_1= NY-1
#endif 
  ! Computational domain boundaries
  integer,parameter   ::  first_active     = nghost + 1
  integer,parameter   ::  last_active      = NX+2*nghost-nghost
#if NY>1
  integer,parameter   ::  first_active_y   = nghost + 1
  integer,parameter   ::  last_active_y    = NY+2*nghost-nghost
#else
  integer,parameter   ::  first_active_y   = 1
  integer,parameter   ::  last_active_y    = 1
#endif 

  ! Variables
  integer,parameter   ::  Ndust       = NDUST  !Number of dust species
  integer, parameter  ::  Nmhd        = MHD*3
  integer,parameter   ::  nvar        = 5 + Ndust * (4) + Nmhd !Number of variables

  real(dp),parameter            :: half = 0.5d0

  logical             ::  static            = .false.
  logical             ::  force_kick        = .false.
  logical             ::  viscosity         = .false.
  integer             ::  freq_out          = 1000    ! Output frequency
  real(dp)            ::  rin = 0.0d0  ! Inner radius boundary for cylindrical geometry
  real(dp)            ::  CFL = half  ! CFL constant
  real(dp)            ::  t21   = 0.0d0
  real(dp)            ::  t32   = 0.0d0
  real(dp)            ::  t43   = 0.0d0
  real(dp)            ::  t54   = 0.0d0
  real(dp)            ::  t65   = 0.0d0
  real(dp)            ::  t76   = 0.0d0
  real(dp)            ::  t87   = 0.0d0
  real(dp)            ::  t98   = 0.0d0
  real(dp)            ::  t109  = 0.0d0
  real(dp)            ::  t1110 = 0.0d0
  integer :: nrestart   = 0 ! For restart
  integer :: restarting = 0 !
 end module parameters

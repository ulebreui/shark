module dust_parameters
  use precision 
  implicit none

  !Flags
  logical             ::  drag                  = .true.      ! Dust growth is activated  
  logical             ::  growth                = .false.     ! Dust growth is activated
  logical             ::  growth_step           = .false.     ! Dust growth is activated with Stepinski solver
  logical             ::  fragmentation         = .false.     ! Fragmentation is activated
  logical             ::  drift_in_growth       = .false.     ! Drift velocity included in growth
  logical             ::  turb_in_growth        = .false.     ! Turb velocity included in growth /!\ this is similar to drift
  logical             ::  brownian_in_growth    = .false.     ! Brownian velocity included in growth
  logical             ::  dust_back_reaction    = .true.      ! Add the dust back-reaction
  real(dp)            ::  sticking_efficiency   = 1.0d0       ! Add the electrostatic barrier for dust growth

  real(dp)            ::  dtcontrol_growth  = -1.0d0
  character (len=60)  ::  dust_distribution = 'mrn'
  
  real(dp):: CFL_growth         = 0.1d0 ! CFL for smoluchowski
  real(dp):: eps_threshold      = 1d-10
  real(dp):: eps_threshold_frag = 1d-10

  real(dp):: rhodust_threshold  = 1d-27
  real(dp):: dust_ratio_min     = 1d-11
  integer :: kernel_type        = 0   ! 0 = physical, 1= constant, 2 = additive
  integer ::  frag_thre         = 0   ! 0 = NRJ, 1= velocity


  logical ::  dust_growth_disk = .false.
  logical             ::  SI = .false.     ! To compute growth within Streaming Instability setup

  !Dust distribution
  real(dp)::  smin            = 1d-7    ! Minimum possible dust size 
  real(dp)::  scutmin         = 5d-7    ! MRN - Minimum dust size of the initial distribution
  real(dp)::  scut            = 250d-7  ! MRN - Maximum dust size of the initial distribution
  real(dp)::  smax            = 1d-2    ! Maximum possible dust size
  real(dp)::  mrn             = 3.5d0   ! MRN power law slope 
  real(dp)::  aO_themis       = 1d-4    ! Themis - Average size
  real(dp)::  acut_themis     = 10d-7   ! Themis - End of monomers
  real(dp)::  awidthcut_themis= 50d-7   ! Themis - End of monomers width
  real(dp)::  themis_slope    = 5.0
  real(dp)::  sigma_themis    = 1.0d0   ! Themis - STDev
  real(dp)::  ice_mantle      = 8.7d-7  ! Ice mantle size
  real(dp)::  rhograin        = 2.7d0   ! Grain intrinsic density
  real(dp)::  dust2gas        = 0.01d0  ! Dust-to-gas ratio

  !Frag parameters (mostly from Ormel's paper)
  real(dp)::  Abr     = 2.8d3 ! Blum & Wurm 2000
  real(dp)::  ksicrit = 2d-7  ! Blum & Wurm 2000
  !real(dp)::  gamma_grains = 370.0d0 ! erg cm-2
  !real(dp)::  estar_grains = 3.7d10 ! dyn cm-2
  real(dp)::  gamma_grains = 25.0d0 ! erg cm-2
  real(dp)::  estar_grains = 2.8d11 ! dyn cm-2  

  real(dp)::  vfrag = 100 ! cm/s
  !Monomer properties
  real(dp):: size_mono   = 1d-5    ! 0.1 micron
  real(dp):: slope_mono  = 3.5d0   ! Index of power law monomer size distribution
  real(dp):: alpha_turb  = 1.5d0
  real(dp):: sminstep    = 1d-5
end module dust_parameters

module dust_commons
  use precision  
  implicit none

 ! Dust related arrays

  real(dp), dimension(:)  , allocatable    :: aplus
  real(dp), dimension(:)  , allocatable    :: aminus
  real(dp), dimension(:)  , allocatable    :: mplus
  real(dp), dimension(:)  , allocatable    :: mminus
  real(dp), dimension(:,:), allocatable    :: epsilondust
  real(dp), dimension(:)  , allocatable    :: sdust
  real(dp), dimension(:)  , allocatable    :: mdust
  real(dp), dimension(:,:,:), allocatable  :: force_dust_x
  real(dp), dimension(:,:,:), allocatable  :: force_dust_y
  real(dp), dimension(:,:,:), allocatable  :: force_dust_z

  real(dp), dimension(:,:,:), allocatable    :: tstop
  real(dp), dimension(:,:,:), allocatable    :: tcoag
  real(dp), dimension(:,:,:), allocatable    :: St



  ! Indices of the variables
  integer,  dimension(:),  allocatable     :: irhod
  integer,  dimension(:),  allocatable     :: ivdx
  integer,  dimension(:),  allocatable     :: ivdy
  integer,  dimension(:),  allocatable     :: ivdz ! Z - velocity component
  integer,  dimension(:,:),  allocatable   :: idust_pscal
  real(dp) :: dt_cfl_dust=1d140


end module dust_commons

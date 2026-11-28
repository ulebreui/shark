module commons
  use precision
  use hydro_commons
#if NDUST>0
  use dust_commons
#endif
  use OMP_LIB
  implicit none

  !Grid
  real(dp), dimension(:,:,:), allocatable :: dx
  real(dp), dimension(:,:,:), allocatable :: position
  real(dp), dimension(:,:), allocatable   :: radii


#if GEOM>=2
  ! GEOM=2: azimuthal angle (0..2pi), GEOM=3: polar angle theta
  real(dp), dimension(:,:), allocatable     :: phi
  real(dp), dimension(:), allocatable       :: vbar   ! FARGO mean azimuthal velocity per ring
#endif
#if GEOM==2
  real(dp), dimension(:,:,:,:), allocatable :: distance
#endif

  real(dp), dimension(:,:), allocatable   :: vol
  real(dp), dimension(:,:,:), allocatable :: Surf
  
  real(dp):: time 
  real(dp):: tend
  real(dp):: dt

end module commons


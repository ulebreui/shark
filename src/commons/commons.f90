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

#if GEOM==2
  ! Disk face-on geometry /!\ Phi is the azimuthal angle here it ranges from 0 to 2pi
  real(dp), dimension(:,:), allocatable     :: phi
  real(dp), dimension(:,:,:,:), allocatable :: distance
#endif

  real(dp), dimension(:,:), allocatable   :: vol
  real(dp), dimension(:,:,:), allocatable :: Surf
  
  real(dp):: time 
  real(dp):: tend
  real(dp):: dt

end module commons


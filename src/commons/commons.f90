module commons
  use precision
  use hydro_commons
#if NDUST>0
  use dust_commons
#endif
  implicit none

  !Grid
  real(dp), dimension(:,:), allocatable :: dx
  real(dp), dimension(:,:), allocatable :: position
#if GEOM==1
  real(dp), dimension(:), allocatable   :: polar_radii
  real(dp), dimension(:), allocatable   :: theta

#endif
  real(dp), dimension(:,:), allocatable :: Surf
  real(dp), dimension(:), allocatable   :: vol
  real(dp), dimension(:), allocatable   :: active_cell
  
  
  real(dp):: time 
  real(dp):: tend
  real(dp):: dt

  
end module commons


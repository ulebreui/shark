module commons
  use precision
  use hydro_commons
#if NDUST>0
  use dust_commons
#endif
#if GRAVITY>0
  use gravity_commons
#endif
  use mhd_commons
  implicit none

  !Grid
  real(dp), dimension(:), allocatable :: r
  real(dp), dimension(:), allocatable :: r_c
  real(dp), dimension(:), allocatable :: vol_cell
  real(dp), dimension(:), allocatable :: dx
  real(dp), dimension(:), allocatable :: dx_l
  real(dp), dimension(:), allocatable :: dx_r
  real(dp), dimension(:), allocatable :: dx_c
  real(dp), dimension(:), allocatable :: dx_r_cell
  real(dp), dimension(:), allocatable :: dx_l_cell
  real(dp), dimension(:), allocatable :: Surf_p
  real(dp), dimension(:), allocatable :: Surf_m
  real(dp), dimension(:), allocatable :: dvol
  
  

  real(dp):: time 
  real(dp):: tend
  real(dp):: dt

  
end module commons


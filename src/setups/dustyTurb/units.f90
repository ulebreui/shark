module units
  use precision
  use parameters
  
  implicit none

    real(dp),parameter :: unit_l    = 1.0d0
    real(dp),parameter :: unit_d    = 1.0d0
    real(dp),parameter :: unit_dcol = 1.0d0
    real(dp),parameter :: unit_t    = 1.0d0
    real(dp),parameter :: unit_m    = 1.0d0
    real(dp),parameter :: unit_v    = 1.0d0
    real(dp),parameter :: unit_P    = 1.0d0
    real(dp),parameter :: unit_nH   = 1.0d0
    real(dp),parameter :: unit_B = dsqrt(4*pi*unit_d)*unit_v 

end module units

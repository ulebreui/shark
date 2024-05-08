module units
  use precision
  use parameters
  
  implicit none

    real(dp),parameter :: unit_l    = au
    real(dp),parameter :: unit_d    = msun/au**3
    real(dp),parameter :: unit_dcol = msun/au**2
    real(dp),parameter :: unit_t    = 1.0/sqrt(grav*msun/unit_l**3)
    real(dp),parameter :: unit_m    = msun
    real(dp),parameter :: unit_v    = unit_l/unit_t
    real(dp),parameter :: unit_P    = unit_v**2
    real(dp),parameter :: unit_B = dsqrt(4*pi*unit_d)*unit_v 
    real(dp),parameter :: unit_nH   = unit_d/unit_m
end module units

module units
  use precision
  use parameters
  
  implicit none

    real(dp),parameter :: unit_l = 1.0d0

    real(dp),parameter :: unit_d =1.0d0
    real(dp),parameter :: unit_P = 1.0d0
    real(dp),parameter :: unit_dcol=1.0d0
    real(dp),parameter :: unit_t = 1.0d0
    real(dp),parameter :: unit_m = unit_l**3*unit_d 

    real(dp),parameter :: unit_v = unit_l/unit_t
    real(dp),parameter :: unit_nH= 1.0d0



    !real(dp),parameter :: unit_t = sqrt(3*pi/(32*Grav*unit_d))
  !real(dp),parameter :: unit_B =1./(4.d0*pi*unit_d*(unit_v)**2)
    real(dp),parameter :: unit_B = dsqrt(4*pi*unit_d)*unit_v !If 4pi appears here, should not appear in the def of B_cell

    real(dp),parameter :: unit_e = dsqrt(unit_l*unit_m) !May be a 4pi missing TO CHECK LATER


end module units

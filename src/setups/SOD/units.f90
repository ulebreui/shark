module units
  use precision
  use parameters
  
  implicit none

    real(dp),parameter :: unit_l = 1.
    real(dp),parameter :: unit_d = 1.
    real(dp),parameter :: unit_P = 1.
    real(dp),parameter :: unit_dcol=1.
    real(dp),parameter :: unit_t = 1.
    real(dp),parameter :: unit_m = 1.
    real(dp),parameter :: unit_v = 1.
    real(dp),parameter :: unit_nH= 1.
    real(dp),parameter :: unit_B =1./(4.d0*pi*unit_d*(unit_v)**2)
end module units

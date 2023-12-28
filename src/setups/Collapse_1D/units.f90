module units
  use precision
  use parameters
  
  implicit none

    real(dp),parameter :: unit_l = pc
    real(dp),parameter :: unit_d = mu_gas*mH
    real(dp),parameter :: unit_dcol=unit_d*unit_l
    real(dp),parameter :: unit_t = 1./sqrt(Grav*unit_d)
    real(dp),parameter :: unit_m = unit_l**3*unit_d
    real(dp),parameter :: unit_v = unit_l/unit_t
    real(dp),parameter :: unit_nH= unit_d/(mu_gas*mH)
    real(dp),parameter :: unit_B =1./(4.d0*pi*unit_d*(unit_v)**2)
    real(dp),parameter :: unit_P = unit_d*unit_v**2
end module units

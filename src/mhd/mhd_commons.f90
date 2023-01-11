module mhd_commons
  use precision
  implicit none

  real(dp) :: B_0_lee     = 3d-5     ! Value of the B field at 10^4
  real(dp) :: B_threshold = 0.1d0 ! Value of the B field threshold
  real(dp) :: B_ext       = 3d-5 ! External B field when MHD>0
  
end module mhd_commons

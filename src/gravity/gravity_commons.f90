module gravity_commons
  use precision
  implicit none

  real(dp) :: l_soft = 0.0d0 ! 0.1 AU
  real(dp) :: M_tot  = 0.0d0
  real(dp) :: M_cent = 0.0d0
#if NDUST>0
  real(dp) :: M_tot_dust = 0.0d0
#endif
  real(dp) :: sink_acc_thre = 1d-11
  logical  :: sink =.false.
  integer  :: ir_sink=4
  
end module gravity_commons

module turb_commons

use precision
use phys_const
implicit none

integer :: nb_turb_modes_driven = 10
real(dp):: k_min = 2.0d0*pi/1.0d0
real(dp):: k_max = 2.0d0*pi/(1.0d0/10.0d0)

logical :: new_seed = .false.
integer :: iseed=0
logical :: phase_drift = .true.
integer :: iseed_phase_drift=0
integer :: count=0
integer :: count_bis=0



logical :: driven_turb = .false.
logical :: turb_compressive = .false.
logical :: turb_solenoidal = .false.

logical :: turb_dust = .false.


logical :: ini_random = .true.

real(dp) , dimension(:)  ,allocatable    :: k_turb_driven 
real(dp) :: Mach_nb_ini = 1.0d0
real(dp) :: Mach_yz_target = 1.0d0
real(dp) :: Mach_x_target = 1.0d0


real(dp) :: corrector = 1.0d0
real(dp) :: corrector_sol = 1.0d-1

real(dp) :: V_rms = 1.0d0

real(dp) :: Vy_rms = 1.0d0
real(dp) :: Vz_rms = 1.0d0	
real(dp) :: Vtot_rms = 1.0d0
real(dp) :: Vyz_rms = 1.0d0

real(dp) :: turnover_time = 1.0d0	




real(dp) , dimension(:)  , allocatable    :: random_array_vx
real(dp) , dimension(:)  , allocatable    :: random_array_vy
real(dp) , dimension(:)  , allocatable    :: random_array_vz
real(dp) , dimension(:)  , allocatable    :: random_array_ax
real(dp) , dimension(:)  , allocatable    :: random_array_ay
real(dp) , dimension(:)  , allocatable    :: random_array_az
real(dp) , dimension(:)  , allocatable    :: random_array_phix
real(dp) , dimension(:)  , allocatable    :: random_array_phiy
real(dp) , dimension(:)  , allocatable    :: random_array_phiz
real(dp) , dimension(:)  , allocatable    :: random_array_phix_d
real(dp) , dimension(:)  , allocatable    :: random_array_phiy_d
real(dp) , dimension(:)  , allocatable    :: random_array_phiz_d


end module turb_commons
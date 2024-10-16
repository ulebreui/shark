  module ionisation_commons
  use precision
  implicit none

  logical             ::  charging              = .false.   ! Charging is activated
  logical             ::  charging_all_the_time = .false.   ! Charging is activated
  logical             ::  analytical_charging = .false.   ! Prescription from (Fujii et. al 2011) and see Lebreuilly 2020.
  logical             :: res_Marchand = .true.   ! To compute charges within Marchand algorithm
  logical             :: dust_inertia = .false.   ! To compute res accounting for dust inertia
  logical             :: electrons = .true.   ! To include electrons in the resistivitiy computation
  logical             :: ions = .true.   
  logical             :: dusty_nonideal_MHD = .false.
  logical             :: dusty_nonideal_MHD_no_electron = .false.   
   





  real(dp), dimension(:), allocatable      :: eta_a
  real(dp), dimension(:), allocatable      :: eta_o
  real(dp), dimension(:), allocatable      :: eta_h

  real(dp), dimension(:), allocatable      :: eta_eff_yy
  real(dp), dimension(:), allocatable      :: eta_eff_yz
  real(dp), dimension(:), allocatable      :: eta_eff_zy
  real(dp), dimension(:), allocatable      :: eta_eff_zz

  real(dp), dimension(:), allocatable      :: eta_eff_ohm
  real(dp), dimension(:), allocatable      :: eta_eff_Hall_y
  real(dp), dimension(:), allocatable      :: eta_eff_Hall_z



  real(dp), dimension(:), allocatable      :: sigma_o
  real(dp), dimension(:), allocatable      :: sigma_p
  real(dp), dimension(:), allocatable      :: sigma_h

  real(dp), dimension(:), allocatable      :: ni
  real(dp), dimension(:), allocatable      :: ne
  real(dp), dimension(:), allocatable      :: Hall_e
  real(dp), dimension(:), allocatable      :: Hall_i

  real(dp), dimension(:), allocatable      :: E_x
  real(dp), dimension(:), allocatable      :: E_y
  real(dp), dimension(:), allocatable      :: E_z

  real(dp), dimension(:), allocatable      :: v_i_x
  real(dp), dimension(:), allocatable      :: v_i_y
  real(dp), dimension(:), allocatable      :: v_i_z

  real(dp), dimension(:), allocatable      :: v_e_x
  real(dp), dimension(:), allocatable      :: v_e_y
  real(dp), dimension(:), allocatable      :: v_e_z

  real(dp), dimension(:), allocatable      :: FLor_x
  real(dp), dimension(:), allocatable      :: FLor_y
  real(dp), dimension(:), allocatable      :: FLor_z



  real(dp), dimension(:), allocatable      :: psi_old

#if NDUST>0
  real(dp), dimension(:,:), allocatable    :: zd, gamma_d
  real(dp), dimension(:,:), allocatable      :: FLor_x_d
  real(dp), dimension(:,:), allocatable      :: FLor_y_d
  real(dp), dimension(:,:), allocatable      :: FLor_z_d
#endif

  real(dp) :: B_0_lee       = 3d-5   ! Value of the B field at 10^4
  real(dp) :: B_threshold   = 0.1d0  ! Value of the B field threshold
  real(dp)::  mu_ions       = 25.0d0 ! Mean ion mass (in units of mh)
  real(dp):: stickeff_el    = 0.5d0  ! Sticking coefficient of electrions
  real(dp):: epsilon_ionis  = 1d-6   ! Tolerance of the ionisation scheme
  integer :: nitermax_ionis = 1000   ! Maximum number of iterations
  real(dp):: x              = 5d-17  ! CR Ionisation rate


  real(dp) :: f_dust_charge = 0.0d0 !Derivative of function to solve in analytical charge scheme
  real(dp) :: fprime_dust_charge = 0.0d0 !function to solve in analytical charge scheme

end module ionisation_commons
  module ionisation_commons
  use precision
  implicit none

  logical             ::  charging              = .false.   ! Charging is activated
  logical             ::  charging_all_the_time = .false.   ! Charging is activated

  real(dp), dimension(:), allocatable      :: eta_a
  real(dp), dimension(:), allocatable      :: eta_o
  real(dp), dimension(:), allocatable      :: eta_h

  real(dp), dimension(:), allocatable      :: sigma_o
  real(dp), dimension(:), allocatable      :: sigma_p
  real(dp), dimension(:), allocatable      :: sigma_h

  real(dp), dimension(:), allocatable      :: ni
  real(dp), dimension(:), allocatable      :: ne

  real(dp), dimension(:), allocatable      :: psi_old

#if NDUST>0
  real(dp), dimension(:,:), allocatable    :: zd, gamma_d
#endif

  real(dp) :: B_0_lee       = 3d-5   ! Value of the B field at 10^4
  real(dp) :: B_threshold   = 0.1d0  ! Value of the B field threshold
  real(dp)::  mu_ions       = 25.0d0 ! Mean ion mass (in units of mh)
  real(dp):: stickeff_el    = 0.5d0  ! Sticking coefficient of electrions
  real(dp):: epsilon_ionis  = 1d-6   ! Tolerance of the ionisation scheme
  integer :: nitermax_ionis = 1000   ! Maximum number of iterations
  real(dp):: x              = 5d-17  ! CR Ionisation rate

end module ionisation_commons
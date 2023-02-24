module hydro_commons
  use precision
  !real(dp), dimension(:,:), allocatable :: unew
  real(dp), dimension(:,:), allocatable :: u_prim
  real(dp), dimension(:),   allocatable :: cs
  real(dp), dimension(:,:), allocatable :: q
  real(dp), dimension(:), allocatable   :: unit_var ! Unit of the variable in the output. 1. by defaut Change in setup.f90

  real(dp), dimension(:,:,:), allocatable :: qm
  real(dp), dimension(:,:,:), allocatable :: qp
  real(dp), dimension(:,:),   allocatable :: force

  real(dp) :: gamma       = 1.66667d0   ! Adiabatic index
  real(dp) :: iso_cs      = -1.0d0 ! Isothermal eos if positiv
  integer :: irho
  integer :: iv
  integer :: ivy
  integer :: ivz
  integer :: iP

end module hydro_commons

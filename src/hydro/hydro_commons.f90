module hydro_commons
  use precision
  real(dp), dimension(:,:,:), allocatable :: q
  real(dp), dimension(:,:,:),   allocatable :: u_prim

  real(dp), dimension(:,:),   allocatable :: phi_grav
  real(dp), dimension(:,:),   allocatable :: grad_phi_sg_x
  real(dp), dimension(:,:),   allocatable :: grad_phi_sg_y


  real(dp), dimension(:,:),     allocatable :: cs

  real(dp), dimension(:,:,:), allocatable :: qm_x
  real(dp), dimension(:,:,:), allocatable :: qp_x
  real(dp), dimension(:,:,:), allocatable :: qm_y
  real(dp), dimension(:,:,:), allocatable :: qp_y


  real(dp), dimension(:,:),   allocatable :: force_x
  real(dp), dimension(:,:),   allocatable :: force_y
  real(dp), dimension(:,:),   allocatable :: force_z

  real(dp), dimension(:,:,:), allocatable  :: flux_x
  real(dp), dimension(:,:,:), allocatable  :: flux_y

  real(dp) :: gamma           = 1.66667d0   ! Adiabatic index
  integer :: iso_cs           = -1          ! Isothermal eos if positiv
  integer :: non_standard_eos = -1          ! Switch to > 0 for barotropic eos or user defined eos
  
  real(dp):: smallr = 1d-27
  real(dp):: smallp = 1d-27
  real(dp):: smallc = 1d-27

  integer :: slope_type=1 ! minmod, 2 VL

  integer :: irho
  integer :: ivx
  integer :: ivy
  integer :: ivz
  integer :: iP

#if NDUST>0
  integer, dimension(1:NDUST,1:2) :: index_vdn
  integer, dimension(1:NDUST,1:2) :: index_vdt
#endif
  integer, dimension(1:2) :: index_vn
  integer, dimension(1:2) :: index_vt

  logical :: no_flux_x_rho_in = .false.
  real(dp) :: rho_sink = 1e11

end module hydro_commons

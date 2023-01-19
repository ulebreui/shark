module hydro_commons
  use precision
  real(dp), dimension(:,:), allocatable :: unew
  real(dp), dimension(:,:), allocatable :: uold
  real(dp), dimension(:),   allocatable :: cs
  real(dp), dimension(:,:), allocatable :: q
  real(dp), dimension(:,:,:), allocatable :: ql
  real(dp), dimension(:,:,:), allocatable :: qr
  real(dp), dimension(:),   allocatable :: force

  real(dp) :: gamma       = 1.66667d0   ! Adiabatic index

  integer :: irho
  integer :: iv
  integer :: ivy
  integer :: iP

end module hydro_commons

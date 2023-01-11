module hydro_commons
  use precision
  real(dp), dimension(:,:), allocatable :: unew
  real(dp), dimension(:,:), allocatable :: uold
  real(dp), dimension(:),   allocatable :: cs
  real(dp), dimension(:),   allocatable :: csl
  real(dp), dimension(:),   allocatable :: csr

  real(dp), dimension(:,:), allocatable :: q
  real(dp), dimension(:,:), allocatable :: ql
  real(dp), dimension(:,:), allocatable :: qr

  real(dp), dimension(:),   allocatable :: force
  
#if GRAVITY==1
  real(dp), dimension(:), allocatable :: Mc
#endif
  real(dp), dimension(:), allocatable :: B_cell

  integer :: irho
  integer :: iv
  integer :: iP
end module hydro_commons

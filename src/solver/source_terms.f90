subroutine Source_terms
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust
  real(dp), dimension(:,:)  , allocatable :: S_U
  if(static) return
  allocate(S_U(1:ncells,1:nvar))
  S_U=0.0d0
  !Update state vector
  unew=unew+S_U
  deallocate(S_U)
  !call kick
  
end subroutine Source_terms

! Force kick (if applicable)
subroutine kick
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust
  
end subroutine kick

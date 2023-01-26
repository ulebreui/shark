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
#if NDUST>0
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
      do idust=1,ndust
        S_U(i,iP)=S_U(i,iP)+dt*(q(i,iv)*(q(i,iv)-q(i,ivd(idust))))
        if(ndim==2) S_U(i,iP)=S_U(i,iP)+dt*(q(i,ivy)*(q(i,ivy)-q(i,ivdy(idust))))
      end do
    endif
  end do
  !$OMP END DO
  !$OMP END PARALLEL
#endif

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

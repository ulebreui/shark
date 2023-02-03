subroutine Source_terms
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust
  real(dp) :: ts
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
    if(dust_back_reaction) then
      do idust=1,ndust
        ts  = sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(i,idust)/(q(i,irho)*cs(i))
        S_U(i,iP) = S_U(i,iP)+dt*(q(i,iv)*(q(i,iv)-q(i,ivd(idust))))/ts
        if(ndim==2) S_U(i,iP) = S_U(i,iP)+dt*(q(i,ivy)*(q(i,ivy)-q(i,ivdy(idust))))/ts
      end do
      endif
    endif
  end do
  !$OMP END DO
  !$OMP END PARALLEL
#endif

  !Update state vector
  unew=unew+S_U
  deallocate(S_U)
  call kick
  
end subroutine Source_terms

! Force kick (if applicable)
subroutine kick
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust
  do i=1,ncells
    if(active_cell(i)==1) then
    unew(i,iv)=unew(i,iv) + q(i,irho)*dt*force(i,1)
    if(ndim==2) unew(i,ivy)=unew(i,ivy) + q(i,irho)*dt*force(i,2)
    endif
  end do
end subroutine kick

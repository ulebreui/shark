subroutine Source_terms
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust
  real(dp), dimension(1:ncells,1:nvar) :: S_U
  if(static) return
  S_U=0.0d0
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do i = first_active,last_active
#if SPHERE==1
     S_U(i,iv)=-dt*(- 2.*q(i,irho)*cs(i)**2./r_c(i))
#endif     
#if GRAVITY==1     
     S_U(i,iv)= S_U(i,iv)-dt*(q(i,irho)*Mc(i)/(r_c(i)**2.+(l_soft/unit_l)**2.))
#if NDUST>0     
     do idust=1,ndust
         S_U(i,ivd(idust))= S_U(i,ivd(idust))-dt*(q(i,irhod(idust))*Mc(i)/(r_c(i)**2.+(l_soft/unit_l)**2.)) 
      end do
#endif
#endif
   end do
  !$OMP END DO
  !$OMP END PARALLEL
  !Update state vector
  unew=unew+S_U
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
    do i = first_active,last_active
     unew(i,iv)=unew(i,iv)+unew(i,irho)*dt*force(i)
#if NDUST>0     
     do idust=1,ndust
     unew(i,ivd(idust))=unew(i,ivd(idust))+unew(i,irhod(idust))*dt*force(i)
     end do
#endif
  end do
end subroutine kick

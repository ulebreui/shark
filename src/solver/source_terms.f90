subroutine Source_terms
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust
  real(dp) :: ts,ekin
  real(dp), dimension(:,:)  , allocatable :: S_U
  if(static) return
  allocate(S_U(1:ncells,1:nvar))
  S_U=0.0d0

#if NDUST>0
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,ts)
  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
      if(dust_back_reaction) then
        do idust = 1,ndust
          ts  = sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(i,idust)/(q(i,irho)*cs(i))
          S_U(i,iP) = S_U(i,iP)+half*dt*(q(i,iv)*(q(i,iv)-q(i,ivd(idust))))/ts
#if NY>1          
          S_U(i,iP) = S_U(i,iP)+half*dt*(q(i,ivy)*(q(i,ivy)-q(i,ivdy(idust))))/ts
#if IVZ==1
          S_U(i,iP) = S_U(i,iP)+half*dt*(q(i,ivz)*(q(i,ivz)-q(i,ivdz(idust))))/ts
#endif    
#endif      
        end do
      endif
    endif
  end do
  !$OMP END DO
  !$OMP END PARALLEL
#endif

  !Update state vector
  u_prim=u_prim+S_U
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
  real(dp) :: energy_old,energy_new
  if(.not. force_kick ) return
#if TURB==1
  call turbulent_force
#endif
  call update_force_setup ! Setup specific update of the force

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,energy_old,energy_new)
  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
      energy_old = 0.5d0*q(i,irho)*q(i,iv)**2
      u_prim(i,iv)=u_prim(i,iv) + q(i,irho)*dt*force(i,1)*half
      energy_new = 0.5d0*u_prim(i,iv)**2/q(i,irho)

#if NY>1  
      energy_old =  energy_old +  0.5d0*q(i,irho)*q(i,ivy)**2
      u_prim(i,ivy)=u_prim(i,ivy) + q(i,irho)*dt*force(i,2)*half
      energy_new=  energy_new +  0.5d0*u_prim(i,ivy)**2/q(i,irho)

#if IVZ==1
      energy_old =  energy_old +  0.5d0*q(i,irho)*q(i,ivz)**2
      u_prim(i,ivz)=u_prim(i,ivz) + q(i,irho)*dt*force(i,3)*half
      energy_new=  energy_new +  0.5d0*u_prim(i,ivz)**2/q(i,irho)
#endif 
#endif    
      u_prim(i,iP) = u_prim(i,iP) +  energy_new - energy_old
#if NDUST>0
  do idust=1,ndust
    u_prim(i,ivd(idust))=u_prim(i,ivd(idust)) + q(i,irhod(idust))*dt*force_dust(i,1,idust)*half
    if(ndim==2) u_prim(i,ivdy(idust))=u_prim(i,ivdy(idust)) + q(i,irhod(idust))*dt*force_dust(i,2,idust)*half
#if IVZ==1
      u_prim(i,ivdz(idust))=u_prim(i,ivdz(idust)) + q(i,irhod(idust))*dt*force_dust(i,3,idust)*half
#endif      
  end do
#endif
  endif
end do
  !$OMP END DO
  !$OMP END PARALLEL 

end subroutine kick

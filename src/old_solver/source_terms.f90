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
  call compute_tstop
#endif
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,ts)
  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
#if NDUST>0
      if(dust_back_reaction) then
        do idust = 1,ndust
          S_U(i,iP) = S_U(i,iP)+dt*(q(i,ivx)*(q(i,ivx)-q(i,ivdx(idust))))/tstop(i,idust)
#if NY>1          
          S_U(i,iP) = S_U(i,iP)+dt*(q(i,ivy)*(q(i,ivy)-q(i,ivdy(idust))))/tstop(i,idust)
          S_U(i,iP) = S_U(i,iP)+dt*(q(i,ivz)*(q(i,ivz)-q(i,ivdz(idust))))/tstop(i,idust)
#endif      
        end do
      endif
#endif
    endif
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !Update state vector
  u_prim=u_prim+S_U
  deallocate(S_U)

end subroutine Source_terms

! Force kick (if applicable)
subroutine kick(coeffdt)
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust
  real(dp) :: energy_old,energy_new,coeffdt

  call update_force_setup
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,energy_old,energy_new)
  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
      energy_old    = half*u_prim(i,ivx)**2/u_prim(i,irho)
      u_prim(i,ivx)  = u_prim(i,ivx) + u_prim(i,irho)*coeffdt*dt*force(i,1)
      energy_new    = half*u_prim(i,ivx)**2/u_prim(i,irho)
#if NY>1  
      energy_old    =  energy_old +  half*u_prim(i,ivy)**2/u_prim(i,irho)
      u_prim(i,ivy) =  u_prim(i,ivy) + u_prim(i,irho)*coeffdt*dt*force(i,2)
      energy_new    =  energy_new +  half*u_prim(i,ivy)**2/u_prim(i,irho)

      energy_old    =  energy_old    + half*u_prim(i,ivz)**2/u_prim(i,irho)
      u_prim(i,ivz) =  u_prim(i,ivz) + u_prim(i,irho)*coeffdt*dt*force(i,3)
      energy_new    =  energy_new    + half*u_prim(i,ivz)**2/u_prim(i,irho)
#endif    
      u_prim(i,iP) = u_prim(i,iP) +  (energy_new - energy_old)
#if NDUST>0
  do idust=1,ndust
    u_prim(i,ivdx(idust))              = u_prim(i,ivdx(idust))  + u_prim(i,irhod(idust))*coeffdt*dt*force_dust(i,1,idust)
#if NY>1
    u_prim(i,ivdy(idust))             = u_prim(i,ivdy(idust)) + u_prim(i,irhod(idust))*coeffdt*dt*force_dust(i,2,idust)
    u_prim(i,ivdz(idust))            = u_prim(i,ivdz(idust)) + u_prim(i,irhod(idust))*coeffdt*dt*force_dust(i,3,idust)
#endif      
  end do
#endif
  endif
end do
  !$OMP END DO
  !$OMP END PARALLEL 

end subroutine kick



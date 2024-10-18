! Force kick (if applicable)
subroutine kick(coeffdt)
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust
  real(dp) :: energy_old,energy_new,coeffdt

  if(static) return
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,energy_old,energy_new)
  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
      if(iso_cs.ne.1 .or. non_standard_eos .ne. 1) then
        energy_old     = half*u_prim(i,ivx)**2/u_prim(i,irho) +  half*u_prim(i,ivy)**2/u_prim(i,irho)  + half*u_prim(i,ivz)**2/u_prim(i,irho)
      endif
        u_prim(i,ivx)  = u_prim(i,ivx)  + u_prim(i,irho)*coeffdt*dt*force(i,1)
        u_prim(i,ivy)  = u_prim(i,ivy)  + u_prim(i,irho)*coeffdt*dt*force(i,2)
        u_prim(i,ivz)  = u_prim(i,ivz)  + u_prim(i,irho)*coeffdt*dt*force(i,3)

      if(iso_cs.ne.1 .or. non_standard_eos .ne. 1) then
        energy_new     = half*u_prim(i,ivx)**2/u_prim(i,irho) +  half*u_prim(i,ivy)**2/u_prim(i,irho) + half*u_prim(i,ivz)**2/u_prim(i,irho)
        u_prim(i,iP)   = u_prim(i,iP) +  (energy_new - energy_old)
      endif
#if NDUST>0
  do idust=1,ndust
    u_prim(i,ivdx(idust))            = u_prim(i,ivdx(idust)) + u_prim(i,irhod(idust))*coeffdt*dt*force_dust(i,1,idust)
    u_prim(i,ivdy(idust))            = u_prim(i,ivdy(idust)) + u_prim(i,irhod(idust))*coeffdt*dt*force_dust(i,2,idust)
    u_prim(i,ivdz(idust))            = u_prim(i,ivdz(idust)) + u_prim(i,irhod(idust))*coeffdt*dt*force_dust(i,3,idust)    
  end do
#endif
  endif
end do
  !$OMP END DO
  !$OMP END PARALLEL 

end subroutine kick


! Force kick (if applicable)
subroutine kick_q(coeffdt)
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust
  real(dp) :: coeffdt

  if(static) return
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
        q(i,ivx)  = q(i,ivx)  + coeffdt*dt*force(i,1)
        q(i,ivy)  = q(i,ivy)  + coeffdt*dt*force(i,2)
        q(i,ivz)  = q(i,ivz)  + coeffdt*dt*force(i,3)

#if NDUST>0
  do idust=1,ndust
        q(i,ivdx(idust))            = q(i,ivdx(idust)) + coeffdt*dt*force_dust(i,1,idust)
        q(i,ivdy(idust))            = q(i,ivdy(idust)) + coeffdt*dt*force_dust(i,2,idust)
        q(i,ivdz(idust))            = q(i,ivdz(idust)) + coeffdt*dt*force_dust(i,3,idust)    
  end do
#endif
  endif
end do
  !$OMP END DO
  !$OMP END PARALLEL 

end subroutine kick_q

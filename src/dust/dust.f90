! This routine performs the dust step
subroutine dust(verbose,outputing)
  use parameters
  use commons
  use units
  implicit none
  logical :: verbose,outputing

  !We account for the gas-dust drag
  call dust_drag
  
end subroutine dust


! Dust drag is computed (implicitely)
subroutine dust_drag
  
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust
  real(dp),dimension(1:ndust):: alphak
  real(dp):: pn,rhon
  if(static) return
  !Re-calc distribution
  call distribution_dust(.false.)
  call compute_tstop
  ! Here we apply the Krapp et al. implict scheme to compute the dust drag source terms

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(pn,rhon,alphak,i,idust)
  !$OMP DO
  do i=1,ncells
   if(active_cell(i)==1) then
     ! Direction - x
     pn   = u_prim(i,iv)
     rhon = u_prim(i,irho)
     do idust=1,ndust
        alphak(idust) = 1.0d0/tstop(i,idust)
        pn   = pn+(alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*q(i,irhod(idust))*q(i,ivd(idust))
     end do
     do idust = 1,ndust
        u_prim(i,ivd(idust)) = q(i,irhod(idust))*(q(i,ivd(idust))/(1.0d0+alphak(idust)*dt)+ (alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*pn/rhon)
     end do 
     if(dust_back_reaction)u_prim(i,iv) = pn/rhon*u_prim(i,irho)
     ! Direction - y
     pn   = u_prim(i,ivy)
     do idust=1,ndust
        pn   = pn+(alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*q(i,irhod(idust))*q(i,ivdy(idust))
     end do
     do idust = 1,ndust
        u_prim(i,ivdy(idust)) = q(i,irhod(idust))*(q(i,ivdy(idust))/(1.0d0+alphak(idust)*dt)+ (alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*pn/rhon)
     end do 
     if(dust_back_reaction)u_prim(i,ivy) = pn/rhon*u_prim(i,irho)   
#if IVZ==1
     !Direction - z
     pn   = u_prim(i,ivz)
     do idust=1,ndust
        ! alpha = 1/ts
        pn   = pn+(alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*q(i,irhod(idust))*q(i,ivdz(idust))
     end do
     do idust = 1,ndust
        u_prim(i,ivdz(idust)) = q(i,irhod(idust))*(q(i,ivdz(idust))/(1.0d0+alphak(idust)*dt)+ (alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*pn/rhon)
     end do 
     if(dust_back_reaction) u_prim(i,ivz) = pn/rhon*u_prim(i,irho) 
#endif      
     end if
  end do
   !$OMP END DO
   !$OMP BARRIER
   !$OMP DO
  ! Regularisation to avoid negative dust densities
  do i=1,ncells
  if(active_cell(i)==1) then
      do idust=1,ndust
         u_prim(i,irhod(idust)) = max(u_prim(i,irho)*dust_ratio_min, u_prim(i,irhod(idust)))
      enddo
  endif
  end do
  !$OMP END PARALLEL 
end subroutine dust_drag


! Dust drag is computed (implicitely)
subroutine dust_drag(coeffdt)
  
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust
  real(dp), dimension(1:ndust):: alphak
  real(dp):: pnx,pny,pnz,rhon,coeffdt
  if(static) return

  ! Here we apply the Krapp et al. implict scheme to compute the dust drag source terms

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(pnx,pny,pnz,rhon,alphak,i,idust)
  !$OMP DO
  do i=1,ncells
   if(active_cell(i)) then

     rhon  = u_prim(i,irho)
     pnx   = u_prim(i,ivx)  
     pny   = u_prim(i,ivy)  
     pnz   = u_prim(i,ivz)  

     do idust=1,ndust

        alphak(idust) = coeffdt * dt / tstop(i,idust) ! Half for half dt

        pnx   = pnx  + alphak(idust)/(1.0d0 + alphak(idust)) * u_prim(i,ivdx(idust))
        pny   = pny  + alphak(idust)/(1.0d0 + alphak(idust)) * u_prim(i,ivdy(idust))
        pnz   = pnz  + alphak(idust)/(1.0d0 + alphak(idust)) * u_prim(i,ivdz(idust))
        rhon  = rhon + alphak(idust)/(1.0d0 + alphak(idust)) * u_prim(i,irhod(idust))

     end do

     do idust = 1,ndust

        u_prim(i,ivdx(idust)) = u_prim(i,ivdx(idust))/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pnx/rhon*u_prim(i,irhod(idust))
        u_prim(i,ivdy(idust)) = u_prim(i,ivdy(idust))/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pny/rhon*u_prim(i,irhod(idust))
        u_prim(i,ivdz(idust)) = u_prim(i,ivdz(idust))/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pnz/rhon*u_prim(i,irhod(idust))

     end do 

     if(dust_back_reaction) then

         u_prim(i,ivx)  = pnx/rhon * u_prim(i,irho) 
         u_prim(i,ivy)  = pny/rhon * u_prim(i,irho)   
         u_prim(i,ivz)  = pnz/rhon * u_prim(i,irho)      

      end if 

     end if
  end do
   !$OMP END DO
   !$OMP BARRIER
   !$OMP DO
   
  ! Regularisation to avoid negative dust densities
  do i=1,ncells
   if(active_cell(i)) then

      do idust=1,ndust

         u_prim(i,irhod(idust)) = max(u_prim(i,irho) * dust_ratio_min, u_prim(i,irhod(idust)))

      enddo

   endif
  end do
  !$OMP END PARALLEL 
end subroutine dust_drag





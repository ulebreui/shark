! Dust drag is computed (implicitely)
subroutine dust_drag(coeffdt)
  
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust,ix,iy,icell
  real(dp), dimension(1:ndust):: alphak
  real(dp):: pnx,pny,pnz,rhon,coeffdt
  if(static) return

  ! Here we apply the Krapp et al. implict scheme to compute the dust drag source terms

  do iy=first_active_y,last_active_y
   do ix=first_active,last_active
   
     i = icell(ix,iy)
     rhon  = u_prim(ix,iy,irho)
     pnx   = u_prim(ix,iy,ivx)  
     pny   = u_prim(ix,iy,ivy)  
     pnz   = u_prim(ix,iy,ivz)  

     do idust=1,ndust

        alphak(idust) = coeffdt * dt / tstop(ix,iy,idust) ! Half for half dt

        pnx   = pnx  + alphak(idust)/(1.0d0 + alphak(idust)) * u_prim(ix,iy,ivdx(idust))
        pny   = pny  + alphak(idust)/(1.0d0 + alphak(idust)) * u_prim(ix,iy,ivdy(idust))
        pnz   = pnz  + alphak(idust)/(1.0d0 + alphak(idust)) * u_prim(ix,iy,ivdz(idust))
        rhon  = rhon + alphak(idust)/(1.0d0 + alphak(idust)) * u_prim(ix,iy,irhod(idust))

     end do

     do idust = 1,ndust

        u_prim(ix,iy,ivdx(idust)) = u_prim(ix,iy,ivdx(idust))/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pnx/rhon*u_prim(ix,iy,irhod(idust))
        u_prim(ix,iy,ivdy(idust)) = u_prim(ix,iy,ivdy(idust))/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pny/rhon*u_prim(ix,iy,irhod(idust))
        u_prim(ix,iy,ivdz(idust)) = u_prim(ix,iy,ivdz(idust))/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pnz/rhon*u_prim(ix,iy,irhod(idust))

     end do 

     if(dust_back_reaction) then

         u_prim(ix,iy,ivx)  = pnx/rhon * u_prim(ix,iy,irho) 
         u_prim(ix,iy,ivy)  = pny/rhon * u_prim(ix,iy,irho)   
         u_prim(ix,iy,ivz)  = pnz/rhon * u_prim(ix,iy,irho)      

      end if 

     end do
  end do

  ! Regularisation to avoid negative dust densities
  do idust=1,ndust
      do iy=first_active_y,last_active_y
         do ix=first_active,last_active
            u_prim(ix,iy,irhod(idust)) = max(u_prim(ix,iy,irho) * dust_ratio_min, u_prim(ix,iy,irhod(idust)))
      enddo
   end do
  end do

  
end subroutine dust_drag





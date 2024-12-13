! Dust drag is computed (implicitely)
subroutine dust_drag(coeffdt)
  
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust
  real(dp), dimension(1:ndust):: alphak
  real(dp):: pnx,pny,pnz,rhon,coeffdt,B_norm
  if(static) return

  ! Here we apply the Krapp et al. implict scheme to compute the dust drag source terms

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(pnx,pny,pnz,rhon,alphak,B_norm,i,idust)
  !$OMP DO
  do i=1,ncells
   if(active_cell(i)==1) then

     rhon  = u_prim(i,irho)
     pnx   = u_prim(i,ivx)  
     pny   = u_prim(i,ivy)  
     pnz   = u_prim(i,ivz)  

#if MHD==1

     B_norm=dsqrt(u_prim(i,iBx)**2+u_prim(i,iBy)**2+u_prim(i,iBz)**2)

#endif


     do idust=1,ndust

        alphak(idust) = coeffdt * dt / tstop(i,idust) ! Half for half dt

#if MHD==1

        if (dusty_nonideal_MHD_no_electron .and. idust==i_coupled_species) then !Works for a single grain only

               !!!Effective alpha due to extra dust/gas coupling caused by ions.

                alphak(idust) = coeffdt * dt * (1.0 / tstop(i,idust) - e_el_stat*zd(i,idust)*B_norm/(clight*Hall_i(i)*mdust(i,idust)))  ! Half for half dt

         endif
#endif

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
   if(active_cell(i)==1) then

      do idust=1,ndust

         u_prim(i,irhod(idust)) = max(u_prim(i,irho) * dust_ratio_min, u_prim(i,irhod(idust)))

      enddo

   endif
  end do
  !$OMP END PARALLEL 
end subroutine dust_drag





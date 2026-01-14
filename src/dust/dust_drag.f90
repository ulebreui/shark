! Dust drag is computed (implicitely)
subroutine dust_drag

   use parameters
   use commons
   use units
   use OMP_LIB

   implicit none
   integer :: idust, ix,iy,ivar
   real(dp), dimension(1:ndust):: alphak
   real(dp), dimension(1:nvar):: uloc

   real(dp):: pnx, pny, pnz, rhon, coeffdt
   if (static) return

   if (dust_back_reaction) then

   ! Here we apply the Krapp et al. implict scheme to compute the dust drag source terms
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy,alphak, pnx, pny, pnz, rhon,uloc,ivar)
   do iy = first_active_y, last_active_y
      !$omp simd
      do ix = first_active, last_active


         uloc(:) =u_prim(:,ix,iy)
         
         rhon = uloc(irho)
         pnx  = uloc(ivx)
         pny  = uloc(ivy)
         pnz  = uloc(ivz)

         do idust = 1, ndust

            alphak(idust) = dt/tstop(idust,ix,iy) ! Half for half dt

            pnx  = pnx  + alphak(idust)/(1.0d0 + alphak(idust))*uloc(ivdx(idust))
            pny  = pny  + alphak(idust)/(1.0d0 + alphak(idust))*uloc(ivdy(idust))
            pnz  = pnz  + alphak(idust)/(1.0d0 + alphak(idust))*uloc(ivdz(idust))
            rhon = rhon + alphak(idust)/(1.0d0 + alphak(idust))*uloc(irhod(idust))

         end do

         do idust = 1, ndust

            u_prim(ivdx(idust),ix,iy) = uloc(ivdx(idust))/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pnx/rhon*uloc(irhod(idust))
            u_prim(ivdy(idust),ix,iy) = uloc(ivdy(idust))/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pny/rhon*uloc(irhod(idust))
            u_prim(ivdz(idust),ix,iy) = uloc(ivdz(idust))/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pnz/rhon*uloc(irhod(idust))
            
            ! Regularisation of the dust to avoid to small dtg ratio
            u_prim(irhod(idust),ix,iy) = max(uloc(irho)*dust_ratio_min, uloc(irhod(idust)))
         end do

         

         u_prim(ivx,ix,iy) = pnx/rhon*uloc(irho)
         u_prim(ivy,ix,iy) = pny/rhon*uloc(irho)
         u_prim(ivz,ix,iy) = pnz/rhon*uloc(irho)

               
      end do
   end do
   
   else

   ! Here we apply the Krapp et al. implict scheme to compute the dust drag source terms
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy,alphak, pnx, pny, pnz, rhon)
   do iy = first_active_y, last_active_y
      !$omp simd
      do ix = first_active, last_active

         rhon = u_prim(irho,ix,iy)
         pnx  = u_prim(ivx,ix,iy)
         pny  = u_prim(ivy,ix,iy)
         pnz  = u_prim(ivz,ix,iy)

         do idust = 1, ndust

            alphak(idust) = dt/tstop(idust,ix,iy) ! Half for half dt

            pnx  = pnx  + alphak(idust)/(1.0d0 + alphak(idust))*u_prim(ivdx(idust),ix,iy)
            pny  = pny  + alphak(idust)/(1.0d0 + alphak(idust))*u_prim(ivdy(idust),ix,iy)
            pnz  = pnz  + alphak(idust)/(1.0d0 + alphak(idust))*u_prim(ivdz(idust),ix,iy)
            rhon = rhon + alphak(idust)/(1.0d0 + alphak(idust))*u_prim(irhod(idust),ix,iy)

         end do

         do idust = 1, ndust

            u_prim(ivdx(idust),ix,iy) = u_prim(ivdx(idust),ix,iy)/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pnx/rhon*u_prim(irhod(idust),ix,iy)
            u_prim(ivdy(idust),ix,iy) = u_prim(ivdy(idust),ix,iy)/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pny/rhon*u_prim(irhod(idust),ix,iy)
            u_prim(ivdz(idust),ix,iy) = u_prim(ivdz(idust),ix,iy)/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pnz/rhon*u_prim(irhod(idust),ix,iy)
            
            ! Regularisation of the dust to avoid to small dtg ratio
            u_prim(irhod(idust),ix,iy) = max(u_prim(irho,ix,iy)*dust_ratio_min, u_prim(irhod(idust),ix,iy))
         end do

         

            ! u_prim(ivx,ix,iy) = pnx/rhon*u_prim(irho,ix,iy)
            ! u_prim(ivy,ix,iy) = pny/rhon*u_prim(irho,ix,iy)
            ! u_prim(ivz,ix,iy) = pnz/rhon*u_prim(irho,ix,iy)

               
         end do
      end do

   endif

end subroutine dust_drag


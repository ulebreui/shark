! Dust drag is computed (implicitely)
subroutine dust_drag

   use parameters
   use commons
   use units
   use OMP_LIB

   implicit none
   integer :: idust, ix,iy
   real(dp), dimension(1:ndust):: alphak
   real(dp):: pnx, pny, pnz, rhon, coeffdt
   if (static) return

   ! Here we apply the Krapp et al. implict scheme to compute the dust drag source terms
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy,alphak, pnx, pny, pnz, rhon)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active

         rhon = u_prim(irho,ix,iy)
         pnx  = u_prim(ivx,ix,iy)
         pny  = u_prim(ivy,ix,iy)
         pnz  = u_prim(ivz,ix,iy)

         do idust = 1, ndust

            alphak(idust) = dt/tstop(idust,ix,iy) ! Half for half dt

            pnx  = pnx + alphak(idust)/(1.0d0 + alphak(idust))*u_prim(ivdx(idust),ix,iy)
            pny  = pny + alphak(idust)/(1.0d0 + alphak(idust))*u_prim(ivdy(idust),ix,iy)
            pnz  = pnz + alphak(idust)/(1.0d0 + alphak(idust))*u_prim(ivdz(idust),ix,iy)
            rhon = rhon + alphak(idust)/(1.0d0 + alphak(idust))*u_prim(irhod(idust),ix,iy)

         end do

         do idust = 1, ndust

            u_prim(ivdx(idust),ix,iy) = u_prim(ivdx(idust),ix,iy)/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pnx/rhon*u_prim(irhod(idust),ix,iy)
            u_prim(ivdy(idust),ix,iy) = u_prim(ivdy(idust),ix,iy)/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pny/rhon*u_prim(irhod(idust),ix,iy)
            u_prim(ivdz(idust),ix,iy) = u_prim(ivdz(idust),ix,iy)/(1.0d0 + alphak(idust)) + (alphak(idust)/(1.0d0 + alphak(idust)))*pnz/rhon*u_prim(irhod(idust),ix,iy)

         end do

         if (dust_back_reaction) then

            u_prim(ivx,ix,iy) = pnx/rhon*u_prim(irho,ix,iy)
            u_prim(ivy,ix,iy) = pny/rhon*u_prim(irho,ix,iy)
            u_prim(ivz,ix,iy) = pnz/rhon*u_prim(irho,ix,iy)

         end if

      end do
   end do
   
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy)
   ! Regularisation to avoid negative dust densities
      do iy = first_active_y, last_active_y
         do ix = first_active, last_active
            do idust = 1, ndust
               u_prim(irhod(idust),ix,iy) = max(u_prim(irho,ix,iy)*dust_ratio_min, u_prim(irhod(idust),ix,iy))
         end do
      end do
   end do

end subroutine dust_drag


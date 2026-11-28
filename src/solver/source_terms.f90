subroutine Source_terms
   use parameters
   use commons
   use OMP_LIB
   use units
   implicit none
   integer :: idust, ix,iy,ivar
   real(dp) :: cs_eos, barotrop,S_r,S_l,P_l,P_r,rho_l,rho_r,vx_l,vx_r,vy_l,vy_r,r_l,r_r
   real(dp) :: cot_th, rr, dAx, dAy
   if (static) then
      return
   end if

#if GEOM==2
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy,S_r,S_l,P_l,P_r,rho_l,rho_r,vx_l,vx_r,vy_l,vy_r,r_l,r_r)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active

        u_prim(ivx,ix,iy) = u_prim(ivx,ix,iy) + dt*(qh(iP,ix,iy)+qh(irho,ix,iy)*qh(ivy,ix,iy)**2)/radii(ix,iy)
        u_prim(ivy,ix,iy) = u_prim(ivy,ix,iy) - dt*(qh(irho,ix,iy)*qh(ivy,ix,iy)*qh(ivx,ix,iy))/radii(ix,iy)

#if NDUST>0
         do idust = 1, ndust
            u_prim(ivdx(idust),ix,iy)   = u_prim(ivdx(idust),ix,iy) + dt*(qh(irhod(idust),ix,iy)*qh(ivdy(idust),ix,iy)**2/radii(ix,iy))
            u_prim(ivdy(idust),ix,iy)   = u_prim(ivdy(idust),ix,iy)-dt*(qh(irhod(idust),ix,iy)*qh(ivdy(idust),ix,iy)*qh(ivdx(idust),ix,iy)/radii(ix,iy))
         end do
#endif

      end do
   end do
#endif

#if GEOM==3
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust,ix,iy,cot_th,rr,dAx,dAy)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         rr     = radii(ix,iy)
         cot_th = cos(phi(ix,iy))/sin(phi(ix,iy))
         dAx    = (surf(ix+1,iy,1) - surf(ix,iy,1))/vol(ix,iy)      ! -> 2/r
         dAy    = (surf(ix,iy+1,2) - surf(ix,iy,2))/vol(ix,iy)      ! -> cot(theta)/r
         u_prim(ivx,ix,iy) = u_prim(ivx,ix,iy) + dt*( qh(iP,ix,iy)*dAx + qh(irho,ix,iy)*(qh(ivy,ix,iy)**2 + qh(ivz,ix,iy)**2)/rr )
         u_prim(ivy,ix,iy) = u_prim(ivy,ix,iy) + dt*( qh(iP,ix,iy)*dAy + qh(irho,ix,iy)*(qh(ivz,ix,iy)**2*cot_th - qh(ivx,ix,iy)*qh(ivy,ix,iy))/rr )
#if NDUST>0
         do idust = 1, ndust
            u_prim(ivdx(idust),ix,iy) = u_prim(ivdx(idust),ix,iy) + dt*qh(irhod(idust),ix,iy)*(qh(ivdy(idust),ix,iy)**2 + qh(ivdz(idust),ix,iy)**2)/rr
            u_prim(ivdy(idust),ix,iy) = u_prim(ivdy(idust),ix,iy) + dt*qh(irhod(idust),ix,iy)*(qh(ivdz(idust),ix,iy)**2*cot_th - qh(ivdx(idust),ix,iy)*qh(ivdy(idust),ix,iy))/rr
         end do
#endif
      end do
   end do
#endif



end subroutine Source_terms


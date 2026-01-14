subroutine Source_terms
   use parameters
   use commons
   use OMP_LIB
   use units
   implicit none
   integer :: idust, ix,iy,ivar
   real(dp) :: cs_eos, barotrop,S_r,S_l,P_l,P_r,rho_l,rho_r,vx_l,vx_r,vy_l,vy_r,r_l,r_r

   if (static) then
      return
   end if

#if GEOM==2
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy,S_r,S_l,P_l,P_r,rho_l,rho_r,vx_l,vx_r,vy_l,vy_r,r_l,r_r)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active

        u_prim(ivx,ix,iy) = u_prim(ivx,ix,iy) + dt*(q(iP,ix,iy)+q(irho,ix,iy)*q(ivy,ix,iy)**2)/radii(ix,iy)
        u_prim(ivy,ix,iy) = u_prim(ivy,ix,iy) - dt*(q(irho,ix,iy)*q(ivy,ix,iy)*q(ivx,ix,iy))/radii(ix,iy)

#if NDUST>0
         do idust = 1, ndust
            u_prim(ivdx(idust),ix,iy)   = u_prim(ivdx(idust),ix,iy) + dt*(q(irhod(idust),ix,iy)*q(ivdy(idust),ix,iy)**2/radii(ix,iy))
            u_prim(ivdy(idust),ix,iy)   = u_prim(ivdy(idust),ix,iy)-dt*(q(irhod(idust),ix,iy)*q(ivdy(idust),ix,iy)*q(ivdx(idust),ix,iy)/radii(ix,iy))
         end do
#endif

      end do
   end do
#endif




end subroutine Source_terms


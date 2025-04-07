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


        ! P_l   = half*(qm_x(iP,ix,iy)+qp_x(iP,ix,iy))
        ! P_r   = half*(qm_x(iP,ix,iy)+qp_x(iP,ix,iy))
        ! rho_l = half*(qm_x(irho,ix,iy)+qp_x(irho,ix,iy))
        ! rho_r = half*(qm_x(irho,ix,iy)+qp_x(irho,ix,iy))

        ! vx_l = half*(qm_x(ivx,ix,iy)+qp_x(ivx,ix,iy))
        ! vx_r = half*(qm_x(ivx,ix,iy)+qp_x(ivx,ix,iy))

        ! vy_l = half*(qm_x(ivy,ix,iy)+qp_x(ivy,ix,iy))
        ! vy_r = half*(qm_x(ivy,ix,iy)+qp_x(ivy,ix,iy))

        ! S_l = P_l+rho_l*vy_l**2
        ! S_r = P_r+rho_r*vy_r**2

        ! r_l = radii(ix,iy)-0.5d0*dx(ix,iy,1)
        ! r_r = radii(ix,iy)+0.5d0*dx(ix,iy,1)

        ! u_prim(ivx,ix,iy) = u_prim(ivx,ix,iy) + dt*half*(r_l*S_l+S_r*r_r)/radii(ix,iy)**2

        ! !S_l = -rho_l*vy_l*vx_l
        ! !S_r = -rho_r*vy_r*vx_r

        ! !u_prim(ivy,ix,iy) = u_prim(ivy,ix,iy)+dt*half*(S_l*r_l+S_r*r_r)/radii(ix,iy)**2


        ! P_l   = half*(qm_y(iP,ix,iy)+qp_y(iP,ix,iy))
        ! P_r   = half*(qm_y(iP,ix,iy)+qp_y(iP,ix,iy))
        ! rho_l = half*(qm_y(irho,ix,iy)+qp_y(irho,ix,iy))
        ! rho_r = half*(qm_y(irho,ix,iy)+qp_y(irho,ix,iy))

        ! vx_l = half*(qm_y(ivx,ix,iy)+qp_y(ivx,ix,iy))
        ! vx_r = half*(qm_y(ivx,ix,iy)+qp_y(ivx,ix,iy))

        ! vy_l = half*(qm_y(ivy,ix,iy)+qp_y(ivy,ix,iy))
        ! vy_r = half*(qm_y(ivy,ix,iy)+qp_y(ivy,ix,iy))

        ! !u_prim(ivx,ix,iy) = u_prim(ivx,ix,iy) + dt*half*(r_l*S_l+S_r*r_r)/radii(ix,iy)**2

        ! S_l = -rho_l*vy_l*vx_l
        ! S_r = -rho_r*vy_r*vx_r

        ! u_prim(ivy,ix,iy) = u_prim(ivy,ix,iy)+dt*half*(S_l+S_r)/radii(ix,iy)


        u_prim(ivx,ix,iy) = u_prim(ivx,ix,iy) + dt*(q(iP,ix,iy)+q(irho,ix,iy)*q(ivy,ix,iy)**2)/radii(ix,iy)
        u_prim(ivy,ix,iy) = u_prim(ivy,ix,iy) - dt*(q(irho,ix,iy)*q(ivy,ix,iy)*q(ivx,ix,iy))/radii(ix,iy)

#if NDUST>0
         print *, 'Careful source terms for the dust not correctly implemented'
         stop

         do idust = 1, ndust
            u_prim(ivdx(idust),ix,iy)   = u_prim(ivdx(idust),ix,iy) + dt*(q(irhod(idust),ix,iy)*q(ivdy(idust),ix,iy)**2/radii(ix,iy))
            u_prim(ivdy(idust),ix,iy)   = u_prim(ivdy(idust),ix,iy)-dt*(q(irhod(idust),ix,iy)*q(ivdy(idust),ix,iy)*q(ivdx(idust),ix,iy)/radii(ix,iy))
         end do
#endif

      end do
   end do
#endif

#if GEOM==4
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active

        ! vr is still vx but vphi is now vz !
        !No source term to vy (which is vz here)
        u_prim(ivx,ix,iy) = u_prim(ivx,ix,iy) + dt*(q(irho,ix,iy)*q(ivz,ix,iy)**2 + q(iP,ix,iy))/radii(ix,iy)
        u_prim(ivz,ix,iy) = u_prim(ivz,ix,iy) - dt*(q(irho,ix,iy)*q(ivz,ix,iy)*q(ivx,ix,iy)/radii(ix,iy))

#if NDUST>0
         do idust = 1, ndust
            u_prim(ivdx(idust),ix,iy) = u_prim(ivdx(idust),ix,iy) + dt*(q(irhod(idust),ix,iy)*q(ivdz(idust),ix,iy)**2/radii(ix,iy))
            u_prim(ivdz(idust),ix,iy)=u_prim(ivdz(idust),ix,iy)-dt*(q(irhod(idust),ix,iy)*q(ivdz(idust),ix,iy)*q(ivdx(idust),ix,iy)/radii(ix,iy))
         end do
#endif
      end do
   end do
#endif



end subroutine Source_terms


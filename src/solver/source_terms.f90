subroutine Source_terms
   use parameters
   use commons
   use OMP_LIB
   use units
   implicit none
   integer :: idust, ix,iy,ivar
   real(dp) :: cs_eos, barotrop

   if (static) then
      return
   end if

#if GEOM==2
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
        u_prim(ivx,ix,iy) = u_prim(ivx,ix,iy) + dt*(q(irho,ix,iy)*q(ivy,ix,iy)**2 + q(iP,ix,iy))/radii(ix,iy)
        u_prim(ivy,ix,iy) = u_prim(ivy,ix,iy) - dt*(q(irho,ix,iy)*q(ivy,ix,iy)*q(ivx,ix,iy)/radii(ix,iy))

#if NDUST>0
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


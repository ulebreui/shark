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
   !$omp parallel do default(shared) private(idust, ix,iy)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
        u_prim(ix,iy,ivx) = u_prim(ix,iy,ivx) + dt*(q(ix,iy,irho)*q(ix,iy,ivy)**2 + q(ix,iy,iP))/radii(ix,iy)
        u_prim(ix,iy,ivy) = u_prim(ix,iy,ivy) - dt*(q(ix,iy,irho)*q(ix,iy,ivy)*q(ix,iy,ivx)/radii(ix,iy))

#if NDUST>0
         do idust = 1, ndust
            u_prim(ix,iy,ivdx(idust)) = u_prim(ix,iy,ivdx(idust)) + dt*(q(ix,iy,irhod(idust))*q(ix,iy,ivdy(idust))**2/radii(ix,iy))
            u_prim(ix,iy,ivdy(idust))   = u_prim(ix,iy,ivdy(idust))-dt*(q(ix,iy,irhod(idust))*q(ix,iy,ivdy(idust))*q(ix,iy,ivdx(idust))/radii(ix,iy))
         end do
#endif

      end do
   end do
#endif

#if GEOM==4
   !$omp parallel do default(shared) private(idust, ix,iy)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active

        ! vr is still vx but vphi is now vz !
        !No source term to vy (which is vz here)
        u_prim(ix,iy,ivx) = u_prim(ix,iy,ivx) + dt*(q(ix,iy,irho)*q(ix,iy,ivz)**2 + q(ix,iy,iP))/radii(ix,iy)
        u_prim(ix,iy,ivz) = u_prim(ix,iy,ivz) - dt*(q(ix,iy,irho)*q(ix,iy,ivz)*q(ix,iy,ivx)/radii(ix,iy))

#if NDUST>0
         do idust = 1, ndust
            u_prim(ix,iy,ivdx(idust)) = u_prim(ix,iy,ivdx(idust)) + dt*(q(ix,iy,irhod(idust))*q(ix,iy,ivdz(idust))**2/radii(ix,iy))
            u_prim(ix,iy,ivdz(idust))=u_prim(ix,iy,ivdz(idust))-dt*(q(ix,iy,irhod(idust))*q(ix,iy,ivdz(idust))*q(ix,iy,ivdx(idust))/radii(ix,iy))
         end do
#endif
      end do
   end do
#endif



end subroutine Source_terms


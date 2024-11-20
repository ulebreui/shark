subroutine Source_terms
   use parameters
   use commons
   use OMP_LIB
   use units
   implicit none
   integer :: idust, ix, iy
   integer :: ixx, iyy, icell, ivar
   real(dp) :: ts, cs_eos, barotrop
   real(dp), dimension(:, :, :), allocatable :: S_U

   if (static) then
      return
   end if

   allocate (S_U(1:nx_max, 1:ny_max, 1:nvar))
   S_U = 0.0d0

   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
#if GEOM==2
        S_U(ix, iy, ivx) = S_U(ix, iy, ivx) + dt*(q(ix, iy, irho)*q(ix, iy, ivy)**2 + q(ix, iy, iP))/radii(ix, iy)
        S_U(ix, iy, ivy) = S_U(ix, iy, ivy) - dt*(q(ix, iy, irho)*q(ix, iy, ivy)*q(ix, iy, ivx)/radii(ix, iy))

#if NDUST>0
         do idust = 1, ndust
            S_U(ix, iy, ivdx(idust)) = S_U(ix, iy, ivdx(idust)) + dt*(q(ix, iy, irhod(idust))*q(ix, iy, ivdy(idust))**2/radii(ix, iy))
            S_U(ix,iy,ivdy(idust))=S_U(ix,iy,ivdy(idust))-dt*(q(ix,iy,irhod(idust))*q(ix,iy,ivdy(idust))*q(ix,iy,ivdx(idust))/radii(ix,iy))
         end do
#endif
#endif
#if GEOM==4
        ! vr is still vx but vphi is now vz !
        !No source term to vy (which is vz here)
        S_U(ix, iy, ivx) = S_U(ix, iy, ivx) + dt*(q(ix, iy, irho)*q(ix, iy, ivz)**2 + q(ix, iy, iP))/radii(ix, iy)
        S_U(ix, iy, ivz) = S_U(ix, iy, ivz) - dt*(q(ix, iy, irho)*q(ix, iy, ivz)*q(ix, iy, ivx)/radii(ix, iy))

#if NDUST>0
         do idust = 1, ndust
            S_U(ix, iy, ivdx(idust)) = S_U(ix, iy, ivdx(idust)) + dt*(q(ix, iy, irhod(idust))*q(ix, iy, ivdz(idust))**2/radii(ix, iy))
            S_U(ix,iy,ivdz(idust))=S_U(ix,iy,ivdz(idust))-dt*(q(ix,iy,irhod(idust))*q(ix,iy,ivdz(idust))*q(ix,iy,ivdx(idust))/radii(ix,iy))
         end do
#endif
#endif
      end do
   end do

   !Update state vector
   u_prim = u_prim + S_U
   deallocate (S_U)

end subroutine Source_terms


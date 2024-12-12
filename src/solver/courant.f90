!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine computes the timestep according to the CFL conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine courant
   use parameters
   use commons
   use units
   implicit none

   integer  :: idust, ix,iy
   real(dp) :: vmax, dxx, force_max, vv,fratio

   if (static) then
      return
   end if

   dt = 2d44
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy,vmax, dxx, force_max, vv, fratio) reduction(min: dt)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         !Cas 1D
         dxx  = min(dx(ix,iy,1), radii(ix,iy)*dx(ix,iy,2))
         vv   = abs(q(ivx,ix,iy)) + abs(q(ivy,ix,iy)) + abs(q(ivz,ix,iy))

         vmax = cs(ix,iy) + vv
#if NDUST>0
         do idust = 1, ndust
            vmax = max(vmax, abs(q(ivdx(idust),ix,iy)) + abs(q(ivdy(idust),ix,iy)) + abs(q(ivdz(idust),ix,iy)))
         end do
#endif
         !print(vmax)
         dt = min(dt, CFL*dxx/abs(vmax))
      end do
   end do

if (force_kick) then
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy,vmax, dxx, force_max, vv, fratio) reduction(min: dt)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         !Cas 1D
         dxx  = min(dx(ix,iy,1), radii(ix,iy)*dx(ix,iy,2))
         vv   = abs(q(ivx,ix,iy)) + abs(q(ivy,ix,iy)) + abs(q(ivz,ix,iy))+1d-30

         force_max = sqrt(force_x(ix,iy)**2 + force_y(ix,iy)**2 + force_z(ix,iy)**2)

         fratio = max(force_max*dxx/vv**2, 1d-3)
         dt = min(dt, CFL*dxx/vv*(sqrt(1.0d0 + 2.0d0*CFL*fratio) - 1.0d0)/fratio)
#if NDUST>0
            do idust = 1, ndust

               force_max = sqrt(force_dust_x(idust,ix,iy)**2 + force_dust_y(idust,ix,iy)**2 + force_dust_z(idust,ix,iy)**2)
               vv = abs(q(ivdx(idust),ix,iy)) + abs(q(ivdy(idust),ix,iy)) + abs(q(ivdz(idust),ix,iy))+1d-30
               fratio = max(force_max*dxx/vv**2, 1d-3)
               dt = min(dt, CFL*dxx/vv*(sqrt(1.0d0 + 2.0d0*CFL*fratio) - 1.0d0)/fratio)
            end do
#endif

      end do
   end do
end if

end subroutine courant


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

   integer  :: idust, ix,iy,ivar
   real(dp) :: vmax, dxx, force_max, vv,fratio,dt_loc
   real(dp), dimension(1:nvar) :: qloc,floc
   if (static) then
      return
   end if

   dt = 2d44
   if (force_kick) then

   !$omp parallel do default(shared) schedule(static) private(idust, ix,iy,vmax, dxx,force_max, vv, fratio,dt_loc,qloc,ivar) reduction(min: dt)
   do iy = first_active_y, last_active_y
      !$omp simd
      do ix = first_active, last_active
         !Cas 1D
         dxx  = min(dx(ix,iy,1), radii(ix,iy)*dx(ix,iy,2))

         do ivar =1,nvar
            qloc(ivar) = q(ivar,ix,iy)
         end do

         vv   = sqrt((qloc(ivx))**2 + (qloc(ivy))**2 + (qloc(ivz))**2)

         vmax = cs(ix,iy) + vv

#if NDUST>0
         do idust = 1, ndust
            vmax = max(vmax, sqrt(qloc(ivdx(idust))**2 + qloc(ivdy(idust))**2 + qloc(ivdz(idust))**2))
         end do
#endif
         dt_loc = 2d44
         !print(vmax)
         dt_loc = min(dt_loc, CFL*dxx/abs(vmax))

         
         force_max = sqrt(force_x(ix,iy)**2 + force_y(ix,iy)**2 + force_z(ix,iy)**2)
         fratio = max(force_max*dxx/vv**2, 1d-3)
         dt_loc = min(dt_loc, CFL*dxx/vv*(sqrt(1.0d0 + 2.0d0*CFL*fratio) - 1.0d0)/fratio)
#if NDUST>0
         do idust = 1, ndust

            force_max = sqrt(force_dust_x(idust,ix,iy)**2 + force_dust_y(idust,ix,iy)**2 + force_dust_z(idust,ix,iy)**2)
            vv = sqrt(qloc(ivdx(idust))**2 + qloc(ivdy(idust))**2 + qloc(ivdz(idust))**2)
            fratio = max(force_max*dxx/vv**2, 1d-3)
            dt_loc = min(dt_loc, CFL*dxx/vv*(sqrt(1.0d0 + 2.0d0*CFL*fratio) - 1.0d0)/fratio)
         end do
#endif
         dt = min(dt, dt_loc)
      end do
   end do

   else



   !$omp parallel do default(shared) schedule(static) private(idust, ix,iy,vmax, dxx,force_max, vv, fratio,dt_loc,qloc,ivar) reduction(min: dt)
   do iy = first_active_y, last_active_y
      !$omp simd
      do ix = first_active, last_active
         !Cas 1D
         dxx  = min(dx(ix,iy,1), radii(ix,iy)*dx(ix,iy,2))

         do ivar =1,nvar
            qloc(ivar) = q(ivar,ix,iy)
         end do

         vv   = sqrt((qloc(ivx))**2 + (qloc(ivy))**2 + (qloc(ivz))**2)

         vmax = cs(ix,iy) + vv

#if NDUST>0
         do idust = 1, ndust
            vmax = max(vmax, sqrt(qloc(ivdx(idust))**2 + qloc(ivdy(idust))**2 + qloc(ivdz(idust))**2))
         end do
#endif
         dt_loc = 2d44
         !print(vmax)
         dt_loc = min(dt_loc, CFL*dxx/abs(vmax))
         dt = min(dt, dt_loc)
      end do
   end do


   endif


end subroutine courant


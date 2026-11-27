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

   integer  :: idust, ix, iy
   real(dp) :: dr, dl, inv_dt, dt_loc, u, v, vb
   real(dp), parameter :: tiny_v = 1d-30

   if (static) return

   dt = 2d44
   !$omp parallel do default(shared) schedule(static) &
   !$omp private(idust,ix,iy,dr,dl,inv_dt,dt_loc,u,v,vb) reduction(min:dt)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         dr = dx(ix,iy,1)                   ! radial (or x) width
         dl = radii(ix,iy)*dx(ix,iy,2)      ! azimuthal (or y) width, radii=1 in Cartesian
         vb = 0.0d0                         ! FARGO: azimuthal velocity handled by the shift
#if GEOM==2
         if (fargo) vb = vbar(ix)
#endif
         ! Unsplit CFL: sum of the inverse crossing times over the directions
         u = abs(q(ivx,ix,iy)); v = abs(q(ivy,ix,iy) - vb)
         inv_dt = (u + cs(ix,iy))/dr
#if NY>1
         inv_dt = inv_dt + (v + cs(ix,iy))/dl
#endif
#if NDUST>0
         do idust = 1, ndust
            u = abs(q(ivdx(idust),ix,iy)); v = abs(q(ivdy(idust),ix,iy) - vb)
#if NY>1
            inv_dt = max(inv_dt, u/dr + v/dl)
#else
            inv_dt = max(inv_dt, u/dr)
#endif
         end do
#endif
         dt_loc = CFL/inv_dt

         ! Force constraint, per direction: |v| dt + |F| dt^2/2 = CFL dx
         if (force_kick) then
            dt_loc = min(dt_loc, dt_force(abs(q(ivx,ix,iy)), abs(force_x(ix,iy)), CFL*dr))
#if NY>1
            dt_loc = min(dt_loc, dt_force(abs(q(ivy,ix,iy) - vb), abs(force_y(ix,iy)), CFL*dl))
#endif
#if NDUST>0
            do idust = 1, ndust
               dt_loc = min(dt_loc, dt_force(abs(q(ivdx(idust),ix,iy)), abs(force_dust_x(idust,ix,iy)), CFL*dr))
#if NY>1
               dt_loc = min(dt_loc, dt_force(abs(q(ivdy(idust),ix,iy) - vb), abs(force_dust_y(idust,ix,iy)), CFL*dl))
#endif
            end do
#endif
         end if

#if GEOM==2
         ! FARGO: shear between adjacent rings (Masset 2000)
         if (fargo) dt_loc = min(dt_loc, CFL*dl/max(abs(vbar(ix+1)-vbar(ix)), abs(vbar(ix)-vbar(ix-1)), tiny_v))
#endif
         dt = min(dt, dt_loc)
      end do
   end do
   !$omp end parallel do

   dt = min(dt, tend - time)   ! land exactly on tend

contains

   ! Positive root of v*t + F*t**2/2 = L, written without cancellation (valid for F=0 or v=0)
   pure real(dp) function dt_force(vel, force, L)
      real(dp), intent(in) :: vel, force, L
      dt_force = 2.0d0*L/max(vel + sqrt(vel**2 + 2.0d0*force*L), tiny_v)
   end function dt_force

end subroutine courant
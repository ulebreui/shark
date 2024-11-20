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

   integer :: i, idust, ix, iy, icell
   real(dp) :: vmax, dxx, force_max, ca, magnetosonic_fast, vv, fratio, D_max

   if (static) then
      return
   end if

   dt = 2d44

   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         i = icell(ix, iy)
         !Cas 1D
         dxx = min(dx(ix, iy, 1), radii(ix, iy)*dx(ix, iy, 2))

         vmax = cs(ix, iy) + abs(q(ix, iy, ivx)) + abs(q(ix, iy, ivy)) + abs(q(ix, iy, ivz))

#if NDUST>0
         do idust = 1, ndust
            vmax = max(vmax, abs(q(ix, iy, ivdx(idust))) + abs(q(ix, iy, ivdy(idust))) + abs(q(ix, iy, ivdz(idust))))
         end do
#endif
         !print(vmax)
         dt = min(dt, CFL*dxx/abs(vmax))
         if (force_kick) then

            force_max = sqrt(force_x(ix, iy)**2 + force_y(ix, iy)**2 + force_z(ix, iy)**2)

#if NDUST>0
            do idust = 1, ndust

               force_max = max(force_max, sqrt(force_dust_x(ix, iy, idust)**2 + force_dust_y(ix, iy, idust)**2 + force_dust_z(ix, iy, idust)**2))

            end do
#endif

            if (vv .ne. 0.0d0) then
               fratio = max(force_max*dxx/vv**2, 1d-3)
               dt = min(dt, CFL*dxx/vv*(sqrt(1.0d0 + 2.0d0*CFL*fratio) - 1.0d0)/fratio)
            end if
         end if
      end do
   end do
end subroutine courant


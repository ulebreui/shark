
! This routine performs the predictor operation !
! Reference for cylindrical geom :
! The Athena Astrophysical MHD Code in Cylindrical Geometry
! A Proper Discretization of Hydrodynamic Equations in Cylindrical Coordinates for Astrophysical Simulations
subroutine predictor
   use parameters
   use commons
   use units
   use OMP_LIB
   use slope_limiter
   implicit none
   integer  :: ix,iy,idust
   integer  :: ivar, idim
   real(dp) :: barotrop, cs_eos, dx_loc,dx_l,dx_r,dy_l,dy_r
   real(dp) :: drx,dry,dpx,dpy,dux,duy,dvx,dvy,dwx,dwy,r_rho,u,v,w,p,sr0,sp0,su0,sv0,sw0,radius_polar
   integer  :: irho_spe, ivx_spe, ivy_spe, ivz_spe, ipscal


   if (static) return

   ! Initialise to zero



   !$omp parallel do default(shared) schedule(RUNTIME) private(ix,iy,idust,ivar, idim,dx_l,dx_r,dy_l,dy_r, dx_loc,drx,dry,dpx,dpy,dux,duy,dvx,dvy,dwx,dwy,r_rho,u,v,w,p,sr0,sp0,su0,sv0,sw0,radius_polar,irho_spe, ivx_spe, ivy_spe, ivz_spe, ipscal)
   do iy = 2, ny_max - 1
      do ix = 2, nx_max - 1
         radius_polar = radii(ix,iy)

         r_rho = q(irho,ix,iy)
         u = q(ivx,ix,iy)
         v = q(ivy,ix,iy)
         w = q(ivz,ix,iy)
         p = q(iP,ix,iy)
         dx_l = (dx(ix,iy,1)+dx(ix-1,iy,1))*half
         dx_r = (dx(ix+1,iy,1)+dx(ix,iy,1))*half
         dy_l = (dx(ix,iy,2)+dx(ix,iy-1,2))*half
         dy_r = (dx(ix,iy+1,2)+dx(ix,iy,2))*half

         drx = slope_limit((q(irho,ix,iy) - q(irho,ix-1,iy))/dx_l,(q(irho,ix+1,iy) - q(irho,ix,iy))/dx_r)
         dux = slope_limit((q(ivx,ix,iy)  - q(ivx,ix-1,iy)) /dx_l,(q(ivx,ix+1,iy)  - q(ivx,ix,iy))/dx_r)
         dvx = slope_limit((q(ivy,ix,iy)  - q(ivy,ix-1,iy)) /dx_l,(q(ivy,ix+1,iy)  - q(ivy,ix,iy))/dx_r)
         dwx = slope_limit((q(ivz,ix,iy)  - q(ivz,ix-1,iy)) /dx_l,(q(ivz,ix+1,iy)  - q(ivz,ix,iy))/dx_r)
         dPx = slope_limit((q(iP,ix,iy)   - q(iP,ix-1,iy))  /dx_l,(q(iP,ix+1,iy)   - q(iP,ix,iy))/dx_r)
         if (iso_cs == 1) then
            P = q(irho,ix,iy)*cs(ix,iy)**2
            dPx = slope_limit((q(irho,ix,iy) - q(irho,ix-1,iy))/dx_l,(q(irho,ix+1,iy) - q(irho,ix,iy))/dx_r)*cs(ix,iy)**2
         end if

         dry = slope_limit((q(irho,ix,iy) - q(irho,ix,iy-1))/dy_l,(q(irho,ix,iy+1) - q(irho,ix,iy))/dy_r)/radius_polar
         duy = slope_limit((q(ivx,ix,iy) - q(ivx,ix,iy-1))/dy_l,(q(ivx,ix,iy+1) - q(ivx,ix,iy))/dy_r)/radius_polar
         dvy = slope_limit((q(ivy,ix,iy) - q(ivy,ix,iy-1))/dy_l,(q(ivy,ix,iy+1) - q(ivy,ix,iy))/dy_r)/radius_polar
         dwy = slope_limit((q(ivz,ix,iy) - q(ivz,ix,iy-1))/dy_l,(q(ivz,ix,iy+1) - q(ivz,ix,iy))/dy_r)/radius_polar
         dpy = slope_limit((q(iP,ix,iy) - q(iP,ix,iy-1))/dy_l,(q(iP,ix,iy+1) - q(iP,ix,iy))/dy_r)/radius_polar

         if (iso_cs == 1) then
            dPy = slope_limit((q(irho,ix,iy) - q(irho,ix,iy-1))/dy_l,(q(irho,ix,iy+1) - q(irho,ix,iy))/dy_r)/radius_polar*cs(ix,iy)**2
         end if

         if (force_kick) then
            u = u + force_x(ix,iy)*half*dt
            v = v + force_y(ix,iy)*half*dt
            w = w + force_z(ix,iy)*half*dt
         end if

!Cartesian geometry
#if GE0M==0
         sr0 = -u*drx - v*dry - (dux + dvy)*r_rho
         sp0 = -u*dpx - v*dpy - (dux + dvy)*gamma*p
         su0 = -u*dux - v*duy - (dpx)/r_rho
         sv0 = -u*dvx - v*dvy - (dpy)/r_rho
         sw0 = -u*dwx - v*dwy
#endif

!Disk (face-on) geometry
#if GEOM==2
         sr0 = -u*drx - v*dry - (dux + dvy)*r_rho - r_rho*u/radius_polar
         sp0 = -u*dpx - v*dpy - (dux + dvy)*gamma*p - gamma*p*u/radius_polar
         su0 = -u*dux - v*duy - (dpx)/r_rho + (v**2.)/radius_polar
         sv0 = -u*dvx - v*dvy - (dpy)/r_rho - u*v/radius_polar
         sw0 = -u*dwx - v*dwy
#endif

!Disk (edge-on) geometry
#if GEOM==4
         sr0 = -u*drx - v*dry - (dux + dvy)*r_rho - r_rho*u/radii(ix,iy)
         sp0 = -u*dpx - v*dpy - (dux + dvy)*gamma*p - gamma*p*u/radii(ix,iy)
         su0 = -u*dux - v*duy - (dpx)/r_rho + (w**2.)/radii(ix,iy)
         sv0 = -u*dvx - v*dvy - (dpy)/r_rho
         sw0 = -u*dwx - v*dwy - u*w/radii(ix,iy)
#endif

         !direction x
         dx_loc = dx(ix,iy,1)
         qm_x(irho,ix,iy) = r_rho + half*dt*sr0 + half*drx*dx_loc
         qm_x(ivx,ix,iy) = u + half*dt*su0 + half*dux*dx_loc
         qm_x(iP,ix,iy) = max(p + half*dt*sP0 + half*dpx*dx_loc, smallP)
         qm_x(ivy,ix,iy) = v + half*dt*sv0 + half*dvx*dx_loc
         qm_x(ivz,ix,iy) = w + half*dt*sw0 + half*dwx*dx_loc
         qp_x(irho,ix,iy) = r_rho + half*dt*sr0 - half*drx*dx_loc
         qp_x(ivx,ix,iy) = u + half*dt*su0 - half*dux*dx_loc
         qp_x(iP,ix,iy) = max(p + half*dt*sP0 - half*dpx*dx_loc, smallP)
         qp_x(ivy,ix,iy) = v + half*dt*sv0 - half*dvx*dx_loc
         qp_x(ivz,ix,iy) = w + half*dt*sw0 - half*dwx*dx_loc

         ! direction y
         dx_loc = radius_polar*dx(ix,iy,2)
         qm_y(irho,ix,iy) = r_rho + half*dt*sr0 + half*dry*dx_loc
         qm_y(ivx,ix,iy) = u + half*dt*su0 + half*duy*dx_loc
         qm_y(iP,ix,iy) = max(p + half*dt*sP0 + half*dpy*dx_loc, smallP)
         qm_y(ivy,ix,iy) = v + half*dt*sv0 + half*dvy*dx_loc
         qm_y(ivz,ix,iy) = w + half*dt*sw0 + half*dwy*dx_loc
         qp_y(irho,ix,iy) = r_rho + half*dt*sr0 - half*dry*dx_loc
         qp_y(ivx,ix,iy) = u + half*dt*su0 - half*duy*dx_loc
         qp_y(iP,ix,iy) = max(p + half*dt*sP0 - half*dpy*dx_loc, smallP)
         qp_y(ivy,ix,iy) = v + half*dt*sv0 - half*dvy*dx_loc
         qp_y(ivz,ix,iy) = w + half*dt*sw0 - half*dwy*dx_loc

         ! Dust terms: same remark as for the gas
#if NDUST>0

         do idust = 1, ndust
            irho_spe = irhod(idust)
            ivx_spe = ivdx(idust)
            ivy_spe = ivdy(idust)
            ivz_spe = ivdz(idust)

            r_rho = q(irho_spe,ix,iy)
            u = q(ivx_spe,ix,iy)
            v = q(ivy_spe,ix,iy)
            w = q(ivz_spe,ix,iy)

            drx = slope_limit((q(irho_spe,ix,iy) - q(irho_spe,ix-1,iy))/dx_l,(q(irho_spe,ix+1,iy) - q(irho_spe,ix,iy))/dx_r)
            dux = slope_limit((q(ivx_spe,ix,iy) - q(ivx_spe,ix-1,iy))/dx_l,(q(ivx_spe,ix+1,iy) - q(ivx_spe,ix,iy))/dx_r)
            dvx = slope_limit((q(ivy_spe,ix,iy) - q(ivy_spe,ix-1,iy))/dx_l,(q(ivy_spe,ix+1,iy) - q(ivy_spe,ix,iy))/dx_r)
            dwx = slope_limit((q(ivz_spe,ix,iy) - q(ivz_spe,ix-1,iy))/dx_l,(q(ivz_spe,ix+1,iy) - q(ivz_spe,ix,iy))/dx_r)

            dry = slope_limit((q(irho_spe,ix,iy) - q(irho_spe,ix,iy-1))/dy_l,(q(irho_spe,ix,iy+1) - q(irho_spe,ix,iy))/dy_r)/radius_polar
            duy = slope_limit((q(ivx_spe,ix,iy) - q(ivx_spe,ix,iy-1))/dy_l,(q(ivx_spe,ix,iy+1) - q(ivx_spe,ix,iy))/dy_r)/radius_polar
            dvy = slope_limit((q(ivy_spe,ix,iy) - q(ivy_spe,ix,iy-1))/dy_l,(q(ivy_spe,ix,iy+1) - q(ivy_spe,ix,iy))/dy_r)/radius_polar
            dwy = slope_limit((q(ivz_spe,ix,iy) - q(ivz_spe,ix,iy-1))/dy_l,(q(ivz_spe,ix,iy+1) - q(ivz_spe,ix,iy))/dy_r)/radius_polar

            if (force_kick) then

               u = u + force_dust_x(idust,ix,iy)*half*dt
               v = v + force_dust_y(idust,ix,iy)*half*dt
               w = w + force_dust_z(idust,ix,iy)*half*dt

            end if

            sr0 = -u*drx - v*dry - (dux + dvy)*r_rho
            su0 = -u*dux - v*duy
            sv0 = -u*dvx - v*dvy
            sw0 = -u*dwx - v*dwy

#if GEOM==2
            !Polar geometry source terms
            sr0 = sr0 - r_rho*u/radius_polar
            su0 = su0 + v**2./radius_polar
            sv0 = sv0 - u*v/radius_polar
#endif

#if GEOM==4
            !Polar geometry source terms -- TODO add the missing source terms
            sr0 = sr0 - r_rho*u/radii(ix,iy)
            su0 = su0 + w**2./radii(ix,iy)
            sw0 = sw0 - u*w/radii(ix,iy)
#endif

            !Direction x
            dx_loc = dx(ix,iy,1)
            qm_x(irho_spe,ix,iy) = r_rho + half*dt*sr0 + half*drx*dx_loc
            qm_x(ivx_spe, ix,iy) = u + half*dt*su0 + half*dux*dx_loc
            qm_x(ivy_spe, ix,iy) = v + half*dt*sv0 + half*dvx*dx_loc
            qm_x(ivz_spe, ix,iy) = w + half*dt*sw0 + half*dwx*dx_loc
            qp_x(irho_spe,ix,iy) = r_rho + half*dt*sr0 - half*drx*dx_loc
            qp_x(ivx_spe, ix,iy) = u + half*dt*su0 - half*dux*dx_loc
            qp_x(ivy_spe, ix,iy) = v + half*dt*sv0 - half*dvx*dx_loc
            qp_x(ivz_spe, ix,iy) = w + half*dt*sw0 - half*dwx*dx_loc

            !direction y
            dx_loc = radius_polar*dx(ix,iy,2)
            qm_y(irho_spe,ix,iy) = r_rho + half*dt*sr0 + half*dry*dx_loc
            qm_y(ivx_spe, ix,iy) = u + half*dt*su0 + half*duy*dx_loc
            qm_y(ivy_spe, ix,iy) = v + half*dt*sv0 + half*dvy*dx_loc
            qm_y(ivz_spe, ix,iy) = w + half*dt*sw0 + half*dwy*dx_loc
            qp_y(irho_spe,ix,iy) = r_rho + half*dt*sr0 - half*dry*dx_loc
            qp_y(ivx_spe, ix,iy) = u + half*dt*su0 - half*duy*dx_loc
            qp_y(ivy_spe, ix,iy) = v + half*dt*sv0 - half*dvy*dx_loc
            qp_y(ivz_spe, ix,iy) = w + half*dt*sw0 - half*dwy*dx_loc

#if NDUSTPSCAL>0
            do ipscal = 1, ndustpscal

               r_rho = q(idust_pscal(idust,ipscal),ix,iy)
               drx   = slope_limit((q(idust_pscal(idust,ipscal),ix,iy) - q(idust_pscal(idust,ipscal),ix-1,iy))/dx_l,(q(idust_pscal(idust,ipscal),ix+1,iy) - q(idust_pscal(idust,ipscal),ix,iy))/dx_r)
               dry   = slope_limit((q(idust_pscal(idust,ipscal),ix,iy) - q(idust_pscal(idust,ipscal),ix,iy-1))/dy_l,(q(idust_pscal(idust,ipscal),ix,iy+1) - q(idust_pscal(idust,ipscal),ix,iy))/dy_r)/radius_polar

               sr0 = -u*drx - v*dry - (dux + dvy)*r_rho
#if GEOM==2
               !Polar geometry source terms
               sr0 = sr0 - r_rho*u/radius_polar
#endif
#if GEOM==4
               !Polar geometry source terms
               sr0 = sr0 - r_rho*u/radii(ix,iy)
#endif
               ! Direction x

               dx_loc = dx(ix,iy,1)
               qm_x(idust_pscal(idust, ipscal),ix,iy) = r_rho + half*dt*sr0 + half*drx*dx_loc
               qp_x(idust_pscal(idust, ipscal),ix,iy) = r_rho + half*dt*sr0 - half*drx*dx_loc

               ! Direction y
               dx_loc = radius_polar*dx(ix,iy,2)
               qm_y(idust_pscal(idust, ipscal),ix,iy) = r_rho + half*dt*sr0 + half*dry*dx_loc
               qp_y(idust_pscal(idust, ipscal),ix,iy) = r_rho + half*dt*sr0 - half*dry*dx_loc

            end do
#endif
         end do !dust loop
#endif
         ! We recompute the thermal pressure if we don't solve for the NRJ equation
         if (iso_cs == 1) then
               qp_x(iP, ix,iy) = cs(ix,iy)**2*qp_x(irho, ix,iy)
               qm_x(iP, ix,iy) = cs(ix,iy)**2*qm_x(irho, ix,iy)
               qp_y(iP, ix,iy) = cs(ix,iy)**2*qp_y(irho, ix,iy)
               qm_y(iP, ix,iy) = cs(ix,iy)**2*qm_y(irho, ix,iy)
         end if
         if (non_standard_eos == 1) then
               qp_x(iP, ix,iy) = cs_eos(barotrop(qp_x(irho, ix,iy)))**2*qp_x(irho, ix,iy)
               qm_x(iP, ix,iy) = cs_eos(barotrop(qm_x(irho, ix,iy)))**2*qm_x(irho, ix,iy)
               qp_y(iP, ix,iy) = cs_eos(barotrop(qp_y(irho, ix,iy)))**2*qp_y(irho, ix,iy)
               qm_y(iP, ix,iy) = cs_eos(barotrop(qm_y(irho, ix,iy)))**2*qm_y(irho, ix,iy)
         end if
      end do
   end do

end subroutine predictor

! This is the Riemmann solver : llf + Godunov scheme
subroutine add_delta_u
   use parameters
   use commons
   use units
   use OMP_LIB

   implicit none
   integer :: idust, ivar, ix,iy

   real(dp), dimension(1:nvar) :: qleft, qright, flx
   real(dp) :: csr, csl, barotrop, cs_eos

   if (static) return



   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ivar, ix,iy, qleft, qright, flx, csr, csl)
   do iy = 2, ny_max - 1
      do ix = 2, nx_max - 1

         !First direction x
         do ivar = 1, nvar
            qleft(ivar)  = qm_x(ivar,ix - 1, iy)
            qright(ivar) = qp_x(ivar,ix,iy)
            flx(ivar) = 0.0d0
         end do

         if (iso_cs == 1) then
            csl = cs(ix - 1, iy)
            csr = cs(ix,iy)
         else if (non_standard_eos == 1) then
            csl = cs_eos(barotrop(qleft(irho)))
            csr = cs_eos(barotrop(qright(irho)))
         else
            csl = sqrt(gamma*qleft(iP)/qleft(irho))
            csr = sqrt(gamma*qright(iP)/qright(irho))
         end if

         call solve_wrapper(qleft, qright, flx, csl, csr, 1)

         do ivar = 1, nvar
            flux_x(ivar,ix,iy) = flx(ivar)
         end do

         ! Then direction y

         do ivar = 1, nvar
            qleft(ivar)  = qm_y(ivar,ix,iy - 1)
            qright(ivar) = qp_y(ivar,ix,iy)
            flx(ivar)    = 0.0d0
         end do

         if (iso_cs == 1) then
            csl = cs(ix,iy - 1)
            csr = cs(ix,iy)
         else if (non_standard_eos == 1) then

            csl = cs_eos(barotrop(qleft(irho)))
            csr = cs_eos(barotrop(qright(irho)))

         else

            csl = sqrt(gamma*qleft(iP)/qleft(irho))
            csr = sqrt(gamma*qright(iP)/qright(irho))

         end if

         call solve_wrapper(qleft, qright, flx, csl, csr, 2)

         do ivar = 1, nvar
            flux_y(ivar,ix,iy) = flx(ivar)
         end do

      end do
   end do

   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ivar, ix,iy)
      do iy = first_active_y, last_active_y
         do ix = first_active, last_active
               do ivar = 1, nvar
                  u_prim(ivar,ix,iy)=u_prim(ivar,ix,iy) + (flux_x(ivar,ix,iy)*surf(ix,iy,1)-flux_x(ivar,ix+1,iy) *surf(ix+1,iy,1))  /vol(ix,iy)*dt&
&                             + (flux_y(ivar,ix,iy)*surf(ix,iy,2) - flux_y(ivar,ix,iy + 1)*surf(ix,iy + 1, 2))/vol(ix,iy)*dt
         end do
      end do
   end do
end subroutine add_delta_u

subroutine solve_wrapper(qleft, qright, flx, csl, csr, idim)
   use hydro_solvers
   use parameters
   use commons

   implicit none

   real(dp), dimension(1:nvar), intent(in) :: qright, qleft
   real(dp), dimension(1:nvar), intent(inout) :: flx
   real(dp) :: csl, csr

   integer  :: idim

   ! First the gas
#if SOLVER==0
   call solver_llf(qleft, qright, flx, csl, csr, idim)
#endif
#if SOLVER==1
   call solver_hll(qleft, qright, flx, csl, csr, idim)
#endif
#if SOLVER==2
   call solver_hllc(qleft, qright, flx, csl, csr, idim)
#endif

   ! Then the dust

#if NDUST>0

#if SOLVERDUST==0

   call solver_dust_Huang_Bai(qleft, qright, flx, idim)

#endif

#if SOLVERDUST==1
   call solver_dust_llf(qleft, qright, flx, idim)
#endif

#if SOLVERDUST==2
   call solver_dust_hll(qleft, qright, flx, idim)
#endif

#if SOLVERDUST==3
   call solver_hllc_dust(qleft, qright, flx, csl, csr, idim)
#endif
#endif

end subroutine solve_wrapper


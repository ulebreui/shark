!Get the index of a cell according to ix and iy
function icell(ix,iy)
   use parameters
   use commons
   use units
   implicit none
   integer  :: ix,iy
   integer  :: icell
#if NY==1
   icell = ix
#else
   icell = 1 + (ix - 1) + (iy - 1)*nx_max
#endif
end function icell

!Get the x index of a cell according to i
function ixx(i)
   use parameters
   use commons
   use units
   implicit none
   integer  :: i
   integer  :: ixx
#if NY==1
   ixx = i
#else
   ixx = mod(i - 1, nx_max) + 1
#endif

end function ixx

!Get the y index of a cell according to i
function iyy(i)
   use parameters
   use commons
   use units
   implicit none
   integer  :: i
   integer  :: iyy

   iyy = (i - 1)/nx_max + 1

end function iyy

!Grid initialisation routine
#if GEOM==0
subroutine gridinit(rmax_x, rmax_y)
#endif
#if GEOM>1
   subroutine gridinit(rmax_x, rmax_y, inner_r)
#endif
      use parameters
      use commons
      use units
      implicit none

      real(dp):: rmax_x, rmax_y, inner_r
      integer :: i, ix,iy,icell

      print *, 'Number of cells        = ', Ncells
      print *, 'Number of active cells = ', Ncells_active

#if GEOM==0
      print *, 'You are using a linear cartesian grid.'
      do ix = 1, nx_max
         do iy = 1, ny_max
            dx(ix,iy,1) = rmax_x/DBLE(nx)
            dx(ix,iy,2) = rmax_y/DBLE(ny)

            surf(ix,iy,1) = dx(ix,iy,1)
            surf(ix,iy,2) = dx(ix,iy,2)

            vol(ix,iy) = dx(ix,iy,1)*dx(ix,iy,2)
            position(ix,iy,1) = (DBLE(ix - first_active) + half)*rmax_x/DBLE(nx)
            position(ix,iy,2) = (DBLE(iy - first_active_y) + half)*rmax_y/DBLE(ny)
         end do
      end do
#endif

! Polar grid (face-on)
#if GEOM==2
      print *, 'You are using a linear cylindrical grid.'

      do ix = 1, nx_max
         do iy = 1, ny_max
            dx(ix,iy,1) = rmax_x/DBLE(nx) ! d_r
            dx(ix,iy,2) = rmax_y/DBLE(ny) ! d_Phi
            radii(ix,iy) = inner_r + (DBLE(ix - first_active) + half)*rmax_x/DBLE(nx)
            phi(ix,iy) = (DBLE(iy - first_active_y) + half)*rmax_y/DBLE(ny)

            position(ix,iy,1) = radii(ix,iy)*cos(phi(ix,iy))
            position(ix,iy,2) = radii(ix,iy)*sin(phi(ix,iy))
            surf(ix,iy,1) = (radii(ix, iy) - half*dx(ix,iy,1))*(dx(ix,iy,2))  ! r dpho
            surf(ix,iy,2) = dx(ix,iy,1)! dr
            vol(ix,iy) = radii(ix,iy)*dx(ix,iy,1)*dx(ix,iy,2)

         end do
      end do
#endif

   end subroutine gridinit

!Grid initialisation routine

#if GEOM==2
#if NY>1
   subroutine gridinit_disk_log(rmax_x, inner_r)
      use parameters
      use commons
      use units
      implicit none

      real(dp):: rmax_x, rmax_y, inner_r
      integer :: i, ix,iy,icell
#if GRIDSPACE==1
      !real(dp), dimension(1,nx+1):: radii_edges
      real(dp) :: zeta_r
#endif
      print *, 'Number of cells        =', Ncells
      print *, 'Number of active cells =', Ncells_active

! Polar grid (face-on)
#if GRIDSPACE==0
      print *, 'You are using a linear space cylindrical grid.'

      do ix = 1, nx_max
         do iy = 1, ny_max
            dx(ix,iy,1) = rmax_x/DBLE(nx) ! d_r
            dx(ix,iy,2) = 2.0d0*pi/DBLE(ny) ! d_Phi
            radii(ix,iy) = inner_r + (DBLE(ix - first_active) + half)*rmax_x/DBLE(nx)
            phi(ix,iy) = (DBLE(iy - first_active_y) + half)*2.0d0*pi/DBLE(ny)

            position(ix,iy,1) = radii(ix,iy)*cos(phi(ix,iy))
            position(ix,iy,2) = radii(ix,iy)*sin(phi(ix,iy))
            surf(ix,iy,1) = (radii(ix,iy) - half*dx(ix,iy,1))*(dx(ix,iy,2))  ! r dpho
            surf(ix,iy,2) = dx(ix,iy,1)! dr
            vol(ix,iy) = radii(ix,iy)*dx(ix,iy,1)*dx(ix,iy,2)

         end do
      end do
#endif
#if GRIDSPACE==1
      print *, 'You are using a log space cylindrical grid.'

      zeta_r = (rmax_x/(inner_r))**(1.0d0/(nx - 1))
      do iy = 1, ny_max
         do ix = 1, nx_max
            radii(ix,iy) = inner_r*half*(zeta_r**(ix - first_active) + zeta_r**(ix + 1 - first_active))
            dx(ix,iy,1) = inner_r*(zeta_r**(ix + 1 - first_active) - zeta_r**(ix - first_active))
            dx(ix,iy,2) = 2.0d0*pi/DBLE(ny) ! d_Phi
            phi(ix,iy) = (DBLE(iy - first_active_y) + half)*2.0d0*pi/DBLE(ny)
            position(ix,iy,1) = radii(ix,iy)*cos(phi(ix,iy))
            position(ix,iy,2) = radii(ix,iy)*sin(phi(ix,iy))
            surf(ix,iy,1) = inner_r*(zeta_r**(ix - first_active))*(dx(ix,iy,2))  ! r dphi
            surf(ix,iy,2) = dx(ix,iy,1)! dr
            vol(ix,iy) = radii(ix,iy)*dx(ix,iy,1)*dx(ix,iy,2)
         end do
         do ix = 1, nghost
            radii(ix,iy) = inner_r*0.5d0
            radii(nx_max + 1 - ix,iy) = rmax_x + 0.5d0*dx(last_active, iy,1)
            dx(ix,iy,1) = dx(first_active, iy,1)
            dx(nx_max + 1 - ix,iy,1) = dx(last_active, iy,1)
         end do
      end do
#endif

   end subroutine gridinit_disk_log
#endif

#endif

#if GEOM==4
#if NY>1
   subroutine gridinit_disk_log(rmax_x, inner_r, rmax_y)
      use parameters
      use commons
      use units
      implicit none

      real(dp):: rmax_x, rmax_y, inner_r
      integer :: i, ix,iy,icell
#if GRIDSPACE==1
      !real(dp), dimension(1,nx+1):: radii_edges
      real(dp) :: zeta_r
#endif
      print *, 'Number of cells        =', Ncells
      print *, 'Number of active cells =', Ncells_active

! Polar grid (face-on)
#if GRIDSPACE==0
      print *, 'You are using a linear space 2D cylindrical grid - R,Z configuration'

      do ix = 1, nx_max
         do iy = 1, ny_max
            dx(ix,iy,1) = rmax_x/DBLE(nx) ! d_r
            dx(ix,iy,2) = rmax_y/DBLE(ny) ! d_y
            radii(ix,iy) = inner_r + (DBLE(ix - first_active) + half)*rmax_x/DBLE(nx)
            position(ix,iy,1) = radii(ix,iy)
            position(ix,iy,2) = (DBLE(iy - first_active_y) + half)*rmax_y/DBLE(ny)
            surf(ix,iy,1) = (dx(ix,iy,2))  ! dy
            surf(ix,iy,2) = dx(ix,iy,1)! dr
            vol(ix,iy) = dx(ix,iy,1)*dx(ix,iy,2)

         end do
      end do
#endif
#if GRIDSPACE==1
      print *, 'You are using a log space 2D cylindrical grid - R,Z configuration'

      zeta_r = (rmax_x/(inner_r))**(1.0d0/(nx))
      do iy = 1, ny_max
         do ix = 1, nx_max
            radii(ix,iy) = inner_r*half*(zeta_r**(ix - first_active) + zeta_r**(ix + 1 - first_active))
            dx(ix,iy,1) = inner_r*(zeta_r**(ix + 1 - first_active) - zeta_r**(ix - first_active))
            dx(ix,iy,2) = rmax_y/DBLE(ny) ! d_y
            position(ix,iy,1) = radii(ix,iy)
            position(ix,iy,2) = (DBLE(iy - first_active_y) + half)*rmax_y/DBLE(ny)
            surf(ix,iy,1) = (dx(ix,iy,2))   ! dy
            surf(ix,iy,2) = dx(ix,iy,1)     ! dr
            vol(ix,iy) = dx(ix,iy,1)*dx(ix,iy,2)
         end do
         do ix = 1, nghost
            radii(ix,iy) = inner_r*0.5d0
            radii(nx_max + 1 - ix,iy)) = rmax_x + 0.5d0*dx(last_active, iy,1)
            dx(ix,iy,1) = dx(first_active, iy,1)
            dx(nx_max + 1 - ix,iy,1) = dx(last_active, iy,1)
         end do
      end do
#endif

   end subroutine gridinit_disk_log
#endif
#endif


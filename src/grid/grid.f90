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
#if GRIDSPACE==0
   subroutine gridinit_disk_log(rmax_x, inner_r,lsoft)
#endif
#if GRIDSPACE==1
   subroutine gridinit_disk_log(rmax_x, inner_r,lsoft)
#endif
#if GRIDSPACE==2
   subroutine gridinit_disk_log(rmax_x, inner_r,rcut,nxcut,lsoft)

#endif
      use parameters
      use commons
      use units
      implicit none

      real(dp):: rmax_x, rmax_y, inner_r,rplus,rminus,lsoft,xx,yy
      integer :: i, ix,iy,icell,ixx,iyy
      real(dp), dimension(1:nx_max+1):: radii_left
#if GRIDSPACE==1
      !real(dp), dimension(1,nx+1):: radii_edges
      real(dp) :: zeta_r
#endif
#if GRIDSPACE==2
      !real(dp), dimension(1,nx+1):: radii_edges
      real(dp) :: zeta_r
      real(dp) :: rcut
      integer :: nxcut
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
      radii_left=0.0d0
      zeta_r = (rmax_x/(inner_r))**(1.0d0/(nx+1)) 
      do ix = 1, nghost
         radii_left(ix) = inner_r/((zeta_r)**(nghost-ix+1))
      end do
      radii_left(1)=0.0d0
      radii_left(first_active)= inner_r
      do ix = first_active+1, nx_max+1
         radii_left(ix)= radii_left(ix-1)*zeta_r
      end do      
      do iy = 1, ny_max
         do ix = 1, nx_max
            radii(ix,iy) = 0.5d0*(radii_left(ix)+radii_left(ix+1))
            dx(ix,iy,1)  = (radii_left(ix+1)-radii_left(ix))
            dx(ix,iy,2)  = 2.0d0*pi/DBLE(ny) ! d_Phi
            phi(ix,iy) = (DBLE(iy - first_active_y)+half)*2.0d0*pi/DBLE(ny)
            position(ix,iy,1) = radii(ix,iy)*cos(phi(ix,iy))
            position(ix,iy,2) = radii(ix,iy)*sin(phi(ix,iy))
            surf(ix,iy,1) = radii_left(ix)*(dx(ix,iy,2))  ! r dphi
            surf(ix,iy,2) = dx(ix,iy,1)! dr
            vol(ix,iy) = pi/DBLE(ny)*(radii_left(ix+1)**2-radii_left(ix)**2)
         end do

      end do
      print*, radii_left
     ! We make sure of the azimuthal periodicity
     do ix = 1,nx_max
      do iy = 1,nghost   
         position(ix,iy,1)          = position(ix,last_active_y  - nghost+iy,1) 
         position(ix,ny_max+1-iy,1) = position(ix,first_active_y + nghost-iy,1)
         position(ix,iy,2)          = position(ix,last_active_y  - nghost+iy,2) 
         position(ix,ny_max+1-iy,2) = position(ix,first_active_y + nghost-iy,2)

         ! surf(ix,iy,1)          = surf(ix,last_active_y  - nghost+iy,1) 
         ! surf(ix,ny_max+1-iy,1) = surf(ix,first_active_y + nghost-iy,1)
         ! surf(ix,iy,2)          = surf(ix,last_active_y  - nghost+iy,2) 
         ! surf(ix,ny_max+1-iy,2) = surf(ix,first_active_y + nghost-iy,2)        

         ! radii(ix,iy)          = radii(ix,last_active_y  - nghost+iy) 
         ! radii(ix,ny_max+1-iy) = radii(ix,first_active_y + nghost-iy)

         phi(ix,iy)          = phi(ix,last_active_y  - nghost+iy) 
         phi(ix,ny_max+1-iy) = phi(ix,first_active_y + nghost-iy)
         end do
   end do
#endif

#if GRIDSPACE==2
      print *, 'You are using a lin-log space cylindrical grid.'
      radii_left=0.0d0

      do ix = 1,nxcut
         radii_left(ix)=DBLE(ix-1)*rcut/nxcut
      end do

      zeta_r = (rmax_x/(rcut))**(1.0d0/(nx-nxcut+1)) 
      do ix = nxcut+1, nx_max+1
         radii_left(ix)= radii_left(ix-1)*zeta_r
      end do      
      do iy = 1, ny_max
         do ix = 1, nx_max
            radii(ix,iy) = 0.5d0*(radii_left(ix)+radii_left(ix+1))
            dx(ix,iy,1)  = (radii_left(ix+1)-radii_left(ix))
            dx(ix,iy,2)  = 2.0d0*pi/DBLE(ny) ! d_Phi
            phi(ix,iy) = (DBLE(iy - first_active_y)+half)*2.0d0*pi/DBLE(ny)
            position(ix,iy,1) = radii(ix,iy)*cos(phi(ix,iy))
            position(ix,iy,2) = radii(ix,iy)*sin(phi(ix,iy))
            surf(ix,iy,1) = radii_left(ix)*(dx(ix,iy,2))  ! r dphi
            surf(ix,iy,2) = dx(ix,iy,1)! dr
            vol(ix,iy) = pi/DBLE(ny)*(radii_left(ix+1)**2-radii_left(ix)**2)
         end do

      end do
      print*, radii_left
     ! We make sure of the azimuthal periodicity
     do ix = 1,nx_max
      do iy = 1,nghost   
         position(ix,iy,1)          = position(ix,last_active_y  - nghost+iy,1) 
         position(ix,ny_max+1-iy,1) = position(ix,first_active_y + nghost-iy,1)
         position(ix,iy,2)          = position(ix,last_active_y  - nghost+iy,2) 
         position(ix,ny_max+1-iy,2) = position(ix,first_active_y + nghost-iy,2)

         ! surf(ix,iy,1)          = surf(ix,last_active_y  - nghost+iy,1) 
         ! surf(ix,ny_max+1-iy,1) = surf(ix,first_active_y + nghost-iy,1)
         ! surf(ix,iy,2)          = surf(ix,last_active_y  - nghost+iy,2) 
         ! surf(ix,ny_max+1-iy,2) = surf(ix,first_active_y + nghost-iy,2)        

         ! radii(ix,iy)          = radii(ix,last_active_y  - nghost+iy) 
         ! radii(ix,ny_max+1-iy) = radii(ix,first_active_y + nghost-iy)

         phi(ix,iy)          = phi(ix,last_active_y  - nghost+iy) 
         phi(ix,ny_max+1-iy) = phi(ix,first_active_y + nghost-iy)
         end do
   end do

#endif
   ! do iy = 1, ny_max
   !    do ix = 1, nx_max
   !             xx = position(ix,iy,1)
   !             yy = position(ix,iy,2)
   !             do iyy = 1, ny_max
   !                do ixx = 1, nx_max
   !                   distance(ix,iy,ixx,iyy)= sqrt((xx-position(ixx,iyy,1))**2+(yy-position(ixx,iyy,2))**2+radii(ix,iy)**2*lsoft**2)      
   !                end do 
   !             end do
   !          enddo
   !       enddo
   end subroutine gridinit_disk_log
#endif

#endif


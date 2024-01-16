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
  icell = 1+(ix-1)+(iy-1)*nx_max
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
  ixx = mod(i-1,nx_max)+1
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

  iyy=(i-1)/nx_max+1

end function iyy

!Get the x index of a cell according to position /!\ Works only for regular grid
function x_to_ix(xx)
  use parameters
  use commons
  use units
  implicit none
  integer   :: x_to_ix
  real(dp)  :: xx
  x_to_ix = 1 + floor(xx/dx(1,1))
end function x_to_ix

subroutine get_active_cells
  use parameters
  use commons
  use units
  implicit none
  integer  :: ix,iy,icell,ii,icount,ixx,iyy
  allocate(active_cell(1:ncells))
  allocate(active_cell_predictor(1:ncells))

  active_cell=0
  active_cell_predictor=0
  if(ndim==1) then
    icount=0
  	do ix=1,nx_max
  		if(ix.ge.first_active.and.ix.le.last_active) then 
      active_cell(ix)=1
      icount=icount+1
      endif
      if(ix.ge.2 .and.ix.le.nx_max-1) active_cell_predictor(ix)=1

  	end do
  else
    icount=0
  	do ix=1,nx_max
  		do iy=1,ny_max
  		ii=icell(ix,iy)
      !print *, ii, ix, iy,ixx(ii),iyy(ii)
  		if(ix.ge.first_active.and.ix.le.last_active) then
  			if(iy.ge.first_active_y.and.iy.le.last_active_y) then
  		 	  active_cell(ii)=1
          icount=icount+1
  		 	endif
  		  endif
        if(ix.ge.2 .and.ix.le.nx_max-1) then
          if(iy.ge.2 .and.iy.le.ny_max-1) then
            active_cell_predictor(ii)=1
          endif
        endif
  		end do
  	end do	
  endif
  print *, icount, ' active cells were prepared.'
end subroutine get_active_cells


!Grid initialisation routine
#if GEOM==0
subroutine gridinit(rmax_x,rmax_y)
#endif
#if GEOM==1
subroutine gridinit(rmax_x,rmax_y)
#endif
#if GEOM==2
subroutine gridinit(rmax_x,rmax_y,inner_r)
#endif
  use parameters
  use commons
  use units
  implicit none

  real(dp):: rmax_x,rmax_y,inner_r
  integer :: i, ix, iy,icell

  print *, 'Number of cells        = ', Ncells
  print *, 'Number of active cells = ', Ncells_active

#if GEOM==0
  print *,'You are using a linear cartesian grid.'

  if(ndim==1) then
    do ix = 1, ncells
      dx (ix,1)  = rmax_x/float(nx)
      vol(ix)    = dx(ix,1)
      surf(ix,1) = 1
      if(active_cell(ix)==1) position(ix,1)=position(ix-1,1)+dx (ix,1)
    end do
  else
      do ix = 1, nx_max
        do iy = 1, ny_max
          dx  (icell(ix,iy),1) = rmax_x/float(nx)
          dx  (icell(ix,iy),2) = rmax_y/float(ny)

          surf(icell(ix,iy),1) = dx(icell(ix,iy),1)
          surf(icell(ix,iy),2) = dx(icell(ix,iy),2)

          vol (icell(ix,iy))       = dx(icell(ix,iy),1)*dx(icell(ix,iy),2)
          position(icell(ix,iy),1) = (float(ix-first_active)+half)   * rmax_x/float(nx)
          position(icell(ix,iy),2) = (float(iy-first_active_y)+half) * rmax_y/float(ny)
      end do
    end do
  endif
#endif


! Polar grid (face-on)
#if GEOM==2
      print *,'You are using a linear cylindrical grid.'

      do ix = 1, nx_max
        do iy = 1, ny_max
          dx  (icell(ix,iy),1)     = rmax_x/float(nx) ! d_r
          dx  (icell(ix,iy),2)     = rmax_y/float(ny) ! d_Phi
          radii_c(icell(ix,iy))    = inner_r + (float(ix-first_active)+half)   * rmax_x/float(nx)
          phi(icell(ix,iy))        = (float(iy-first_active_y)+half) * rmax_y/float(ny)


          position(icell(ix,iy),1) = radii_c(icell(ix,iy))*cos(phi(icell(ix,iy)))
          position(icell(ix,iy),2) = radii_c(icell(ix,iy))*sin(phi(icell(ix,iy)))  
          surf(icell(ix,iy),1)     = (radii_c(icell(ix,iy))-half* dx (icell(ix,iy),1))*(dx  (icell(ix,iy),2))  ! r dpho
          surf(icell(ix,iy),2)     = dx  (icell(ix,iy),1)! dr
          vol (icell(ix,iy))       = radii_c(icell(ix,iy))*dx  (icell(ix,iy),1)*dx  (icell(ix,iy),2)

        end do
      end do
#endif

end subroutine gridinit

!Grid initialisation routine

#if GEOM==2
subroutine gridinit_disk_log(rmax_x,inner_r)
  use parameters
  use commons
  use units
  implicit none

  real(dp):: rmax_x,rmax_y,inner_r
  integer :: i, ix, iy,icell
#if GRIDSPACE==1
  !real(dp), dimension(1,nx+1):: radii_edges
  real(dp) :: zeta_r
#endif
  print *, 'Number of cells        =', Ncells
  print *, 'Number of active cells =', Ncells_active

! Polar grid (face-on)
#if GRIDSPACE==0
      print *,'You are using a linear space cylindrical grid.'

      do ix = 1, nx_max
        do iy = 1, ny_max
          dx  (icell(ix,iy),1)     = rmax_x/float(nx) ! d_r
          dx  (icell(ix,iy),2)     = 2.0d0*pi/float(ny) ! d_Phi
          radii_c(icell(ix,iy))    = inner_r + (float(ix-first_active)+half)   * rmax_x/float(nx)
          phi(icell(ix,iy))        = (float(iy-first_active_y)+half) * 2.0d0*pi/float(ny)


          position(icell(ix,iy),1) = radii_c(icell(ix,iy))*cos(phi(icell(ix,iy)))
          position(icell(ix,iy),2) = radii_c(icell(ix,iy))*sin(phi(icell(ix,iy)))  
          surf(icell(ix,iy),1)     = (radii_c(icell(ix,iy))-half* dx (icell(ix,iy),1))*(dx  (icell(ix,iy),2))  ! r dpho
          surf(icell(ix,iy),2)     = dx  (icell(ix,iy),1)! dr
          vol (icell(ix,iy))       = radii_c(icell(ix,iy))*dx  (icell(ix,iy),1)*dx  (icell(ix,iy),2)

        end do
      end do
#endif
#if GRIDSPACE==1
  print *,'You are using a log space cylindrical grid.'

  zeta_r=(rmax_x/(inner_r))**(1.0d0/(nx))
  do iy = 1, ny_max
    do ix = 1, nx_max
      radii_c(icell(ix,iy))    = inner_r*half*(zeta_r**(ix-first_active)+zeta_r**(ix+1-first_active))
      dx  (icell(ix,iy),1)     = inner_r*(zeta_r**(ix+1-first_active)-zeta_r**(ix-first_active))
      dx  (icell(ix,iy),2)     = 2.0d0*pi/float(ny) ! d_Phi
      phi(icell(ix,iy))        = (float(iy-first_active_y)+half) *2.0d0*pi/float(ny)
      position(icell(ix,iy),1) = radii_c(icell(ix,iy))*cos(phi(icell(ix,iy)))
      position(icell(ix,iy),2) = radii_c(icell(ix,iy))*sin(phi(icell(ix,iy))) 
      surf(icell(ix,iy),1)     = inner_r*(zeta_r**(ix-first_active))*(dx  (icell(ix,iy),2))  ! r dphi
      surf(icell(ix,iy),2)     = dx  (icell(ix,iy),1)! dr
      vol (icell(ix,iy))       = radii_c(icell(ix,iy))*dx  (icell(ix,iy),1)*dx  (icell(ix,iy),2)
    end do
    do ix=1,nghost
      radii_c(icell(ix,iy))             = inner_r*0.5d0
      radii_c(icell(nx_max+1-ix,iy))    = rmax_x+0.5d0*dx  (icell(last_active,iy),1)
      dx  (icell(ix,iy),1)              = dx  (icell(first_active,iy),1)
      dx  (icell(nx_max+1-ix,iy),1)     = dx  (icell(last_active,iy),1)
    end do
  end do
#endif

! #if GRIDSPACE==1
!   print *,'You are using a log space cylindrical grid.'

!   zeta_r=(rmax_x/inner_r)**(1.0d0/nx)
!   do iy = 1, ny_max
!     do ix = first_active, last_active
!       rplus =inner_r*zeta_r**(ix+1-first_active)
!       rminus=inner_r*zeta_r**(ix-first_active)
!       radii_c(icell(ix,iy))    = half*(rplus+rminus)
!       dx  (icell(ix,iy),1)     = rplus-rminus
!       dx  (icell(ix,iy),2)     = 2.0d0*pi/float(ny) ! d_Phi
!       phi(icell(ix,iy))        = (float(iy-first_active_y)+half) *2.0d0*pi/float(ny)
!       position(icell(ix,iy),1) = radii_c(icell(ix,iy))*cos(phi(icell(ix,iy)))
!       position(icell(ix,iy),2) = radii_c(icell(ix,iy))*sin(phi(icell(ix,iy))) 
!       surf(icell(ix,iy),1)     = rminus*(dx  (icell(ix,iy),2))  ! r dphi
!       surf(icell(ix,iy),2)     = dx  (icell(ix,iy),1)! dr
!       vol (icell(ix,iy))       = radii_c(icell(ix,iy))*dx  (icell(ix,iy),1)*dx  (icell(ix,iy),2)
!     end do
!     do ix=1,nghost
!       radii_c(icell(ix,iy))             = inner_r*0.5d0
!       radii_c(icell(nx_max+1-ix,iy))    = rmax_x+0.5d0*dx  (icell(last_active,iy),1)
!       dx  (icell(ix,iy),1)              = dx  (icell(first_active,iy),1)
!       dx  (icell(nx_max+1-ix,iy),1)     = dx  (icell(last_active,iy),1)
!     end do
!   end do
! #endif
end subroutine gridinit_disk_log
#endif

#if NX>1
#if GEOM==1
!Grid initialisation routine
subroutine gridinit_sphere1D(rmax)
  use parameters
  use commons
  use units
  implicit none

  real(dp):: rmax,zeta_r,ms,inner_r
  integer :: i, ix, iy,icell

  print *, 'Number of cells        =', Ncells
  print *, 'Number of active cells =', Ncells_active

  print *,'You are using a spherical grid in 1D'
  print *,' Logarithmic grid'


 !  zeta_r=(rmax/(rin*au))**(1.0d0/(ncells_active-1))

 !  radii   = 0.0d0
 !  radii_c = 0.0d0
 !  do i=first_active,last_active
 !     radii(i)= (rin*au/unit_l)*zeta_r**(i-first_active)
 !  enddo

 !  do i=first_active,last_active
 !     dx(i,1)=radii(i)-radii(i-1)
 !  end do

 !  do i=1,first_active
 !     dx(i,1)=dx(first_active,1)
 !  end do

 !  !Cell center
 !  do i=first_active,ncells
 !     radii_c(i)     = ( (radii(i)**3 + radii(i-1)**3) / 2.)**(1./3.)
 !  end do

 ! !Cell volume
 !  do i=first_active,last_active
 !     vol(i) = (radii(i)**3.-radii(i-1)**3.)/3.0d0
 !  end do

 !  !Cell volume from center to the right 
 !  do i=first_active,last_active
 !     dvol(i)  =(radii(i)**3-(radii_c(i)**3) )/3.0d0
 !   end do

 !  do i=first_active,last_active
 !     Surf(i,1) = radii(i-1)**2.
 !  end do

 !  !Surf(first_active-1,1) = 0.0
 !  !Surf(last_active+1,1)  = 0.0
 !  do i=1,ncells
 !  position(i,1)=radii(i)
 !  end do

 ! do i=first_active,ncells-1!-nghost
 !     dx_c(i)      = radii_c(i+1) - radii_c(i-1)
 !     dx_r(i)      = radii_c(i+1) - radii_c(i)
 !     dx_l(i)      = radii_c(i)   - radii_c(i-1)
 !     dx_r_cell(i) = radii_c(i+1) - radii(i)
 !     dx_l_cell(i) = radii(i)     - radii_c(i)
 !  end do
 !  do i=1,nghost
 !     dx_l(i)      = radii_c(first_active)
 !     dx_l_cell(i) = radii_c(first_active)
 !     dx_r(i)      = dx_r(first_active)
 !     dx_r_cell(i) = dx_r_cell(first_active)
 !  end do
  
  
    zeta_r =(rmax/(rin*au))**(1.0d0/(nx-1))
    inner_r= rin*au/unit_l
  
    do i = first_active,last_active
      radii(i)= (rin*au/unit_l)*zeta_r**(i-first_active)
    end do
    do i=first_active,last_active
     dx(i,1)=radii(i)-radii(i-1)
    end do
    do i=1,first_active
     dx(i,1)=dx(first_active,1)
    end do
    !Cell center
    do i=first_active,ncells
       radii_c(i)     = ( (radii(i)+ radii(i-1)) / 2.)
    end do
    !Cell volume
    do i=first_active,last_active
       vol(i) = (radii(i)**3.-radii(i-1)**3.)/3.0d0
    end do

    !Cell volume from center to the right 
    do i=first_active,last_active
       dvol(i)  =(radii(i)**3-(radii_c(i)**3) )/3.0d0
     end do

    do i=first_active,last_active
       Surf(i,1) = radii(i-1)**2.
    end do
    !Surf(first_active-1,1) = 0.0
    !Surf(last_active+1,1)  = 0.0
    do i=1,ncells
    position(i,1)=radii_c(i)
    end do


end subroutine gridinit_sphere1D

#endif
#endif





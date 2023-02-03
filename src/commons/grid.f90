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


!Get the x index of a cell according to i
function iyy(i)
  use parameters
  use commons
  use units
  implicit none
  integer  :: i
  integer  :: iyy

  iyy (i-1)/nx_max+1
#endif

end function iyy

subroutine get_active_cells
  use parameters
  use commons
  use units
  implicit none
  integer  :: ix,iy,icell,ii,icount,ixx,iyy
  allocate(active_cell(1:ncells))
  active_cell=0
  if(ndim==1) then
  	do ix=1,nx_max
  		if(ix.ge.first_active.and.ix.le.last_active) active_cell(ix)=1
  	end do
  else
    icount=0
  	do ix=1,nx_max
  		do iy=1,ny_max
  		ii=icell(ix,iy)
      print *, ii, ix, iy,ixx(ii),iyy(ii)
  		if(ix.ge.first_active.and.ix.le.last_active) then
  			if(iy.ge.first_active_y.and.iy.le.last_active_y) then
  		 	  active_cell(ii)=1
          icount=icount+1
  		 	endif
  		  endif
  		end do
  	end do	
  endif
  print *, icount, ' active cells were prepared.'
end subroutine get_active_cells


!Grid initialisation routine
subroutine gridinit(rmax_x,rmax_y)
  use parameters
  use commons
  use units
  implicit none

  real(dp):: rmax_x,rmax_y
  integer :: i, ix, iy,icell

  print *, 'Number of cells        =', Ncells
  print *, 'Number of active cells =', Ncells_active
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

          vol (icell(ix,iy))   = dx(icell(ix,iy),1)*dx(icell(ix,iy),2)
          position(icell(ix,iy),1)= (float(ix-first_active)+0.5d0) * rmax_x/float(nx)
          position(icell(ix,iy),2)= (float(iy-first_active_y)+0.5d0) * rmax_y/float(ny)
      end do
    end do
  endif
#endif

#if GEOM==1 
! Polar geometry. Only in 2D
      do ix = 1, nx_max
        do iy = 1, ny_max
          dx  (icell(ix,iy),1) = rmax_x/float(nx)
          dx  (icell(ix,iy),2) = rmax_y/float(ny) ! This becomes theta

          polar_radii(icell(ix,iy)) = (float(ix-first_active)+0.5d0)   * rmax_x/float(nx) ! R 
          theta(icell(ix,iy))       = (float(iy-first_active_y)+0.5d0) * rmax_y/float(ny! Theta

          surf(icell(ix,iy),1) = dx(icell(ix,iy),1)
          surf(icell(ix,iy),2) = dx(icell(ix,iy),2)

          vol (icell(ix,iy))   = dx(icell(ix,iy),1)*dx(icell(ix,iy),2)


          position(icell(ix,iy),1) =  polar_radii(icell(ix,iy))*cos(theta(icell(ix,iy)))! x
          position(icell(ix,iy),2) =  polar_radii(icell(ix,iy))*sin(theta(icell(ix,iy)))! y
      end do
    end do
#endif
end subroutine gridinit
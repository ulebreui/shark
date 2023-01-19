!Get the index of a cell according to ix and iy
function icell(ix,iy)
  use parameters
  use commons
  use units
  implicit none
  integer  :: ix,iy
  integer  :: icell
  if(ndim==1) icell = ix
  if(ndim==2) icell = 1+(ix-1)+(iy-1)*nx_max
end function icell

subroutine get_active_cells
  use parameters
  use commons
  use units
  implicit none
  integer  :: ix,iy,icell,ii,icount
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
subroutine gridinit(rmax)
  use parameters
  use commons
  use units
  implicit none

  real(dp):: rmax
  integer :: i, ix, iy,icell

  print *,'You are using a linear cartesian grid.'
  print *, 'Number of cells        =', Ncells
  print *, 'Number of active cells =', Ncells_active

  if(ndim==1) then
    do ix = 1, ncells
      dx (ix,1)  = rmax/float(nx)
      vol(ix)    = dx(ix,1)
      surf(ix,1) = 1
      if(active_cell(ix)==1) position(ix,1)=position(ix-1,1)+dx (ix,1)
    end do
  else
      do ix = 1, nx_max
        do iy = 1, ny_max
          dx  (icell(ix,iy),1) = rmax/float(nx)
          dx  (icell(ix,iy),2) = rmax/float(ny)

          surf(icell(ix,iy),1) = dx(icell(ix,iy),1)
          surf(icell(ix,iy),2) = dx(icell(ix,iy),2)

          vol (icell(ix,iy))   = dx(icell(ix,iy),1)*dx(icell(ix,iy),2)
          position(icell(ix,iy),1)= float(ix-first_active+1) * rmax/float(nx)
          position(icell(ix,iy),2)= float(iy-first_active_y+1) * rmax/float(ny)
      end do
    end do
  endif
end subroutine gridinit
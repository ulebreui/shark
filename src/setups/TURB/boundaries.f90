!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This applies the boundaries either to uold or unew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine apply_boundaries(who_app,uu,nn,nn2)
  use parameters
  use commons
  use units
  implicit none
  integer :: who_app,idust,ighost,nn,nn2,ix,iy,icell,ii,icount
  real(dp), dimension (1:nn,1:nn2) :: uu

  if(ndim==1) then
   do ix = 1,nghost
      uu(ix,:)          = uu(first_active,:)
      uu(nx_max+1-ix,:) = uu(last_active,:)
   end do
  else
  !icount=0
  do ix = first_active,last_active 
   do iy = 1,nghost
      !icount=icount+1
      !uu(icell(ix,iy),:)          = uu(icell(ix,first_active_y),:) ! Zero gradient along y
      !uu(icell(ix,ny_max+1-iy),:) = uu(icell(ix,last_active_y),:)
      !print *, iy, last_active_y-nghost+iy, ny_max ,'first bound'
      !print *, ny_max+1-iy, first_active_y+nghost-iy, ny_max, '2nd bound'

      uu(icell(ix,iy),:)          = uu(icell(ix,last_active_y-nghost+iy),:) ! Periodic along y
      uu(icell(ix,ny_max+1-iy),:) = uu(icell(ix,first_active_y+nghost-iy),:)

   end do
   !stop
  end do
   do ix=1,nghost ! Periodic along x
     do iy=1,ny_max
         !icount=icount+1
         !print *, ix,ny_max+1-ix, first_active,last_active_y         
         uu(icell(ix,iy),:)          = uu(icell(last_active-nghost+ix,iy),:)
         uu(icell(nx_max+1-ix,iy),:) = uu(icell(first_active+nghost-ix,iy),:)
      end do
      !stop
  end do
  endif
  !print *, icount, ' ghost cells have been counted'
  if(who_app .eq. 1) then
     uold=uu
  else if(who_app .eq. 2) then
     unew=uu
  endif
end subroutine apply_boundaries



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This applies the boundaries either to u_prim or unew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine apply_boundaries(who_app,u_prim,nn,nn2)
  use parameters
  use commons
  use units
  implicit none
  integer :: who_app,idust,ighost,nn,nn2,ix,iy,icell,ii,icount
  real(dp), dimension (1:nn,1:nn2) :: u_prim

  if(ndim==1) then
   do ix = 1,nghost
      u_prim(ix,:)          = u_prim(first_active,:)
      u_prim(nx_max+1-ix,:) = u_prim(last_active,:)
   end do
  else
  !icount=0
  do ix = first_active,last_active 
   do iy = 1,nghost
      !icount=icount+1
      !u_prim(icell(ix,iy),:)          = u_prim(icell(ix,first_active_y),:) ! Zero gradient along y
      !u_prim(icell(ix,ny_max+1-iy),:) = u_prim(icell(ix,last_active_y),:)
      !print *, iy, last_active_y-nghost+iy, ny_max ,'first bound'
      !print *, ny_max+1-iy, first_active_y+nghost-iy, ny_max, '2nd bound'

      u_prim(icell(ix,iy),:)          = u_prim(icell(ix,last_active_y-nghost+iy),:) ! Periodic along y
      u_prim(icell(ix,ny_max+1-iy),:) = u_prim(icell(ix,first_active_y+nghost-iy),:)

   end do
   !stop
  end do
   do ix=1,nghost ! Periodic along x
     do iy=1,ny_max
         !icount=icount+1
         !print *, ix,ny_max+1-ix, first_active,last_active_y         
         u_prim(icell(ix,iy),:)          = u_prim(icell(last_active-nghost+ix,iy),:)
         u_prim(icell(nx_max+1-ix,iy),:) = u_prim(icell(first_active+nghost-ix,iy),:)
      end do
      !stop
  end do
  endif
  !print *, icount, ' ghost cells have been counted'

end subroutine apply_boundaries



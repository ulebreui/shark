!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This applies the boundaries either to u_prim or unew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine apply_boundaries
  use parameters
  use commons
  use units
  implicit none
  integer :: who_app,idust,ighost,nn,nn2,ix,iy,icell,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar
  real(dp) :: total_energy,ekin_old,ekin_new,pm

  if(ndim==1) then
   do ix = 1,nghost
      u_prim(ix,:)          = u_prim(first_active,:)
      u_prim(nx_max+1-ix,:) = u_prim(last_active,:)
   end do
  else
  !icount=0
  do ix = first_active,last_active  
   do iy = 1,nghost
      ibound_left   = icell(ix,iy)
      ibound_right  = icell(ix,ny_max+1-iy)
      !i_active_left = icell(ix,last_active_y-nghost+iy)
      !i_active_right= icell(ix,first_active_y+nghost-iy)

      i_active_left = first_active_y
      i_active_right= last_active_y      
      do ivar =1,nvar
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! Periodic along y
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
      end do
   end do
   !stop
  end do
   do ix=1,nghost 
     do iy=1,ny_max
         ibound_left    = icell(ix,iy)
         ibound_right   = icell(nx_max+1-ix,iy)
         !i_active_left  = icell(last_active-nghost+ix,iy)
         !i_active_right = icell(first_active+nghost-ix,iy)
         i_active_left = first_active!
         i_active_right= last_active!    
         do ivar =1,nvar
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! Periodic conditions along x
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
         end do     

      end do
      !stop
  end do
  endif

end subroutine apply_boundaries





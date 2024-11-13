module boundary_types
contains

subroutine disk_boundaries
  use parameters
  use commons
  use units
  implicit none
  integer :: idust,ighost,ix,iy,icell,ivar


  do ix = first_active,last_active  
   do iy = 1,nghost
      do ivar = 1,nvar
         u_prim(ix,iy,ivar)  = u_prim(ix,last_active_y -nghost+iy,ivar) 
         u_prim(ix,ny_max+1-iy,ivar) = u_prim(ix,first_active_y+nghost-iy,ivar)
      end do

   end do

  end do
   do ix= 1,nghost 
     do iy= 1,ny_max  
         do ivar = 1,nvar
            u_prim(ix,iy,ivar)  = u_prim(first_active,iy,ivar) 
            u_prim(nx_max+1-ix,iy,ivar) = u_prim(last_active,iy,ivar)
         end do     
         u_prim(nx_max+1-ix,iy,ivx)  = min(u_prim(nx_max+1-ix,iy,ivx),0.0d0)
         u_prim(ix,iy,ivx)   = max(u_prim(ix,iy,ivx) ,0.0d0)
 
      end do
  end do

end subroutine disk_boundaries


subroutine periodic_boundaries
  use parameters
  use commons
  use units
  implicit none
  integer  :: idust,ighost,ix,iy,icell, ibound_left, ibound_right, i_active_left, i_active_right,ivar

  do ix = first_active,last_active  
   do iy = 1,nghost
      do ivar =1,nvar
         u_prim(ix,iy,ivar)  = u_prim(ix,last_active_y-nghost+iy,ivar) ! Periodic along y
         u_prim(ix,ny_max+1-iy,ivar) = u_prim(ix,first_active_y+nghost-iy,ivar)
      end do
   end do
  end do
   do ix= 1,nghost 
     do iy= 1,ny_max   
         do ivar = 1,nvar
            u_prim(ix,iy,ivar)  = u_prim(last_active-nghost+ix,iy,ivar) ! Periodic conditions along x
            u_prim(nx_max+1-ix,iy,ivar) = u_prim(first_active+nghost-ix,iy,ivar)
         end do     

      end do
      !stop
  end do

end subroutine periodic_boundaries


subroutine shear_boundaries(vel_shear)
  use parameters
  use commons
  use units
  implicit none
  integer  :: idust,ighost,ix,iy,icell,ivar
  real(dp) :: pm,mom_new,vel_shear

   do ix  = first_active,last_active
   do iy  = 1, nghost
      do ivar =1,nvar
         u_prim(ix,iy,ivar)          = u_prim(ix,last_active_y-nghost+iy,ivar) ! Periodic along y
         u_prim(ix,ny_max+1-iy,ivar) = u_prim(ix,first_active_y+nghost-iy,ivar)
      end do
   end do
  end do

   do ix=1,nghost ! Shear along x
     do iy= 1,ny_max
         do ivar =1,nvar
            u_prim(ix,iy,ivar)          = u_prim(last_active-nghost+ix,iy,ivar) ! We first apply the periodic conditions along x
            u_prim(nx_max+1-ix,iy,ivar) = u_prim(first_active+nghost-ix,iy,ivar)
         end do     
     
         ! Then, we compensate the shear in the z velocity

         u_prim(ix,iy,ivz)           = u_prim(last_active-nghost+ix,iy,ivz) - u_prim(last_active-nghost+ix,iy,irho)*vel_shear
         u_prim(nx_max+1-ix,iy,ivz)  = u_prim(first_active+nghost-ix,iy,ivz) + u_prim(first_active+nghost-ix,iy,irho)*vel_shear     
#if NDUST>0         
         do idust=1,ndust
            u_prim(ix,iy ,ivdz(idust))             = u_prim(last_active-nghost+ix,iy,ivdz(idust)) - u_prim(last_active-nghost+ix,iy,irhod(idust))*vel_shear
            u_prim(nx_max+1-ix,iy,ivdz(idust))     = u_prim(first_active+nghost-ix,iy ,ivdz(idust)) + u_prim(first_active+nghost-ix,iy,irhod(idust)) *vel_shear
         end do
#endif  
      end do
      !stop
  end do



end subroutine shear_boundaries

end module boundary_types

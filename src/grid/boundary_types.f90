module boundary_types
contains

subroutine boundary_collapse_1D
  use parameters
  use commons
  use units
  implicit none
  real(dp) :: cs_eos,barotrop
  integer :: who_app,idust,ighost,nn,nn2,ix,iy,icell,ii,icount,ibound_left, ibound_right, i_active_left, i_active_right,ivar

   do ix=1,nghost 
     do iy=1,ny_max
         ibound_left    = icell(ix,iy)
         ibound_right   = icell(nx_max+1-ix,iy)

         i_active_left  = icell(first_active,iy)
         i_active_right = icell(last_active,iy)

        do ivar =1,nvar
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! Zero gradient
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
        end do 

        u_prim(ibound_left,ivx)   = 0.0d0
        u_prim(first_active,ivx)  = 0.0d0
        u_prim(ibound_right,ivx)  = min(u_prim(ibound_right,ivx),0.0d0)
        u_prim(last_active,ivx)   = min(u_prim(last_active,ivx),0.0d0)
#if NDUST>0
        do idust=1,ndust
            u_prim(ibound_left ,ivdx(idust))   = 0.0d0
            u_prim(first_active,ivdx(idust))   = 0.0d0
            u_prim(ibound_right,ivdx(idust))   = min(u_prim(ibound_right,ivdx(idust)),0.0d0)
            u_prim(last_active ,ivdx(idust))   = min(u_prim(last_active ,ivdx(idust)),0.0d0)
        end do
#endif

      end do
  
  end do


end subroutine boundary_collapse_1D



subroutine zero_gradients
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
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! Zero grad along y
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
      end do

       !  u_prim(ibound_left,ivy)  = max(u_prim(ibound_left,ivy),0.0d0)
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
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! Zero grad along x
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
         end do     

      end do
      !stop
  end do
  endif

end subroutine zero_gradients

subroutine disk_boundaries
  use parameters
  use commons
  use units
  implicit none
  integer :: who_app,idust,ighost,nn,nn2,ix,iy,icell,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar


  !icount=0
  do ix = first_active,last_active  
   do iy = 1,nghost
      ibound_left    = icell(ix,iy)
      ibound_right   = icell(ix,ny_max+1-iy)
      i_active_left  = icell(ix,last_active_y -nghost+iy)
      i_active_right = icell(ix,first_active_y+nghost-iy)
   
      do ivar = 1,nvar
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) 
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
      end do

   end do

   !stop
  end do
   do ix= 1,nghost 
     do iy= 1,ny_max
         ibound_left    = icell(ix,iy)
         ibound_right   = icell(nx_max+1-ix,iy)
         !i_active_left  = icell(last_active-nghost+ix,iy)
         !i_active_right = icell(first_active+nghost-ix,iy)
         i_active_left  = first_active!
         i_active_right = last_active!    
         do ivar = 1,nvar
            u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) 
            u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
         end do     
         u_prim(ibound_right,ivx)  = min(u_prim(ibound_right,ivx),0.0d0)
         u_prim(ibound_left,ivx)   = max(u_prim(ibound_left,ivx) ,0.0d0)
 
      end do
      !stop
  end do

end subroutine disk_boundaries


subroutine periodic_boundaries
  use parameters
  use commons
  use units
  implicit none
  integer :: who_app,idust,ighost,nn,nn2,ix,iy,icell,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar
  real(dp) :: total_energy,ekin_old,ekin_new,pm

  do ix = first_active,last_active  
   do iy = 1,nghost
      ibound_left    = icell(ix,iy)
      ibound_right   = icell(ix,ny_max+1-iy)
      i_active_left  = icell(ix,last_active_y-nghost+iy)
      i_active_right = icell(ix,first_active_y+nghost-iy)
      do ivar =1,nvar
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! Periodic along y
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
      end do
   end do
   !stop
  end do
   do ix= 1,nghost 
     do iy= 1,ny_max
         ibound_left    = icell(ix,iy)
         ibound_right   = icell(nx_max+1-ix,iy)
         i_active_left  = icell(last_active-nghost+ix,iy)
         i_active_right = icell(first_active+nghost-ix,iy)
   
         do ivar = 1,nvar
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! Periodic conditions along x
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
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
  integer  :: idust,ighost,nn,nn2,ix,iy,icell,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar
  real(dp) :: total_energy,ekin_old,ekin_new,pm,mom_new,vel_shear

!vel_shear=q_shear*Omega_shear*box_l
   do ix  = first_active,last_active
   do iy  = 1, nghost
        !print *,'is it a boundary ? ', active_cell(icell(ix,iy))

      ibound_left    = icell(ix,iy)
      ibound_right   = icell(ix,ny_max+1-iy)
      i_active_left  = icell(ix,last_active_y-nghost+iy)
      i_active_right = icell(ix,first_active_y+nghost-iy)
      do ivar =1,nvar
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! Periodic along y
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
      end do
      ! We must apply the shear to the corner cells 
      if(ix==first_active .or. ix==last_active) then
         if(ix==first_active) pm = -1.0d0
         if(ix==last_active)  pm = 1.0d0

         mom_new  = u_prim(i_active_left,ivz)  + pm * u_prim(i_active_left,irho)*vel_shear
         mom_new  = u_prim(i_active_right,ivz) + pm * u_prim(i_active_right,irho)*vel_shear

      end if
   end do
  end do

   do ix=1,nghost ! Shear along x
     do iy= 1,ny_max
         ibound_left    = icell(ix,iy)
         ibound_right   = icell(nx_max+1-ix,iy)
         i_active_left  = icell(last_active-nghost+ix,iy)
         i_active_right = icell(first_active+nghost-ix,iy)

         do ivar =1,nvar
            u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! We first apply the periodic conditions along x
            u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
         end do     
         
         ! Then, we compensate the shear in the z velocity

         mom_new = u_prim(i_active_left,ivz) - u_prim(i_active_left,irho)*vel_shear
         u_prim(ibound_left,ivz)        = mom_new
         mom_new  =  u_prim(i_active_right,ivz) + u_prim(i_active_right,irho)*vel_shear
         u_prim(ibound_right,ivz)        = mom_new         
#if NDUST>0         
         do idust=1,ndust
            u_prim(ibound_left ,ivdz(idust))     = u_prim(i_active_left,ivdz(idust)) - u_prim(i_active_left,irhod(idust))*vel_shear
            u_prim(ibound_right,ivdz(idust))     = u_prim(i_active_right ,ivdz(idust)) + u_prim(i_active_right,irhod(idust)) *vel_shear
         end do
#endif  
      end do
      !stop
  end do



end subroutine shear_boundaries

end module boundary_types

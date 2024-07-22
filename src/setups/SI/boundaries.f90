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
  integer  :: idust,ighost,nn,nn2,ix,iy,icell,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar
  real(dp) :: total_energy,ekin_old,ekin_new,pm,mom_new


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

         mom_new= u_prim(i_active_left,ivz) + pm * u_prim(i_active_left,irho)*q_shear*Omega_shear*box_l
         !u_prim(ibound_left,ivz)        = mom_new
         !u_prim(ibound_left,iP) =  u_prim(i_active_left,iP) +  half*(mom_new**2-u_prim(i_active_left,ivz)**2.)/u_prim(i_active_left,irho)

         mom_new  = u_prim(i_active_right,ivz) + pm * u_prim(i_active_right,irho)*q_shear*Omega_shear*box_l
         !u_prim(ibound_right,ivz)        = mom_new
         !u_prim(ibound_right,iP) = u_prim(i_active_right,iP) +  half*(mom_new**2-u_prim(i_active_right,ivz)**2.)/u_prim(i_active_right,irho)
#if NDUST>0         
         do idust=1,ndust
           !u_prim(ibound_left,ivdz(idust))     = u_prim(i_active_left,ivdz(idust))   + pm*u_prim(i_active_left,irhod(idust))*q_shear*Omega_shear*box_l
           !u_prim(ibound_right,ivdz(idust))    = u_prim(i_active_right,ivdz(idust))  + pm*u_prim(i_active_right,irhod(idust))*q_shear*Omega_shear*box_l
         end do
#endif  
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

         mom_new = u_prim(i_active_left,ivz) - u_prim(i_active_left,irho)*q_shear*Omega_shear*box_l
         u_prim(ibound_left,ivz)        = mom_new
         !u_prim(ibound_left,iP) = u_prim(i_active_left,iP) +  half*(mom_new**2-u_prim(i_active_left,ivz)**2.)/u_prim(i_active_left,irho)

         mom_new  =  u_prim(i_active_right,ivz) + u_prim(i_active_right,irho)*q_shear*Omega_shear*box_l
         u_prim(ibound_right,ivz)        = mom_new         
         !u_prim(ibound_right,iP) = u_prim(i_active_right,iP) +  half*(mom_new**2-u_prim(i_active_right,ivz)**2.)/u_prim(i_active_right,irho)
#if NDUST>0         
         do idust=1,ndust
            u_prim(ibound_left ,ivdz(idust))     = u_prim(i_active_left,ivdz(idust)) - u_prim(i_active_left,irhod(idust))*q_shear*Omega_shear*box_l
            u_prim(ibound_right,ivdz(idust))     = u_prim(i_active_right ,ivdz(idust)) + u_prim(i_active_right,irhod(idust)) *q_shear*Omega_shear*box_l
         end do
#endif  
      end do
      !stop
  end do



end subroutine apply_boundaries



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
  use boundary_types
  use parameters
  use commons
  use units
  implicit none
  integer  :: idust,ighost,nn,nn2,ix,iy,icell,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar
  real(dp) :: total_energy,ekin_old,ekin_new,pm,mom_new

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

         mom_new = u_prim(last_active-nghost+ix,iy,ivz) - u_prim(last_active-nghost+ix,iy,irho)*q_shear*Omega_shear*box_l
         u_prim(ix,iy,ivz)        = mom_new

         mom_new  =  u_prim(first_active+nghost-ix,iy,ivz) + u_prim(first_active+nghost-ix,iy,irho)*q_shear*Omega_shear*box_l
         u_prim(nx_max+1-ix,iy,ivz)        = mom_new         
#if NDUST>0         
         do idust=1,ndust
            u_prim(ix,iy ,ivdz(idust))     = u_prim(last_active-nghost+ix,iy,ivdz(idust)) - u_prim(last_active-nghost+ix,iy,irhod(idust))*q_shear*Omega_shear*box_l
            u_prim(nx_max+1-ix,iy,ivdz(idust))     = u_prim(first_active+nghost-ix,iy ,ivdz(idust)) + u_prim(first_active+nghost-ix,iy,irhod(idust)) *q_shear*Omega_shear*box_l
         end do
#endif  
      end do
      !stop
  end do



end subroutine apply_boundaries



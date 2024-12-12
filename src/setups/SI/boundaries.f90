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
  use OMP_LIB
  implicit none
  integer  :: idust,ighost,nn,nn2,ix,iy,icell,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar
  real(dp) :: total_energy,ekin_old,ekin_new,pm,mom_new

   !$omp parallel do default(shared) schedule(RUNTIME) private(ivar, ix,iy)
      do ix  = first_active,last_active
         do iy  = 1, nghost
            do ivar =1,nvar
               u_prim(ivar,ix,iy)          = u_prim(ivar,ix,last_active_y-nghost+iy) ! Periodic along y
               u_prim(ivar,ix,ny_max+1-iy) = u_prim(ivar,ix,first_active_y+nghost-iy)
            end do
         end do
  end do
   !$omp parallel do default(shared) schedule(RUNTIME) private(ivar, ix, iy, idust)
   do iy= 1,ny_max
      do ix=1,nghost ! Shear along x
         do ivar =1,nvar
            u_prim(ivar,ix,iy)          = u_prim(ivar,last_active-nghost+ix,iy) ! We first apply the periodic conditions along x
            u_prim(ivar,nx_max+1-ix,iy) = u_prim(ivar,first_active+nghost-ix,iy)
         end do     
         ! Then, we compensate the shear in the z velocity

         u_prim(ivz,ix,iy) =  u_prim(ivz,last_active-nghost+ix,iy) - u_prim(irho,last_active-nghost+ix,iy)*q_shear*Omega_shear*box_l
         u_prim(ivz,nx_max+1-ix,iy)  = u_prim(ivz,first_active+nghost-ix,iy) + u_prim(irho,first_active+nghost-ix,iy)*q_shear*Omega_shear*box_l        
#if NDUST>0         
         do idust=1,ndust
            u_prim(ivdz(idust),ix,iy)     = u_prim(ivdz(idust),last_active-nghost+ix,iy) - u_prim(irhod(idust),last_active-nghost+ix,iy)*q_shear*Omega_shear*box_l
            u_prim(ivdz(idust),nx_max+1-ix,iy)     = u_prim(ivdz(idust),first_active+nghost-ix,iy) + u_prim(irhod(idust),first_active+nghost-ix,iy) *q_shear*Omega_shear*box_l
         end do
#endif  
      end do
      !stop
  end do



end subroutine apply_boundaries



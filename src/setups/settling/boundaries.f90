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
  integer  :: idust,ighost,nn,nn2,ix,iy,icell,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar,inner,outer
  real(dp) :: total_energy,ekin_old,ekin_new,pm,mom_new

   !$omp parallel do default(shared) schedule(RUNTIME) private(ivar, ix,iy)
      do ix  = first_active,last_active
         do iy  = 1, nghost
            do ivar =1,nvar
               u_prim(ivar,ix,iy)          = u_prim(ivar,ix,first_active_y)
               u_prim(ivar,ix,ny_max+1-iy) = u_prim(ivar,ix,last_active_y)
            end do

            ! u_prim(ivy,ix,iy)          = max(u_prim(ivy,ix,first_active_y),0.d0)
            ! u_prim(ivy,ix,ny_max+1-iy) = min(u_prim(ivy,ix,last_active_y),0.0d0)
            ! do idust=1,ndust
            !    u_prim(ivdy(idust),ix,iy)          = max(u_prim(ivdy(idust),ix,first_active_y),0.d0)
            !    u_prim(ivdy(idust),ix,ny_max+1-iy) = min(u_prim(ivdy(idust),ix,last_active_y),0.0d0)
            ! end do

         end do
  end do
   !$omp parallel do default(shared) schedule(RUNTIME) private(ivar, ix, iy, idust,inner,outer)
   do iy= 1,ny_max
      do ix=1,nghost 
         do ivar =1,nvar
            u_prim(ivar,ix,iy)          = u_prim(ivar,last_active-nghost+ix,iy) ! We first apply the periodic conditions along x
            u_prim(ivar,nx_max+1-ix,iy) = u_prim(ivar,first_active+nghost-ix,iy)
         end do     
      end do
      !stop
  end do



end subroutine apply_boundaries



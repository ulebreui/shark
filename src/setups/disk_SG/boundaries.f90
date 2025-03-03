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
    integer :: i,idust,icell,ipscal
    real(dp):: unit_lout,sigma_out,t_rel,vrout,maccreted,rhod_old

    integer :: who_app,ighost,nn,nn2,ix,iy,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar


  !icount=0
  do ix = first_active,last_active  
   do iy = 1,nghost   
      do ivar = 1,nvar
         u_prim(ivar,ix,iy)          = u_prim(ivar,ix,last_active_y  - nghost+iy) 
         u_prim(ivar,ix,ny_max+1-iy) = u_prim(ivar,ix,first_active_y + nghost-iy)
      end do

   end do

   !stop
  end do
   do ix= 1,nghost 
     do iy= 1,ny_max
        !u_prim(ivx,last_active,iy)  = 0.0d0
        u_prim(ivx,first_active,iy) = 0.0d0

         do ivar = 1,nvar
            u_prim(ivar,ix,iy)          = u_prim(ivar,first_active,iy) 
            u_prim(ivar,nx_max+1-ix,iy) = u_prim(ivar,last_active,iy)
         end do 
        ! ! Prevents any inflow
        ! u_prim(ivx,nx_max+1-ix,iy)  = min(u_prim(ivx,nx_max+1-ix,iy),0.0d0)
        !u_prim(ivx,ix,iy)           = max(u_prim(ivx,ix,iy) ,0.0d0)  
        ! u_prim(ivx,first_active,iy) = min(u_prim(ivx,first_active,iy),0.0d0)
        !u_prim(ivx,last_active,iy)  = max(u_prim(ivx,last_active,iy) ,0.0d0)  
#if NDUST>0
         print *, 'Careful dust boundary not correctly implemented'
         stop  
        do idust=1,ndust
           u_prim(ivdx(idust),nx_max+1-ix,iy)  = min(u_prim(ivdx(idust),nx_max+1-ix,iy),0.0d0)
           u_prim(ivdx(idust),ix,iy)           = max(u_prim(ivdx(idust),ix,iy) ,0.0d0)      
        end do
#endif
      end do
      !stop
  end do

end subroutine apply_boundaries

subroutine apply_boundaries_phi
    use boundary_types
    use parameters
    use commons
    use units
    implicit none
    integer :: i,idust,icell,ipscal
    real(dp):: unit_lout,sigma_out,t_rel,vrout,maccreted,rhod_old

    integer :: who_app,ighost,nn,nn2,ix,iy,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar


  !icount=0
  do ix = first_active,last_active  
   do iy = 1,nghost   

          phi_grav(ix,iy)          =  phi_grav(ix,last_active_y  - nghost+iy) 
          phi_grav(ix,ny_max+1-iy) =  phi_grav(ix,first_active_y + nghost-iy)

   end do

   !stop
  end do
   do ix= 1,nghost 
     do iy= 1,ny_max
        phi_grav(ix,iy)          =  phi_grav(first_active,iy) 
        phi_grav(nx_max+1-ix,iy) =  phi_grav(last_active,iy)


      end do
      !stop
  end do

end subroutine apply_boundaries_phi


subroutine apply_boundaries_force
    use boundary_types
    use parameters
    use commons
    use units
    implicit none
    integer :: i,idust,icell,ipscal
    real(dp):: unit_lout,sigma_out,t_rel,vrout,maccreted,rhod_old

    integer :: who_app,ighost,nn,nn2,ix,iy,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar


  !icount=0
  do ix = first_active,last_active  
   do iy = 1,nghost   

          force_x(ix,iy)          =  force_x(ix,last_active_y  - nghost+iy) 
          force_x(ix,ny_max+1-iy) =  force_x(ix,first_active_y + nghost-iy)
          force_y(ix,iy)          =  force_y(ix,last_active_y  - nghost+iy) 
          force_y(ix,ny_max+1-iy) =  force_y(ix,first_active_y + nghost-iy)
   end do

   !stop
  end do
   do ix= 1,nghost 
     do iy= 1,ny_max
        force_x(ix,iy)          =  force_x(first_active,iy) 
        force_x(nx_max+1-ix,iy) =  force_x(last_active,iy)
        force_y(ix,iy)          =  force_y(first_active,iy) 
        force_y(nx_max+1-ix,iy) =  force_y(last_active,iy)

      end do
      !stop
  end do

end subroutine apply_boundaries_force



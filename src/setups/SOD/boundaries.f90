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
  integer :: who_app,idust,ighost,nn,nn2,ix,iy,icell,ii,icount,ibound_left, ibound_right, i_active_left, i_active_right,ivar

#if NY>1
   do ix = first_active,last_active  
   do iy = 1,nghost
      ibound_left   = icell(ix,iy)
      ibound_right  = icell(ix,ny_max+1-iy)


      i_active_left = icell(ix,last_active_y-nghost+iy)
      i_active_right= icell(ix,first_active_y+nghost-iy)

     if(direction_shock==1) then 
         u_prim(ibound_left,irho)   = rho_l
         u_prim(ibound_right,irho)  = rho_r
         u_prim(ibound_left,iP)     = P_l/(gamma-1.0d0)
         u_prim(ibound_right,iP)    = P_r/(gamma-1.0d0)
         u_prim(ibound_left,ivx)     = 0.0d0
         u_prim(ibound_right,ivx)    = 0.0d0       
         u_prim(ibound_left,ivy)    = 0.0d0
         u_prim(ibound_right,ivy)   = 0.d0
         u_prim(ibound_left,ivz)    = 0.0d0
         u_prim(ibound_right,ivz)   = 0.d0
    else
        do ivar =1,nvar
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! Periodic along y
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
      end do 
    endif 

   end do
  end do

#endif
   do ix=1,nghost 
     do iy=1,ny_max
         ibound_left    = icell(ix,iy)
         ibound_right   = icell(nx_max+1-ix,iy)

        i_active_left  = icell(last_active-nghost+ix,iy)
        i_active_right = icell(first_active+nghost-ix,iy)
        if(direction_shock==0) then 
         u_prim(ibound_left,irho)   = rho_l
         u_prim(ibound_right,irho)  = rho_r
         u_prim(ibound_left,iP)     = P_l/(gamma-1.0d0)
         u_prim(ibound_right,iP)    = P_r/(gamma-1.0d0)
         u_prim(ibound_left,ivx)    = 0.0d0
         u_prim(ibound_right,ivx)   = 0.0d0       
         u_prim(ibound_left,ivy)    = 0.0d0
         u_prim(ibound_right,ivy)   = 0.d0
         u_prim(ibound_left,ivz)    = 0.0d0
         u_prim(ibound_right,ivz)   = 0.d0
    else 
        do ivar =1,nvar
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) ! Periodic along x
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
      end do 
    endif 
  

      end do
  
  end do


end subroutine apply_boundaries



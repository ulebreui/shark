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
      ibound_left    = icell(ix,iy)
      ibound_right   = icell(ix,ny_max+1-iy)
      !i_active_left  = icell(ix,last_active_y -nghost+iy)
      !i_active_right = icell(ix,first_active_y+nghost-iy)
      i_active_left  = first_active_y
      i_active_right = last_active_y   
      do ivar = 1,nvar
         u_prim(ibound_left,ivar)  = u_prim(i_active_left,ivar) 
         u_prim(ibound_right,ivar) = u_prim(i_active_right,ivar)
      end do
!       u_prim(ibound_right,ivy)  = min(u_prim(ibound_right,ivy),0.0d0)
!       u_prim(ibound_left,ivy)   = max(u_prim(ibound_left,ivy) ,0.0d0)

! #if NDUST>0
!         do idust=1,ndust
!            u_prim(ibound_right,ivdy(idust))  = min(u_prim(ibound_right,ivdy(idust)),0.0d0)
!            u_prim(ibound_left,ivdy(idust))   = max(u_prim(ibound_left,ivdy(idust)) ,0.0d0)      
!         end do
! #endif
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
!          u_prim(ibound_right,ivx)  = min(u_prim(ibound_right,ivx),0.0d0)
!          u_prim(ibound_left,ivx)   = max(u_prim(ibound_left,ivx) ,0.0d0)
! #if NDUST>0
!         do idust=1,ndust
!            u_prim(ibound_right,ivdx(idust))  = min(u_prim(ibound_right,ivdx(idust)),0.0d0)
!            u_prim(ibound_left,ivdx(idust))   = max(u_prim(ibound_left,ivdx(idust)) ,0.0d0)      
!         end do
! #endif
      end do
      !stop
  end do

!   !Accretion in the outer boundary
!    if(accretion) then
!    do ix = 1,nghost+1
!     do iy=first_active_y,last_active_y
!     ibound_right   = icell(nx_max+1-ix,iy)
!         if(phi(icell(last_active,iy))<=phi_mom*pi) then
!             maccreted = M_acc/(365.25*24.*3600.)*unit_t
!             vrout     = cs(icell(last_active,iy))
!             sigma_out = maccreted/vrout/phi_mom/pi/box_l

!             u_prim(ibound_right,irho) =  sigma_out
!             u_prim(ibound_right,ivx)  =- sigma_out*vrout
!             u_prim(ibound_right,ivy)  =  sigma_out*L_out*sqrt(Mstar/100.)*100./box_l
! #if NDUST>0
!            do idust=1,ndust
! ! #if NDUSTPSCAL>0
! !             rhod_old=u_prim(ibound_right,irhod(idust))
! ! #endif           
!              u_prim(ibound_right,irhod(idust))  =  dust2gas*sigma_out
!              u_prim(ibound_right,ivdx(idust))   =  -dust2gas*sigma_out*vrout
!              u_prim(ibound_right,ivdy(idust))   =  dust2gas*sigma_out*L_out*sqrt(Mstar/100.)*100./box_l
! #if NDUSTPSCAL>0
!             !do ipscal=1,ndustpscal
!             u_prim(ibound_right,idust_pscal(idust,1))  = dust2gas*sigma_out*scut/unit_l
!             !end do
! #endif     
!            end do
! #endif
!         endif
!    end do
!    end do
!    endif
end subroutine apply_boundaries


subroutine apply_boundaries_phi
    use boundary_types
    use parameters
    use commons
    use units
    implicit none
    integer :: i,idust,icell
    real(dp):: unit_lout,sigma_out,t_rel,vrout,maccreted

    integer :: who_app,ighost,nn,nn2,ix,iy,ii,icount, ibound_left, ibound_right, i_active_left, i_active_right,ivar


  !icount=0
  do ix = first_active,last_active  
   do iy = 1,nghost
      ibound_left    = icell(ix,iy)
      ibound_right   = icell(ix,ny_max+1-iy)
      i_active_left  = icell(ix,last_active_y -nghost+iy)
      i_active_right = icell(ix,first_active_y+nghost-iy)
   
      phi_sg(ibound_left)       = phi_sg(i_active_left) 
      phi_sg(ibound_right)      = phi_sg(i_active_right)

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
         phi_sg(ibound_left)  = phi_sg(i_active_left) 
         phi_sg(ibound_right) = phi_sg(i_active_right)
      end do
      !stop
  end do


end subroutine apply_boundaries_phi




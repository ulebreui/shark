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
         do ivar = 1,nvar
            u_prim(ivar,ix,iy)          = u_prim(ivar,first_active,iy) 
            u_prim(ivar,nx_max+1-ix,iy) = u_prim(ivar,last_active,iy)
         end do     
         u_prim(ivx,nx_max+1-ix,iy)  = min(u_prim(ivx,nx_max+1-ix,iy),0.0d0)
         u_prim(ivx,ix,iy)           = max(u_prim(ivx,ix,iy) ,0.0d0)
#if NDUST>0
        do idust=1,ndust
           u_prim(ivdx(idust),nx_max+1-ix,iy)  = min(u_prim(ivdx(idust),nx_max+1-ix,iy),0.0d0)
           u_prim(ivdx(idust),ix,iy)           = max(u_prim(ivdx(idust),ix,iy) ,0.0d0)      
        end do
#endif
      end do
      !stop
  end do

  !Accretion in the outer boundary
   if(accretion) then
   do ix = 1,nghost+1
    do iy=first_active_y,last_active_y
        if(phi(last_active,iy)<=phi_mom*pi) then
            maccreted = M_acc/(365.25*24.*3600.)*unit_t
            vrout     = cs(last_active,iy)
            sigma_out = maccreted/vrout/phi_mom/pi/box_l

            u_prim(irho,nx_max+1-ix,iy) =  sigma_out
            u_prim(ivx,nx_max+1-ix,iy)  =- sigma_out*vrout
            u_prim(ivy,nx_max+1-ix,iy)  =  sigma_out*L_out*sqrt(Mstar/100.)*100./box_l
#if NDUST>0
           do idust=1,ndust
    
             u_prim(irhod(idust),nx_max+1-ix,iy)  =  dust2gas*sigma_out
             u_prim(ivdx(idust),nx_max+1-ix,iy)   =  -dust2gas*sigma_out*vrout
             u_prim(ivdy(idust),nx_max+1-ix,iy)   =  dust2gas*sigma_out*L_out*sqrt(Mstar/100.)*100./box_l
#if NDUSTPSCAL>0
            !do ipscal=1,ndustpscal
            u_prim(idust_pscal(idust,1),nx_max+1-ix,iy)  = dust2gas*sigma_out*scut/unit_l
            !end do
#endif     
           end do
#endif
        endif
   end do
   end do
   endif
end subroutine apply_boundaries




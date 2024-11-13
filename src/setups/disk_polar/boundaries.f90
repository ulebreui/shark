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
         u_prim(ix,iy,ivar)          = u_prim(ix,last_active_y -nghost+iy,ivar) 
         u_prim(ix,ny_max+1-iy,ivar) = u_prim(ix,first_active_y+nghost-iy,ivar)
      end do

   end do

   !stop
  end do
   do ix= 1,nghost 
     do iy= 1,ny_max
         do ivar = 1,nvar
            u_prim(ix,iy,ivar)          = u_prim(first_active,iy,ivar) 
            u_prim(nx_max+1-ix,iy,ivar) = u_prim(last_active,iy,ivar)
         end do     
         u_prim(nx_max+1-ix,iy,ivx)  = min(u_prim(nx_max+1-ix,iy,ivx),0.0d0)
         u_prim(ix,iy,ivx)           = max(u_prim(ix,iy,ivx) ,0.0d0)
#if NDUST>0
        do idust=1,ndust
           u_prim(nx_max+1-ix,iy,ivdx(idust))  = min(u_prim(nx_max+1-ix,iy,ivdx(idust)),0.0d0)
           u_prim(ix,iy,ivdx(idust))           = max(u_prim(ix,iy,ivdx(idust)) ,0.0d0)      
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
            vrout     = cs(icell(last_active,iy))
            sigma_out = maccreted/vrout/phi_mom/pi/box_l

            u_prim(nx_max+1-ix,iy,irho) =  sigma_out
            u_prim(nx_max+1-ix,iy,ivx)  =- sigma_out*vrout
            u_prim(nx_max+1-ix,iy,ivy)  =  sigma_out*L_out*sqrt(Mstar/100.)*100./box_l
#if NDUST>0
           do idust=1,ndust
    
             u_prim(nx_max+1-ix,iy,irhod(idust))  =  dust2gas*sigma_out
             u_prim(nx_max+1-ix,iy,ivdx(idust))   =  -dust2gas*sigma_out*vrout
             u_prim(nx_max+1-ix,iy,ivdy(idust))   =  dust2gas*sigma_out*L_out*sqrt(Mstar/100.)*100./box_l
#if NDUSTPSCAL>0
            !do ipscal=1,ndustpscal
            u_prim(nx_max+1-ix,iy,idust_pscal(idust,1))  = dust2gas*sigma_out*scut/unit_l
            !end do
#endif     
           end do
#endif
        endif
   end do
   end do
   endif
end subroutine apply_boundaries




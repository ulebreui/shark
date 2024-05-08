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

    call boundary_collapse_1D

  !Accretion in the outer boundary
   if(accretion) then
   do ix = 1,nghost+1
    do iy=1,ny_max
    ibound_right   = icell(nx_max+1-ix,iy)
        if(phi(icell(last_active,iy))<=phi_mom*pi) then
            maccreted = M_acc/(365.25*24.*3600.)*unit_t
            vrout     = cs(icell(last_active,iy))
            sigma_out = maccreted/vrout/phi_mom/pi/box_l

            u_prim(ibound_right,irho) =  sigma_out
            u_prim(ibound_right,ivx)  =- sigma_out*vrout
            u_prim(ibound_right,ivy)  =  sigma_out*L_out*sqrt(Mstar/100.)*100./box_l
#if NDUST>0
           do idust=1,ndust
! #if NDUSTPSCAL>0
!             rhod_old=u_prim(ibound_right,irhod(idust))
! #endif           
             u_prim(ibound_right,irhod(idust))  =  dust2gas*sigma_out
             u_prim(ibound_right,ivdx(idust))   =  -dust2gas*sigma_out*vrout
             u_prim(ibound_right,ivdy(idust))   =  dust2gas*sigma_out*L_out*sqrt(Mstar/100.)*100./box_l
#if NDUSTPSCAL>0
            !do ipscal=1,ndustpscal
            u_prim(ibound_right,idust_pscal(idust,1))  = dust2gas*sigma_out*scut/unit_l
            !end do
#endif     
           end do
#endif
        endif
   end do
   end do
   endif
end subroutine apply_boundaries



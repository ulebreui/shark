subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud,perturbation,cs0,rr,r_cyl
  real(dp) :: rmax,vol_tot,xx,yy,x_soft,y_soft,Omega,H,cutoff,T_100
  real(dp) :: B_field,vfarg
  integer  :: i,idust,imax,ix,iy,icell,ivar

  call allocate_init
  allocate(uprim_condinit(1:ncells,1:nvar))
  uprim_condinit = 0.0d0

  ! Dimensionless units
  sigma_R0   = sigma_R0/unit_dcol
  tend = tend * (365.25*3600.*24.)/unit_t

  print *, 'tend = ', tend

  call gridinit_disk_log(box_l,smooth_r,box_l/3.)
  q      = 0.0d0
  iso_cs = 1

  do iy = 1,ny_max
    do ix = 1,nx_max
      xx       = position(icell(ix,iy),1)  ! Boxlen already in pc
      yy       = position(icell(ix,iy),2)-box_l/6.0
      rr       = radii_c(icell(ix,iy))
      H        = HoverR*rr
      Omega    = sqrt(Mstar/rr**3)
      cs0      = Omega*H
      vfarg    = 0.0d0
      if(fargo) then
        vfarg = -omega*rr
        fargo_velocity(icell(ix,iy))=vfarg
      endif
      q(icell(ix,iy),irho)                    = (sigma_R0/sqrt(2.0*pi*H))*(R0_disk/rr)*max(exp(-rr/R_out_disk)*exp(-yy**2/(2.0d0*H**2)),1e-3)
      q(icell(ix,iy),ivx)                     = 0.0
      q(icell(ix,iy),ivz)                     = omega*rr*sqrt(1.0d0-HoverR**2.0*(2.0d0+(rr/R_out_disk)))
      q(icell(ix,iy),iP)                      = q(icell(ix,iy),irho)*cs0**2.0
      cs(icell(ix,iy))                        = cs0
#if NDUST>0
     do idust=1,ndust
        q(icell(ix,iy),irhod(idust))    = dust2gas*q(icell(ix,iy),irho)
        epsilondust(icell(ix,iy),idust) = dust2gas
        sdust(icell(ix,iy),idust)       = scut/unit_l ! The 2H term comes from the vertical integration 
        sminstep=scut
#if NDUSTPSCAL>0
        if(growth_step)  q(icell(ix,iy),idust_pscal(idust,1)) = scut/unit_l
#endif
        q(icell(ix,iy),ivdx(idust))    = 0.0d0
        q(icell(ix,iy),ivdz(idust))    = omega*rr+vfarg
     end do
#endif
  eta_visc(icell(ix,iy)) = alpha_visc*H*cs(icell(ix,iy))
  enddo
  enddo
#if NDUST>0  
 ! call distribution_dust(.true.)
#endif

  call primtoc
  call update_force_setup
#if NDUST>0  
  call compute_tstop
#endif
do i=1,ncells
do ivar=1,nvar
    uprim_condinit(i,ivar)=u_prim(i,ivar)
    end do
end do
end subroutine setup



subroutine get_rhoturb(pert,del)

  use random
  use precision
  implicit none

  integer :: i
  integer::iseed=0         
  integer ,dimension(1,1:IRandNumSize)    :: allseed
  integer,dimension(IRandNumSize) :: localseed=-1
  real(dp) ::  pert,del,randno

  if (localseed(1)==-1) then
     call rans(1,iseed,allseed)
     localseed = allseed(1,1:IRandNumSize)
  end if

  call ranf(localseed,randno)

  del = randno*pert - 0.5d0*pert


end subroutine get_rhoturb 

subroutine write_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
end subroutine write_setup_info


subroutine read_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  real(dp):: info1
  
  read(ilun,*) info1
  time=info1
  close(ilun)
end subroutine read_setup_info


subroutine read_setup_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/setup_params/n_rel,box_l,r_relax,accretion,a_smooth_pl, sigma_R0,R_out_disk, R0_disk,Mstar , HoverR ,smooth_r, Mstar, Mplanet, Rplanet,alpha_visc,decrease_density,test_planet_disk

   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "Setup namelist reading  !"

   read(13,setup_params,IOSTAT=io)
   rewind(13)
   if (io/=0) then
      write(*,*) 'Invalid line in the setup namelist'
      stop
   end if

   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"

  
 end subroutine read_setup_params

 subroutine flag_continue(continue_sim)
   use parameters
   use commons
   use units
   implicit none
   logical :: continue_sim
   !Here you can add flags to kill the simulation
   if(time>=tend) continue_sim=.false.

 end subroutine flag_continue

 subroutine check_output(icount,iout,outputing,verbose)
   use parameters
   use commons
   use units
   implicit none
   integer :: icount,iout
   logical :: outputing,verbose

        !We make an output at a certain frequency rate of for specific values of the density. This can be tuned at will
     if(icount.eq.freq_out) then

        print *, "time =",time,' tend = ', tend
        verbose=.true.
        outputing=.true.
        icount=0
     endif

     if(time.eq.0) outputing=.true.   
     if(outputing)call output(iout)
     if(outputing) iout=iout+1
     if(outputing) print *, "Outputing data "
     if(outputing) print *, "Total mass is", sum(u_prim(:,irho))
     if(outputing) print *, "Total momentum is", sum(u_prim(:,ivx)+u_prim(:,ivy))

     outputing=.false.

 end subroutine check_output

 subroutine setup_preloop
   use parameters
   use commons
   use units
   implicit none
  
end subroutine setup_preloop

subroutine setup_inloop
   use parameters
   use commons
   use units
   use OMP_LIB

   implicit none
   integer :: i,ix,iy,idust,icell,ivar,ixx,iyy,ipscal
   real(dp):: unit_lout,sigma_out,t_rel,vrout,maccreted
   real(dp)::rho_new,vx_new,vy_new,vz_new
   real(dp)::rho_old,vx_old,vy_old,vz_old
   real(dp)::rho_init,vx_init,vy_init,vz_init

   return

   !Relaxation in the inner boundary
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,t_rel,rho_init,vx_init,vy_init,vz_init,rho_old,vx_old,vy_old,vz_old,rho_new,vx_new,vy_new,vz_new,ipscal)
  !$OMP DO 
   do i=1,ncells
    if(active_cell(i)==1) then
            if(radii_c(i)<r_relax) then
                t_rel=n_rel*2.0d0*pi/(dsqrt(Mstar/(radii_c(i))**3)) ! Relaxation timescale
                rho_init = uprim_condinit(i,irho)
                vx_init  = uprim_condinit(i,ivx)/rho_init
                vy_init  = uprim_condinit(i,ivy)/rho_init
                vz_init  = uprim_condinit(i,ivz)/rho_init
                rho_old  = u_prim(i,irho)
                vx_old   = u_prim(i,ivx)/rho_old
                vy_old   = u_prim(i,ivy)/rho_old
                vz_old   = u_prim(i,ivz)/rho_old
                rho_new  = (rho_init)  * (1.0d0-exp(-dt/t_rel)) + rho_old*exp(-dt/t_rel)
                vx_new   = ( vx_init)  * (1.0d0-exp(-dt/t_rel)) +  vx_old*exp(-dt/t_rel)
                vy_new   = ( vy_init)  * (1.0d0-exp(-dt/t_rel)) +  vy_old*exp(-dt/t_rel)
                vz_new   = ( vz_init)  * (1.0d0-exp(-dt/t_rel)) +  vz_old*exp(-dt/t_rel)
                u_prim(i,irho) = rho_new
                u_prim(i,ivx)  = rho_new*vx_new
                u_prim(i,ivy)  = rho_new*vy_new
                u_prim(i,ivz)  = rho_new*vz_new

#if NDUST>0
                do idust=1,ndust
                    rho_init= uprim_condinit(i,irhod(idust))
                    vx_init = uprim_condinit(i,ivdx(idust))/rho_init
                    vy_init = uprim_condinit(i,ivdy(idust))/rho_init
                    vz_init = uprim_condinit(i,ivdz(idust))/rho_init

                    rho_old = u_prim(i,irhod(idust))
                    vx_old  = u_prim(i,ivdx(idust))/rho_old
                    vy_old  = u_prim(i,ivdy(idust))/rho_old
                    vz_old  = u_prim(i,ivdx(idust))/rho_old

                    rho_new = (rho_init)  * (1.0d0-exp(-dt/t_rel)) + rho_old*exp(-dt/t_rel)
                    vx_new  = ( vx_init)  * (1.0d0-exp(-dt/t_rel)) +  vx_old*exp(-dt/t_rel)
                    vy_new  = ( vy_init)  * (1.0d0-exp(-dt/t_rel)) +  vy_old*exp(-dt/t_rel)
                    vz_new  = ( vz_init)  * (1.0d0-exp(-dt/t_rel)) +  vz_old*exp(-dt/t_rel)
                    u_prim(i,irhod(idust)) = rho_new
                    u_prim(i,ivdx(idust))  = rho_new*vx_new
                    u_prim(i,ivdy(idust))  = rho_new*vy_new
                    u_prim(i,ivdz(idust))  = rho_new*vz_new

#if NDUSTPSCAL>0
                do ipscal=1,ndustpscal
                    !rho_init= uprim_condinit(i,idust_pscal(idust,ipscal))
                    !rho_old = u_prim(i,idust_pscal(idust,ipscal))
                    !rho_new = (rho_init)  * (1.0d0-exp(-dt/t_rel)) + rho_old*exp(-dt/t_rel)
                    u_prim(i,idust_pscal(idust,ipscal)) = u_prim(i,idust_pscal(idust,ipscal))/rho_old*rho_new ! We don't relax passive scalars
                end do
#endif                
                end do
#endif                
            endif
        endif
    end do
   !$OMP END DO
   !$OMP END PARALLEL

   ! do ix=1,first_active
   !      do iy=first_active_y,last_active_y
   !              do ivar=1,nvar
   !                  u_prim(icell(ix,iy),ivar)=uprim_condinit(icell(first_active,iy),ivar)
   !              end do
   !      end do
   !  end do


end subroutine setup_inloop


 subroutine update_force_setup
   use parameters
   use commons
   use units
   use OMP_LIB
   use gravity_commons
   use slope_limiter
   implicit none

   integer :: i,idust,ix,iy,ixx,iyy,icell
   real(dp) :: xx, yy, rr, xx_p, yy_p, rr_p, theta_p,r_cyl,x_soft,y_soft,d_planet,px,py,p_r,p_theta,phi_loc,phi_p,r_sphe


  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,xx,yy,rr, xx_p, yy_p, rr_p,theta_p,r_cyl,x_soft,y_soft,px,py,p_r,p_theta,r_sphe)
  !$OMP DO
   do i=1,ncells
    if(active_cell(i)==1)then

        yy_p       = position(i,2)-box_l/6.0

        r_sphe = sqrt(radii_c(i)**2+yy_p**2)
        force(i,1)  = - Mstar/r_sphe**3 * radii_c(i)
        force(i,2)  = - Mstar/r_sphe**3 * yy_p
        force(i,3)  = 0.0d0
#if NDUST>0
        do idust=1,ndust
         force_dust(i,1,idust)  = - Mstar/r_sphe**3 * radii_c(i)
         force_dust(i,2,idust)  = - Mstar/r_sphe**3 * yy_p
         force_dust(i,3,idust)  = 0.0d0
        end do
#endif

#if NDUST>0
        do idust=1,ndust
         if(.not. drag) then
          force(i,1)  = force(i,1)  - q(i,irhod(idust))/q(i,irho)*(q(i,ivx) - q(i,ivdx(idust)))/tstop(i,idust)
          force(i,2)  = force(i,2)  - q(i,irhod(idust))/q(i,irho)*(q(i,ivy) - q(i,ivdy(idust)))/tstop(i,idust)
          force(i,3)  = force(i,3)  - q(i,irhod(idust))/q(i,irho)*(q(i,ivz) - q(i,ivdz(idust)))/tstop(i,idust)
          force_dust(i,1,idust)  = force_dust(i,1,idust)  + (q(i,ivx) - q(i,ivdx(idust)))/tstop(i,idust)
          force_dust(i,2,idust)  = force_dust(i,2,idust)  + (q(i,ivy) - q(i,ivdy(idust)))/tstop(i,idust)
          force_dust(i,3,idust)  = force_dust(i,3,idust)  + (q(i,ivz) - q(i,ivdz(idust)))/tstop(i,idust)

         endif
        end do
#endif
    endif

   end do
  !$OMP END DO  
  !$OMP END PARALLEL  
end subroutine update_force_setup



#if NDUST>0
! Dust stopping time
subroutine compute_tstop
  
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust
  real(dp):: xx, yy, H,rr,sloc,x_stokes,f_Stokes,St1,vdrift_turb,sd,nd,t_L
  !Re-calc distribution

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,H,rr,sloc,x_stokes,f_Stokes,St1,vdrift_turb,sd,nd,t_L)
  !$OMP DO
  do i=1,ncells
   if(active_cell(i)==1) then
        rr       = radii_c(i)
        H        = HoverR*rr
        t_L      =  sqrt(radii_c(i)**3) ! 1/Omega

        do idust=1,ndust
        sloc                 = sdust(i,idust)
#if NDUSTPSCAL>0
        if(growth_step) sloc = q(i,idust_pscal(idust,1)) 
#endif    
        tstop(i,idust) = sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(i,idust)/(q(i,irho)/(sqrt(2.0d0*pi)*H)*cs(i))
        x_stokes       = 1.0d0
        f_Stokes       = 3.2-1.0d0-x_stokes+2.0d0/(1.+x_stokes)*(1./2.6+x_stokes**3./(1.6+x_stokes))
        St1            = tstop(i,idust)/t_l
        vdrift_turb    = sqrt(alpha_turb)*cs(i)*dsqrt(f_Stokes*St1)
        sd             = q(i,idust_pscal(idust,1))
        nd             = q(i,irhod(idust))/(4./3.*pi*sd**3*rhograin/unit_d)/(sqrt(2.0d0*pi)*H)
        tcoag(i,idust) = 3.0d0/(pi*sd**2*nd*vdrift_turb)/min(1.0d0,-log10(vdrift_turb*unit_v/vfrag)/log10(5.))
        end do
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine compute_tstop
#endif

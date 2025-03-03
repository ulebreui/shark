subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud,perturbation,cs0,rr,r_cyl
  real(dp) :: rmax,xx,yy,x_soft,y_soft,Omega,H,cutoff,T_100
  real(dp) :: B_field,vfarg,ts_loc,eta_drift,Stokes_num
  integer  :: idust,imax,ix,iy,icell,ivar

  call allocate_init

  mcloud      = M_cloud*Msun
  r0          = 2.0d0/5.0d0*Grav*mcloud/(kB*T0)*mu_gas*mH*alpha_cloud
  rho_cloud   = mcloud/(4.0d0/3.0d0*pi*r0**3.)
  Omega_0     = sqrt(beta_cloud*3.0*Grav*mcloud/r0**3)*unit_t

  cs0         = sqrt(gamma*kB*T0/(mu_gas*mH))/unit_v
  sigma_crit  = sqrt(2.0*rho_crit*pi/Grav)*sqrt(gamma*kB*T0/(mu_gas*mH))

  rho_sink    = sigma_crit/unit_dcol
  print*, 'sigma_crit = ',sigma_crit, 'g/cm2'

  print*, 'r0 = ', r0/au, ' AU '
  print*, 'rho_cloud = ', rho_cloud, ' g/cm3 '
  print*, 'cs0 = ', cs0*unit_v,' cm/s'

  tend = tend * (365.25*3600.*24.)/unit_t

  print *, 'tend = ', tend

  r0        = r0/unit_l
  rho_cloud = rho_cloud/unit_d
#if GRIDSPACE==1
  call gridinit_disk_log(2.0*r0,smooth_r)
#endif
#if GRIDSPACE==2
  call gridinit_disk_log(2.0*r0,smooth_r,r_cut,nx/4)
#endif
  q      = 0.0d0
  !iso_cs = 1
  non_standard_eos =1
  do iy = 1,ny_max
    do ix = 1,nx_max
      xx       = position(ix,iy,1)  ! Boxlen already in pc
      yy       = position(ix,iy,2)
      rr       = radii(ix,iy)

      if(rr<r0) then
        q(irho,ix,iy) = (2.0d0*rho_cloud*sqrt(r0**2.-rr**2.0)*(1.0d0+rho_pert*cos(2*phi(ix,iy)))+rho_cloud*r0/100.0d0)
      else
        q(irho,ix,iy) = (rho_cloud*r0/100.0d0)!*(1.0d0+0.1d0*cos(2*phi))
      endif
      q(ivx,ix,iy)  = 0.d0
      q(ivy,ix,iy)  = Omega_0*rr
      q(ivz,ix,iy)  = 0.d0
      cs(ix,iy)     = cs0
      q(iP,ix,iy)   = q(irho,ix,iy)*cs(ix,iy)**2
      !print*, rr
  enddo
  enddo

  call primtoc
  call update_force_setup
#if NDUST>0  
  call compute_tstop
#endif

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
  namelist/setup_params/r0,rsink,smooth_r,sigma_0,tend,Omega_0,T0,sigma_crit,sigma_sink,Mstar,rout,r_cut,alpha_cloud,beta_cloud,rho_pert

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
   integer :: icount,iout,ix,iy
   logical :: outputing,verbose

        !We make an output at a certain frequency rate of for specific values of the density. This can be tuned at will
     if(icount.eq.freq_out.or.time.eq.0) then

        print *, "time =",time,' tend = ', tend
        verbose=.true.
        outputing=.true.
        icount=0
     endif

     if(time.eq.0) outputing=.true. 
     if(outputing) then
        Mtot = 0.0d0
        vol_tot=0.0d0
        do iy=first_active_y,last_active_y
            do ix=first_active,last_active
                Mtot    = Mtot + q(irho,ix,iy)*vol(ix,iy)
                vol_tot = vol_tot+vol(ix,iy)
            enddo
        enddo
     endif  
     if(outputing)call output(iout)
     if(outputing) iout=iout+1
     if(outputing) print *, "Outputing data "
     if(outputing) print *,'Stellar mass  is ', Mstar*unit_m/2e33, ' Msun'
     if(outputing) print *,'Total mass    is ', Mtot*unit_m/2e33, ' Msun'
     if(outputing) print *,'Total surface is ', vol_tot, ' au^2'

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
   integer :: iy,ix
   real(dp) :: T,delta_rho,cs0


   !do ix = first_active,first_active+3
   ! ix = first_active
   !  do iy=first_active_y,last_active_y
   !      delta_rho = q(irho,ix,iy)-sigma_sink/unit_dcol
   !      if(delta_rho>0) then
   !          Mstar = Mstar +delta_rho*vol(ix,iy)
   !          q(irho,ix,iy) = sigma_sink/unit_dcol !+ 0.9*delta_rho
   !          u_prim(irho,ix,iy) = q(irho,ix,iy)
   !          u_prim(ivx,ix,iy)  = q(irho,ix,iy)*q(ivx,ix,iy)
   !          u_prim(ivy,ix,iy)  = q(irho,ix,iy)*q(ivy,ix,iy)
   !          u_prim(ivz,ix,iy)  = q(irho,ix,iy)*q(ivz,ix,iy)
   !      endif
   !  enddo
   !end do
   do iy=first_active_y,last_active_y
        Mstar = Mstar -flux_x(irho,first_active,iy)*surf(first_active,iy,1)*dt
   enddo
   !call apply_boundaries
   !cs0       = sqrt(gamma*kB*T0/(mu_gas*mH))/unit_v

   ! !$omp parallel do default(shared) schedule(RUNTIME) private(ix,iy)
   !  do iy=1,ny_max
   !   do ix=1,nx_max
   !          !T = T0 *(1.0d0+(q(irho,ix,iy)*unit_dcol/sigma_crit)**(gamma-1.))
   !          cs(ix,iy)          = cs0 *sqrt(1.0d0+(q(irho,ix,iy)*unit_dcol/sigma_crit)**(gamma-1.))
   !      enddo
   !  enddo
end subroutine setup_inloop


 subroutine update_force_setup
   use parameters
   use commons
   use units
   use OMP_LIB
   use slope_limiter
   implicit none

   integer :: idust,ix,iy,ixx,iyy,icell
   real(dp) :: xx, yy, rr, xx_p, yy_p, rr_p, theta_p,r_cyl,x_soft,y_soft,d_planet,px,py,p_r,p_theta,phi_loc,phi_p,phi_lft,phi_rgt,rho_tot,dx_l,dx_r,dy_l,dy_r

   if(self_gravity) then
   
   ! We compute phi the gravitational potential
    phi_grav=0.0d0

   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy, ixx, iyy,rho_tot,xx,yy)
    do iy=first_active_y,last_active_y
    do ix=first_active,last_active
     phi_grav(ix,iy) =0.0d0
     xx = position(ix,iy,1)
     yy = position(ix,iy,2)

     do iyy=first_active_y,last_active_y
            do ixx=first_active,last_active
                rho_tot = q(irho,ixx,iyy)
#if NDUST>0
                do idust=1,ndust
                    rho_tot = rho_tot + q(irhod(idust),ixx,iyy)
                end do
#endif
                phi_grav(ix,iy) = phi_grav(ix,iy) - rho_tot*vol(ixx,iyy)&
                &/sqrt((xx-position(ixx,iyy,1))**2+(yy-position(ixx,iyy,2))**2+(dx(ixx,iyy,1)**2+dx(ixx,iyy,2)**2*radii(ixx,iyy)**2)*l_soft**2)      
            enddo
        enddo
    enddo
    enddo

    call  apply_boundaries_phi

   !$omp parallel do default(shared) schedule(RUNTIME) private(ix,iy,phi_lft,phi_rgt)
    do iy=first_active_y,last_active_y
     do ix=first_active,last_active

            phi_lft=2.0d0*(phi_grav(ix,iy)-phi_grav(ix-1,iy))/(dx(ix,iy,1)+dx(ix-1,iy,1))
            phi_rgt=2.0d0*(phi_grav(ix+1,iy)-phi_grav(ix,iy))/(dx(ix+1,iy,1)+dx(ix,iy,1))

            grad_phi_sg_x(ix,iy) = half*(phi_lft+phi_rgt)

            phi_lft=2.0d0*(phi_grav(ix,iy)-phi_grav(ix,iy-1))/(dx(ix,iy,2)+dx(ix,iy-1,2))
            phi_rgt=2.0d0*(phi_grav(ix,iy+1)-phi_grav(ix,iy))/(dx(ix,iy+1,2)+dx(ix,iy,2))

            grad_phi_sg_y(ix,iy) = half*(phi_lft+phi_rgt)/radii(ix,iy)
        enddo
    enddo
    endif
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy)
   do iy=first_active_y,last_active_y
    do ix=first_active,last_active
        if(self_gravity) then
            force_x(ix,iy)  = - Mstar/radii(ix,iy)**2-grad_phi_sg_x(ix,iy)
            force_y(ix,iy)  = - grad_phi_sg_y(ix,iy)
            force_z(ix,iy)  = 0.0d0
#if NDUST>0
            do idust=1,ndust
                force_dust_x(idust,ix,iy)  =  - Mstar/radii(ix,iy)**2 - grad_phi_sg_x(ix,iy)
                force_dust_y(idust,ix,iy)  =  - grad_phi_sg_y(ix,iy)
                force_dust_z(idust,ix,iy)  = 0.0d0
            end do
#endif
        else
            force_x(ix,iy)  = - Mstar/radii(ix,iy)**2
            force_y(ix,iy)  = 0.0d0
            force_z(ix,iy)  = 0.0d0
#if NDUST>0
            do idust=1,ndust
                force_dust_x(idust,ix,iy)  =  - Mstar/radii(ix,iy)**2
                force_dust_y(idust,ix,iy)  = 0.0d0
                force_dust_z(idust,ix,iy)  = 0.0d0
            end do
#endif
        endif

#if NDUST>0
        do idust=1,ndust
         if(.not. drag) then
          force_x(ix,iy)  = force_x(ix,iy)  - q(irhod(idust),ix,iy)/q(irho,ix,iy)*(q(ivx,ix,iy) - q(ivdx(idust),ix,iy))/tstop(idust,ix,iy)
          force_y(ix,iy)  = force_y(ix,iy)  - q(irhod(idust),ix,iy)/q(irho,ix,iy)*(q(ivy,ix,iy) - q(ivdy(idust),ix,iy))/tstop(idust,ix,iy)
          force_z(ix,iy)  = force_z(ix,iy)  - q(irhod(idust),ix,iy)/q(irho,ix,iy)*(q(ivz,ix,iy) - q(ivdz(idust),ix,iy))/tstop(idust,ix,iy)
          force_dust_x(idust,ix,iy)  = force_dust_x(idust,ix,iy)  + (q(ivx,ix,iy) - q(ivdx(idust),ix,iy))/tstop(idust,ix,iy)
          force_dust_y(idust,ix,iy)  = force_dust_y(idust,ix,iy)  + (q(ivy,ix,iy) - q(ivdy(idust),ix,iy))/tstop(idust,ix,iy)
          force_dust_z(idust,ix,iy)  = force_dust_z(idust,ix,iy)  + (q(ivz,ix,iy) - q(ivdz(idust),ix,iy))/tstop(idust,ix,iy)

         endif
        end do
#endif
    end do

   end do 

   call apply_boundaries_force
end subroutine update_force_setup



#if NDUST>0
! Dust stopping time
subroutine compute_tstop
  
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: idust,ix,iy,icell
  real(dp):: xx, yy, H,rr,sloc,x_stokes,f_Stokes,St1,vdrift_turb,sd,nd,t_L,Hd
  
  ! Re-calc distribution
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy, xx, yy, H,rr,sloc,x_stokes,f_Stokes,St1,vdrift_turb,sd,nd,t_L,Hd)
   do iy=first_active_y,last_active_y
    do ix=first_active,last_active    

         print *, 'Careful dust drag not correctly implemented'
         stop  
        end do
     end do
  end do


end subroutine compute_tstop
#endif

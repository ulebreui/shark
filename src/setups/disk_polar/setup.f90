subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud,cs0,rr,r_cyl
  real(dp) :: rmax,vol_tot,xx,yy,x_soft,y_soft,Omega,H,cutoff,T_100
  real(dp) :: B_field,vfarg,ts_loc,eta_drift,Stokes_num
  integer  :: idust,imax,ix,iy,icell,ivar
  real(dp) :: dphi_b
  real(dp) :: r_in_edge, r_out_edge, xi
  real(dp) :: del
  call allocate_init
  allocate(uprim_condinit(1:nx_max,1:ny_max,1:nvar))
  uprim_condinit = 0.0d0

  ! Dimensionless units
  sigma_R0   = sigma_R0/unit_dcol
  tend = tend * (365.25*3600.*24.)/unit_t

  print *, 'tend = ', tend

  call gridinit_disk_log(box_l,smooth_r,smooth_r)

  ! Damping rate R(r)/tau per ring (edges taken from the actual grid)
  allocate(damp_rate(1:nx_max)); damp_rate = 0.0d0
  r_in_edge  = radii(first_active,first_active_y) - half*dx(first_active,first_active_y,1)
  r_out_edge = radii(last_active, first_active_y) + half*dx(last_active, first_active_y,1)
  do ix = first_active, last_active
     rr = radii(ix,first_active_y)
     if (rr < damp_in_frac*r_in_edge) then
        xi = (damp_in_frac*r_in_edge - rr)/((damp_in_frac - 1.0d0)*r_in_edge)     ! 0 -> 1 towards r_in
        damp_rate(ix) = xi**2*sqrt(Mstar/r_in_edge**3)/(2.0d0*pi*damp_tau)
     else if (rr > damp_out_frac*r_out_edge) then
        xi = (rr - damp_out_frac*r_out_edge)/((1.0d0 - damp_out_frac)*r_out_edge) ! 0 -> 1 towards r_out
        damp_rate(ix) = xi**2*sqrt(Mstar/r_out_edge**3)/(2.0d0*pi*damp_tau)
     end if
  end do

  q      = 0.0d0
  iso_cs = 1

  do iy = 1,ny_max
    do ix = 1,nx_max
      xx       = position(ix,iy,1)  ! Boxlen already in pc
      yy       = position(ix,iy,2)
      rr       = radii(ix,iy)
      H        = HoverR*rr
      Omega    = sqrt(Mstar/rr**3)
      cs0      = Omega*H

      q(irho,ix,iy)                      = sigma_R0*(R0_disk/rr)*exp(-rr/R_out_disk)
      if(test_planet_disk) q(irho,ix,iy) = sigma_R0*(rr/R0_disk)**(-0.5d0)

      !if(rr>R_out_disk) q(irho,ix,iy)  = sigma_R0*(R0_disk/rr)*decrease_density
      q(ivx,ix,iy)                     = 0.0d0
      !q(ivy,ix,iy)                     = omega*rr*sqrt(1.0d0-HoverR**2.0*2.0d0)
      q(ivy,ix,iy)                     = omega*rr*sqrt(1.0d0-HoverR**2.0*(2.0d0+(rr/R_out_disk)))
      if(test_planet_disk)q(ivy,ix,iy) = omega*rr
      if (perturbation > 0.0d0) then
         call get_rhoturb(perturbation, del)
         q(irho,ix,iy) = q(irho,ix,iy)*(1.0d0 + del)
      end if
      q(iP,ix,iy)                      = q(irho,ix,iy)*cs0**2.0
      cs(ix,iy)                        = cs0
#if NDUST>0
     do idust=1,ndust
        q(irhod(idust),ix,iy)   = dust2gas*q(irho,ix,iy)
        if (bump_amp > 0.0d0) then   ! test: Gaussian dust bump advected by the rotation
           dphi_b = modulo(phi(ix,iy) - bump_phi + pi, 2.0d0*pi) - pi
           q(irhod(idust),ix,iy) = q(irhod(idust),ix,iy)*(1.0d0 + bump_amp* &
                &  exp(-((rr - bump_r)**2 + (bump_r*dphi_b)**2)/(2.0d0*bump_w**2)))
        end if     
        !if(rr>R_out_disk) q(irhod(idust),ix,iy)   = dust2gas*q(irho,ix,iy)/100.0d0
        epsilondust(icell(ix,iy),idust)= dust2gas
        sdust(idust)      = scut/unit_l ! The 2H term comes from the vertical integration 
        sminstep=scut
#if NDUSTPSCAL>0
        if(growth_step)  q(idust_pscal(idust,1),ix,iy) = scut/unit_l
#endif
        ts_loc=  sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(idust)/(q(irho,ix,iy)/(sqrt(2.0d0*pi)*H)*cs(ix,iy))
        eta_drift = sqrt(1.0d0-HoverR**2.0*(2.0d0+(rr/R_out_disk)))-1.0d0
        Stokes_num= ts_loc * omega
        q(ivdx(idust),ix,iy) = 2.0d0*eta_drift*omega*rr/(Stokes_num+1.0d0/Stokes_num)
        q(ivdy(idust),ix,iy)    = omega*rr
     end do
#endif
  enddo
  enddo
#if NDUST>0  
 ! call distribution_dust
#endif

  call primtoc
  call update_force_setup
#if NDUST>0  
  call compute_tstop
#endif
  do iy = 1,ny_max
    do ix = 1,nx_max
    do ivar=1,nvar
            uprim_condinit(ix,iy,ivar)=u_prim(ivar,ix,iy)
        end do
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
  namelist/setup_params/perturbation,damping, damp_in_frac, damp_out_frac, damp_tau,indirect_term,freq_torque,n_rel,box_l,r_relax,accretion,a_smooth_pl, sigma_R0,R_out_disk, R0_disk,Mstar , HoverR ,smooth_r, Mstar, Mplanet, Rplanet,alpha_visc,decrease_density,test_planet_disk,bump_amp,bump_r,bump_phi,bump_w,nplanet, phi_planet, nramp_planet

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
   integer, save :: istep_tq = 0
   real(dp) :: Mtot

        !We make an output at a certain frequency rate of for specific values of the density. This can be tuned at will
     if(icount.eq.freq_out) then

        print *, "time =",time,' tend = ', tend
        verbose=.true.
        outputing=.true.
        icount=0
     endif
     istep_tq = istep_tq + 1
     if (nplanet > 0 .and. mod(istep_tq - 1, freq_torque) == 0) call write_planet_torques
     if(time.eq.0) outputing=.true.   
     if(outputing)call output(iout)
     if(outputing) iout=iout+1
     if(outputing) print *, "Outputing data "
     if(outputing) then
        Mtot = 0.0d0
        do iy=first_active_y,last_active_y
            do ix=first_active,last_active
            Mtot = Mtot + q(irho,ix,iy)*vol(ix,iy)
            enddo
        enddo
     endif
     if(outputing) print *, "Total mass is", Mtot
     if(outputing) print *, "Total momentum is", sum(u_prim(ivx,:,:)+u_prim(ivy,:,:))
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
   integer :: ix,iy,idust,icell,ivar,ixx,iyy,ipscal
   real(dp):: unit_lout,sigma_out,t_rel,vrout,maccreted
   real(dp)::rho_new,vx_new,vy_new
   real(dp)::rho_old,vx_old,vy_old
   real(dp)::rho_init,vx_init,vy_init
   real(dp) :: f, f1


   !Relaxation in the inner boundary

   !Relaxation in the inner boundary
   if (damping) then
   !$omp parallel do default(shared) schedule(static) &
   !$omp private(ix,iy,idust,ipscal,f,f1,rho_init,vx_init,vy_init,rho_old,vx_old,vy_old,rho_new,vx_new,vy_new)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         if (damp_rate(ix) <= 0.0d0) cycle
         f1 = dt*damp_rate(ix)/(1.0d0 + dt*damp_rate(ix))   ! implicit relaxation weight
         f  = 1.0d0 - f1

         rho_init = uprim_condinit(ix,iy,irho)
         vx_init  = uprim_condinit(ix,iy,ivx)/rho_init
         vy_init  = uprim_condinit(ix,iy,ivy)/rho_init

         rho_old  = u_prim(irho,ix,iy)
         vx_old   = u_prim(ivx,ix,iy)/rho_old
         vy_old   = u_prim(ivy,ix,iy)/rho_old

         rho_new  = rho_init*f1 + rho_old*f
         vx_new   = vx_init *f1 + vx_old *f
         vy_new   = vy_init *f1 + vy_old *f

         u_prim(irho,ix,iy) = rho_new
         u_prim(ivx,ix,iy)  = rho_new*vx_new
         u_prim(ivy,ix,iy)  = rho_new*vy_new
#if NDUST>0
         do idust = 1, ndust
            rho_init = uprim_condinit(ix,iy,irhod(idust))
            vx_init  = uprim_condinit(ix,iy,ivdx(idust))/rho_init
            vy_init  = uprim_condinit(ix,iy,ivdy(idust))/rho_init

            rho_old  = u_prim(irhod(idust),ix,iy)
            vx_old   = u_prim(ivdx(idust),ix,iy)/rho_old
            vy_old   = u_prim(ivdy(idust),ix,iy)/rho_old

            rho_new  = rho_init*f1 + rho_old*f
            vx_new   = vx_init *f1 + vx_old *f
            vy_new   = vy_init *f1 + vy_old *f

            u_prim(irhod(idust),ix,iy) = rho_new
            u_prim(ivdx(idust),ix,iy)  = rho_new*vx_new
            u_prim(ivdy(idust),ix,iy)  = rho_new*vy_new
#if NDUSTPSCAL>0
            do ipscal = 1, ndustpscal  ! passive scalars are not relaxed
               u_prim(idust_pscal(idust,ipscal),ix,iy) = u_prim(idust_pscal(idust,ipscal),ix,iy)/rho_old*rho_new
            end do
#endif
         end do
#endif
      end do
   end do
   !$omp end parallel do

  endif

end subroutine setup_inloop


 subroutine update_force_setup
   use parameters
   use commons
   use units
   use OMP_LIB
   use slope_limiter
   implicit none

   integer :: idust,ix,iy,ixx,iyy,icell
   real(dp) :: xx, yy, rr, xx_p, yy_p, rr_p, theta_p,r_cyl,x_soft,y_soft,d_planet,px,py,p_r,p_theta,phi_loc,phi_p
   integer :: ipl
   real(dp) :: dphi, d2, fac, fr, fphi
   real(dp), dimension(nplanet_max) :: phip, mp, eps2

    do ipl = 1, nplanet
      call planet_state(ipl, time, phip(ipl), mp(ipl))
      eps2(ipl) = (a_smooth_pl*HoverR*Rplanet(ipl))**2     ! Plummer softening, eps = a_smooth_pl*H(r_p)
    end do
   !$omp parallel do default(shared) schedule(static) &
   !$omp private(ix,iy,idust,ipl,dphi,d2,fac,fr,fphi)
   do iy = 1, ny_max
    do ix = 1, nx_max
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
        do ipl = 1, nplanet                    ! planet gravity (polar components)
           dphi = phi(ix,iy) - phip(ipl)
           d2   = radii(ix,iy)**2 + Rplanet(ipl)**2 - 2.0d0*radii(ix,iy)*Rplanet(ipl)*cos(dphi) + eps2(ipl)
           fac  = mp(ipl)/(d2*sqrt(d2))
           fr   = -fac*(radii(ix,iy) - Rplanet(ipl)*cos(dphi))
           fphi = -fac*Rplanet(ipl)*sin(dphi)
           if (indirect_term) then   ! star accelerated by the planet (frame centred on the star)
              fr   = fr   - mp(ipl)/Rplanet(ipl)**2*cos(dphi)
              fphi = fphi + mp(ipl)/Rplanet(ipl)**2*sin(dphi)
           end if           
           force_x(ix,iy) = force_x(ix,iy) + fr
           force_y(ix,iy) = force_y(ix,iy) + fphi
#if NDUST>0
           do idust = 1, ndust
              force_dust_x(idust,ix,iy) = force_dust_x(idust,ix,iy) + fr
              force_dust_y(idust,ix,iy) = force_dust_y(idust,ix,iy) + fphi
           end do
#endif
        end do

#if NDUST>0
        if (ix >= first_active .and. ix <= last_active .and. iy >= first_active_y .and. iy <= last_active_y) then
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
        end if
#endif

    end do
   end do 
   !$omp end parallel do
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

  !$omp parallel do default(shared) schedule(static) &
  !$omp private(ix,iy,idust,rr,H,t_L,sloc,x_stokes,f_Stokes,St1,vdrift_turb,sd,Hd,nd)
   do iy=first_active_y,last_active_y
    do ix=first_active,last_active
        rr       = radii(ix,iy)
        H        = HoverR*rr
        t_L      =  sqrt(radii(ix,iy)**3) ! 1/Omega

        do idust=1,ndust
        sloc                 = sdust(idust)
#if NDUSTPSCAL>0
        if(growth_step) sloc = q(idust_pscal(idust,1),ix,iy) 
#endif    
        tstop(idust,ix,iy) = sqrt(pi/8.0d0)*(rhograin/unit_d)*sloc/(q(irho,ix,iy)/(sqrt(2.0d0*pi)*H)*cs(ix,iy))
#if NDUSTPSCAL>0
        if(growth_step) then

            x_stokes       = 1.0d0
            f_Stokes       = 3.2-1.0d0-x_stokes+2.0d0/(1.+x_stokes)*(1./2.6+x_stokes**3./(1.6+x_stokes))
            St1            = tstop(idust,ix,iy)/t_l
            vdrift_turb    = sqrt(alpha_turb)*cs(ix,iy)*dsqrt(f_Stokes*St1)
            sd             = q(idust_pscal(idust,1),ix,iy)
            Hd             = H*min(1.0d0,sqrt(alpha_turb/(min(St1,0.5d0/(1.0d0+St1**2)))))
            nd             = q(irhod(idust),ix,iy)/(4./3.*pi*sd**3*rhograin/unit_d)/(sqrt(2.0d0*pi)*Hd)
            tcoag(ix,iy,idust) = 3.0d0/(pi*sd**2*nd*vdrift_turb)/min(1.0d0,-log10(vdrift_turb*unit_v/vfrag)/log10(5.))
            sfrag(ix,iy,idust) = sd*(vfrag/(vdrift_turb*unit_v))**2
        endif
#endif        
        end do
     end do
  end do
  !$omp end parallel do




end subroutine compute_tstop
#endif


subroutine planet_state(ipl, t, phip, mp)
  use parameters
  use commons
  implicit none
  integer,  intent(in)  :: ipl
  real(dp), intent(in)  :: t
  real(dp), intent(out) :: phip, mp
  real(dp) :: omp, tramp
  omp   = sqrt((Mstar + Mplanet(ipl))/Rplanet(ipl)**3)
  phip  = modulo(phi_planet(ipl) + omp*t, 2.0d0*pi)
  tramp = nramp_planet(ipl)*2.0d0*pi/omp
  mp    = Mplanet(ipl)
  if (tramp > 0.0d0) mp = Mplanet(ipl)*sin(0.5d0*pi*min(t/tramp, 1.0d0))**2
end subroutine planet_state

! Torque of the gas and dust on each planet (about the star), appended to planet_torques.dat
subroutine write_planet_torques
  use parameters
  use commons
  use units
  implicit none
  integer  :: ipl, ix, iy, idust
  real(dp) :: phip, mp, eps2, dphi, d2, aphi, tq_gas, tq_dust
  if (nplanet == 0) return
  open(31, file='planet_torques.dat', position='append')
  do ipl = 1, nplanet
     call planet_state(ipl, time, phip, mp)
     eps2 = (a_smooth_pl*HoverR*Rplanet(ipl))**2
     tq_gas = 0.0d0; tq_dust = 0.0d0
     do iy = first_active_y, last_active_y
        do ix = first_active, last_active
           dphi = phi(ix,iy) - phip
           d2   = radii(ix,iy)**2 + Rplanet(ipl)**2 - 2.0d0*radii(ix,iy)*Rplanet(ipl)*cos(dphi) + eps2
           aphi = -mp*Rplanet(ipl)*sin(dphi)/(d2*sqrt(d2))     ! azimuthal acceleration of the cell
           ! central force => torque on the planet = - torque on the cell
           tq_gas = tq_gas - u_prim(irho,ix,iy)*vol(ix,iy)*radii(ix,iy)*aphi
#if NDUST>0
           do idust = 1, ndust
              tq_dust = tq_dust - u_prim(irhod(idust),ix,iy)*vol(ix,iy)*radii(ix,iy)*aphi
           end do
#endif
        end do
     end do
     write(31,'(E23.15,I4,3E23.15)') time, ipl, mp, tq_gas, tq_dust
  end do
  close(31)
end subroutine write_planet_torques

! Gas scale height (used by the vertically integrated coagulation kernel)
function gas_scale_height(ix,iy)
  use parameters
  use commons
  implicit none
  integer  :: ix, iy
  real(dp) :: gas_scale_height
  gas_scale_height = HoverR*radii(ix,iy)
end function gas_scale_height
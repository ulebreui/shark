subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud
  real(dp) :: rmax,vol_tot
  real(dp) :: B_field

  integer :: i,idust,imax
  zero_flux_in=.true. ! Mandatory for collapse calculations

  call allocate_init
  mcloud      = M_cloud*Msun
  r_cloud     = 2.0d0/5.0d0*Grav*mcloud/(kB*T0_cloud)*mu_gas*mH*alpha_cloud
  rmax        = nr_cloud*r_cloud
  rho_cloud   = mcloud/(4.0d0/3.0d0*pi*r_cloud**3.)/unit_d


  call gridinit(rmax)
 print *, r(last_active),rmax/unit_l,r_cloud/unit_l
 !initialisation of the hydro variables
 if(restarting.eq.0) then
  do i=first_active,last_active
     q(i,irho)=rho_cloud
     !if(r(i)>r_cloud/unit_l)q(i,irho)=rho_cloud/100.0d0
     q(i,iv)=0.0d0
  end do
#if NDUST>0
  call distribution_dust(.true.)
  do i=1,ncells
     do idust=1,ndust
        q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho)
     end do
  end do
#endif
else
   call read_output(nrestart)
   call distribution_dust(.true.)
   do i=1,ncells
      do idust=1,ndust
         q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho)
      end do
   end do
   restarting=0
endif

  call primtoc

  !The velocity is initially 0
  tend =nff*sqrt(3.*pi/(32.0*Grav*rho_cloud*unit_d))/unit_t

  call Mtot
  print *, 'Total mass is',M_tot*unit_m/2d33, 'Solar masses'
#if NDUST>0
  if(restarting.eq.0) then 
     call compute_B
  endif
  call charge 
#endif
  
end subroutine setup

subroutine gridinit(rmax)
  use parameters
  use commons
  use units
  implicit none

  real(dp):: rmax,zeta_r,ms
  integer :: i

  !Grid initialisation
  print *,' Logarithmic grid'
  zeta_r=(rmax/(rin*au))**(1.0d0/(ncells-2*nghost-1))
  r=0.0d0
  do i=1,last_active
     r(nghost+i)= (rin*au/unit_l)*zeta_r**(i-1)
  enddo

  do i=first_active,last_active
     dx(i)=r(i)-r(i-1)
  end do

  do i=1,first_active
     dx(i)=dx(first_active)
  end do

  !Cell center
  do i=first_active,ncells!-nghost
     r_c(i)     = ( (r(i)**3 + r(i-1)**3) / 2.)**(1./3.)
  end do
  !r_c(inner_bound)=r_c(first_active)  

 !Cell volume
  do i=first_active,last_active
     vol_cell(i) = (r(i)**3.-r(i-1)**3.)/3.0d0
  end do

  !Cell volume from center to the right 
  do i=first_active,last_active
     dvol(i)  =(r(i)**3-(r_c(i)**3) )/3.0d0
   end do

  do i=first_active,last_active
     Surf_p(i) = r(i)**2.
  end do

  Surf_p(inner_bound)=0.0
  Surf_p(last_active)=0.0

  !Needed for 2nd order
  
  do i=first_active,ncells-1!-nghost
     dx_c(i)      = r_c(i+1)- r_c(i-1)
     dx_r(i)      = r_c(i+1)- r_c(i)
     dx_l(i)      = r_c(i)- r_c(i-1)
     dx_r_cell(i) = r_c(i+1)- r(i)
     dx_l_cell(i) = r(i)- r_c(i)
  end do
  do i=1,nghost
     dx_l(i)      = r_c(first_active)
     dx_l_cell(i) = r_c(first_active)
     dx_r(i)      = dx_r(first_active)
     dx_r_cell(i) = dx_r_cell(first_active)
     !dx_r(ncells+1-i)=dx_r(last_active)
     !dx_r_cell(ncells+1-i)=dx_r(last_active)*0.5d0
     !dx_l_cell(ncells+1-i)=dx_l(last_active)*0.5d0
  end do
  
  
  !Not dx_l etc at the boundary dont really matter, they are not really used
  print *, 'r'
  print *, r
  print *, 'rc'
  print *, r_c
  print *, 'dx'
  print *, dx
  print *, 'vol_cell'
  print *, vol_cell
  print *, 'dvol'
  print *, dvol
  print *, 'dxc'
  print *, dx_c 
  print *, 'dxr'
  print *, dx_r
  print *, 'dxl'
  print *, dx_l
  print *, 'dxrcell'
  print *, dx_r_cell
  print *, 'dxlcell'
  print *, dx_l_cell
  !stop
end subroutine gridinit


subroutine write_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  write(ilun,*) time*unit_t/(365.*24.*3600.*1000.)
  write(ilun,*) alpha_cloud
end subroutine write_setup_info


subroutine read_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  real(dp):: info1
  
  read(ilun,*) info1
  time=info1*(365.*24.*3600.*1000.)/unit_t
  close(ilun)
end subroutine read_setup_info


subroutine read_setup_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/setup_params/M_cloud,alpha_cloud,T0_cloud,gamma,grid_type,ramses_baro,rho_max_sim,more_outputs
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "Setup namelist reading  !"
   read(13,setup_params,IOSTAT=io)
   rewind(13)
   if (io/=0) then
      write(*,*) 'Invalid line in the setup namelist'
      stop
   end if
   
   print *,"Mass of the cloud: "     ,  M_cloud
   print *,"Virial parameter: "      ,  alpha_cloud
   print *,"Isothermal temperature: ",  T0_cloud
   print *,"Adiabatic index: "       ,  gamma
   print*, 'Simulation will end at density', rho_max_sim
   if(ramses_baro)print *,"Barotrop as in RAMSES"
   if(.not.ramses_baro)print *,"Full barotrop"
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
   if(time>=tend.and..not.stop_at_first_core) continue_sim=.false.
   if(stop_at_first_core.and.maxval(uold(:,irho)*unit_d)>rho_max_sim) continue_sim=.false.

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
        print *, "time =",time*unit_t/(365.*24.*3600.*1000.),' tend = ', tend*unit_t/(365.*24.*3600.*1000.), ' kyr'
        print *, "Max density = ",maxval(uold(:,irho)*unit_d)," g / cc"
        icount=0
        print *, 'Total gas mass is',M_tot*unit_m/2d33, 'Solar masses'
#if NDUST>0
        print *, 'Total dust mass is',M_tot_dust*unit_m/2d33, 'Solar masses'
#endif
        if(sink) print *, 'Sink mass is ', M_cent*unit_m/2d33, 'Solar masses'
        verbose=.true.
        if(maxval(uold(:,irho)*unit_d)>rho_max_sim.and..not.stop_at_first_core) more_outputs=.true.
        if(more_outputs)outputing=.true.
     endif
     if(maxval(uold(:,irho)*unit_d)>order_mag) then
        order_mag=order_mag*10.
        outputing=.true.
     endif
     if(time.eq.0) outputing=.true.
#if NDUST>0     
     if(outputing.and.charging) call charge 
#endif     
     if(outputing)call output(iout)
     if(outputing) iout=iout+1
     if(outputing) print *, "Outputing data for Max density = ",maxval(uold(:,irho)*unit_d)," g / cc"
     outputing=.false.

 end subroutine check_output

 subroutine setup_preloop
   use parameters
   use commons
   use units
   implicit none

   order_mag=1d-24
   ! The idea here is to dump an output each time the density gains an order of magnitude. Since we don't know the initial density of the cloud here we find
   ! the next order of magnitude iteratively
  do while(order_mag<maxval(uold(:,irho)*unit_d))
     order_mag=order_mag*10.
  end do
  print *, "max density", maxval(uold(:,irho)*unit_d), 'order mag' ,order_mag
  
end subroutine setup_preloop

subroutine setup_inloop
   use parameters
   use commons
   use units
   implicit none
   !call output(1)
   return
end subroutine setup_inloop

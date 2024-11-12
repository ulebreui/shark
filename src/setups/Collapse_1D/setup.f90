subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: barotrop, cs_eos, tff0,mcloud,rcloud,rho_cloud,r_cloud

  integer :: i,idust,imax,ix,iy,icell,itimes_part


  call allocate_init

#if NX>1
  mcloud      = M_cloud*Msun
  r_cloud     = 2.0d0/5.0d0*Grav*mcloud/(kB*T0_cloud)*mu_gas*mH*alpha_cloud
  rho_cloud   = mcloud/(4.0d0/3.0d0*pi*r_cloud**3.)/unit_d
  tff0        = sqrt(3.*pi/32./grav/(rho_cloud*unit_d))/unit_t

  tend        = tend*tff0
  call gridinit_sphere1D(r_cloud)
  q=0.0d0
  non_standard_eos=1

  do iy = 1,ny_max
    do ix = 1,nx_max
      q(icell(ix,iy),irho) = rho_cloud
      q(icell(ix,iy),ivx)  = 0.0d0
      q(icell(ix,iy),iP)   = rho_cloud*cs_eos(barotrop(rho_cloud))**2
    end do
  end do

#if NDUST>0
  call distribution_dust(.true.)
  do i=1,ncells
     do idust=1,ndust
        q(i,irhod(idust))   = epsilondust(i,idust)*q(i,irho)
        if(single_size) then
          sdust(i,idust)    = smax/unit_l
          q(i,irhod(idust)) = 0.01*q(i,irho)
        endif
     end do
  end do
#endif

  call apply_boundaries
  call primtoc
#if GRAVITY==1  
  call mtot
#endif 

  !if (dust_inertia) call resistivities_with_dust_inertia
  print *, 'Total mass is',M_tot*unit_m/2d33, 'Solar masses'
#else
  non_standard_eos=1
  static=.true.
  open(15,file=trim('single_zone_time.dat'))
  open(16,file=trim('single_zone_density.dat'))

  !read(15,*), ntimes_part
  allocate(times_particles(1:ntimes_part))
  allocate(rho_particles(1:ntimes_part))
  ! Here we need to read the collapse table

  do itimes_part=1,ntimes_part
    read(15,*) times_particles(itimes_part)
    read(16,*) rho_particles(itimes_part)
  end do
  close(15)
  close(16)
  times_particles=times_particles/unit_t
  rho_particles  =rho_particles

  print *, times_particles
  print *, rho_particles
  !stop
  !/!\ Carefull we will need the correct unit for this ! /!\  
  tend = maxval(times_particles)
  time = times_particles(1) ! Set the initial time
  q=0.d0
#if NDUST>0
  call distribution_dust(.true.)
#endif
  do i=1,ncells
        q(i,irho)= rho_particles(ilist)
        q(i,ivx) = 0.0d0
        q(i,ivy) = 0.0d0
        q(i,ivx) = 0.0d0
#if NDUST>0                
     do idust=1,ndust
        q(i,irhod(idust))= q(i,irho)*epsilondust(i,idust)
        q(i,ivdx(idust)) = 0.0d0
        q(i,ivdy(idust)) = 0.0d0
        q(i,ivdx(idust)) = 0.0d0        
     end do
#endif     
  end do
  call primtoc
  !if(charging) call charge
  if(charging) call analytical_charge

  if (dust_inertia) call resistivities_with_dust_inertia
#endif
end subroutine setup


subroutine write_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  write(ilun,*) time
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
  namelist/setup_params/M_cloud,alpha_cloud,more_outputs,single_size,ntimes_part
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
#if NX>1
   !Here you can add flags to kill the simulation
   !if(time>=tend) continue_sim=.false.
   if(maxval(u_prim(:,irho)*unit_d)>rho_max_sim) continue_sim=.false.
#else

   dt=(times_particles(ilist+1)-times_particles(ilist))

   ilist=ilist+1

   if(ilist>ntimes_part) then

      continue_sim = .false.

   endif

#endif
 end subroutine flag_continue

 subroutine check_output(icount,iout,outputing,verbose)
   use parameters
   use commons
   use units
   implicit none
   integer :: icount,iout
   logical :: outputing,verbose

#if NX>1   
    !We make an output at a certain frequency rate of for specific values of the density. This can be tuned at will
     if(icount.eq.freq_out) then
        verbose=.true.
        if(maxval(u_prim(:,irho)*unit_d)>rho_max_sim) more_outputs=.true.
        if(more_outputs)outputing=.true.
     endif
     if(maxval(u_prim(:,irho)*unit_d)>order_mag) then
        order_mag=order_mag*10.
        outputing=.true.
     endif
     if(time.eq.0) outputing=.true.
    
     if(outputing)call output(iout)
     if(outputing) iout=iout+1
     if(outputing) then 
      print *, "time = ",time*unit_t/(365.*24.*3600.*1000.),' tend = ', tend*unit_t/(365.*24.*3600.*1000.), ' kyr'
      print *, "Max density = ",maxval(u_prim(:,irho)*unit_d)," g / cc"
      icount=0
      print *, 'Total gas mass is',M_tot*unit_m/2d33, 'Solar masses'
#if NDUST>0
       print *, 'Total dust mass is',M_tot_dust*unit_m/2d33, 'Solar masses'
#endif
     print *, "Outputing data for Max density = ",maxval(u_prim(:,irho)*unit_d)," g / cc"
     endif
     outputing=.false.
#else
   outputing = .true.
   verbose   = .true.
   call output(iout)
   !iout=2
   iout=iout+1
#endif
 end subroutine check_output

 subroutine setup_preloop
   use parameters
   use commons
   use units
   implicit none
#if NX>1
    order_mag=1d-24
    ! The idea here is to dump an output each time the density gains an order of magnitude. Since we don't know the initial density of the cloud here we find
     ! the next order of magnitude iteratively
    do while(order_mag<maxval(u_prim(:,irho)*unit_d))
      order_mag=order_mag*10.
  end do
  print *, "max density", maxval(u_prim(:,irho)*unit_d), 'order mag' ,order_mag
#endif
end subroutine setup_preloop

subroutine setup_inloop
   use parameters
   use commons
   use units
   implicit none
   !call output(1)
   integer :: icells,idust,i

#if NX==1
    do i=1,ncells 
#if NDUST>0     
    do idust=1,ndust
        epsilondust(i,idust) = u_prim(i,irhod(idust))/u_prim(i,irho)
    end do
#endif      
    u_prim(i,irho) = rho_particles(ilist)
#if NDUST>0
    do idust=1,ndust
        u_prim(i,irhod(idust)) = u_prim(i,irho)*epsilondust(i,idust)
    end do
#endif
  end do
#endif
end subroutine setup_inloop

 subroutine update_force_setup
   use parameters
   use commons
   use units
   implicit none

    integer :: i, idust

     return
#if NDUST>0
    call compute_tstop
#endif

   do i=1,ncells
    if(active_cell(i)==1)then
#if NDUST>0
        do idust=1,ndust
         if(.not. drag) then
          force(i,1)             =  - q(i,irhod(idust))/q(i,irho)*(q(i,ivx) - q(i,ivdx(idust)))/tstop(i,idust)
          force_dust(i,1,idust)  =  (q(i,ivx) - q(i,ivdx(idust)))/tstop(i,idust)
         endif
        end do
#endif
    endif

   end do
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
  real(dp):: pn,rhon,cs_eos,barotrop
  !Re-calc distribution


  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do i=1,ncells
   if(active_cell(i)==1) then
     do idust=1,ndust
        tstop(i,idust) = sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(i,idust)/(q(i,irho)*cs_eos(barotrop(q(i,irho))))
     end do
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine compute_tstop
#endif


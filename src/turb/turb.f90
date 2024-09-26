

subroutine random_acceleration_and_velocity(velocity)

!!To be called once initially to compute amplitude and phase of each turb mode
  use commons
  use parameters
  use units
  use precision
  use random

  implicit none

  logical :: velocity
  integer :: i
  integer ,dimension(1,1:IRandNumSize)    :: allseed
  integer,dimension(IRandNumSize) :: localseed=-1
  real(dp) ::  rand_nb1,rand_nb2,rand_nb3,rand_nb4,rand_nb5,rand_nb6,rand_nb7,rand_nb8,rand_nb9


  rand_nb8=0.0d0
  !Fill array of wavenumbers k_turb_driven between k_min and _k_max
  call linspace(k_min,k_max,k_turb_driven)


  
   if (localseed(1)==-1) then

    if (new_seed) then !Change iseed

      call rans(1,iseed,allseed)
      print *, 'New seed, iseed = ',iseed
      localseed = allseed(1,1:IRandNumSize)
      open(1,file="random_seed_turb.dat",status="replace",action='write')
      write(1, *) localseed
      close(1)
    end if

    if (.not. new_seed) then
      print *, 'Read existing seed, iseed =', iseed
      open(2,file="random_seed_turb.dat",status="old",action="read")
      read(2, *) localseed
      close(2)

    end if  
  end if

  do i = 1,nb_turb_modes_driven

    call gaussdev(localseed,rand_nb1) !Generate random numbers from a Gaussian for acceleration amplitudes
    random_array_ax(i) = corrector*rand_nb1 !Corrector to tune a posteriori the amplitudes to approach desired Mach
    call gaussdev(localseed,rand_nb2)
    random_array_ay(i) = corrector_sol*rand_nb2
    call gaussdev(localseed,rand_nb3)
    random_array_az(i) = corrector_sol*rand_nb3
    call ranf(localseed,rand_nb4) !Generate random numbers from a uniform distrib for phases
    random_array_phix(i) = 2.0d0*pi*rand_nb4
    call ranf(localseed,rand_nb5)
    random_array_phiy(i) = 2.0d0*pi*rand_nb5
    call ranf(localseed,rand_nb6)
    random_array_phiz(i) = 2.0d0*pi*rand_nb6

  if (velocity) then

    call gaussdev(localseed,rand_nb7) !Generate random velocity modes for initial conditions
    random_array_vx(i) = rand_nb7 
    call gaussdev(localseed,rand_nb8)
    random_array_vy(i) = rand_nb8
    call gaussdev(localseed,rand_nb9)
    random_array_vz(i) = rand_nb9
  end if 



  end do
  return

  end subroutine random_acceleration_and_velocity

subroutine initial_velocity_correction 

  !!Tune the amplitudes of each initial velocity mode to 'meet' desired Mach
  !!Velocity modes are taken to be sinusoidal
  !!Corecctions based on total root mean squared velocity

  use commons
  use parameters
  use units
  use precision

  implicit none
  integer :: i,i_turb
    

  real(dp) :: V_rms_ini,vx_mean,first_term,second_term,third_term,vel_corrector,vel_yz_corrector,vel_xyz_corrector



  if (turb_compressive) then
    if(.not. turb_solenoidal) then
      vel_corrector = Mach_nb_ini*cs_0/V_rms
      random_array_vx = vel_corrector*random_array_vx
      print *,'Initial velocity corrector (compressive) estimated to: ',vel_corrector

    end if
  end if



  if (turb_solenoidal) then
    if (.not. turb_compressive) then
      vel_yz_corrector = Mach_nb_ini*cs_0/Vtot_rms
      random_array_vy = vel_yz_corrector*random_array_vy
      random_array_vz = vel_yz_corrector*random_array_vz 
      print *,'Initial velocity corrector (solenoidal-yz) estimated to: ',vel_yz_corrector

    end if
  end if
 
 



  if (turb_solenoidal .and. turb_compressive) then
    vel_xyz_corrector = Mach_nb_ini*cs_0/Vtot_rms
    random_array_vx = vel_xyz_corrector*random_array_vx
    random_array_vy = vel_xyz_corrector*random_array_vy 
    random_array_vz = vel_xyz_corrector*random_array_vz  
    print *,'Initial velocity corrector (solenoidal-xyz) estimated to: ',vel_xyz_corrector



  endif

  


  return



end subroutine initial_velocity_correction




subroutine compute_rms_velocity


  use commons
  use parameters
  use units
  use precision

  implicit none

  real(dp) :: lenght 
  real (dp) :: vx_mean,vy_mean,vz_mean
  integer :: i


  lenght = 0.0d0 !To be reset to 0
  vx_mean=0.0d0
  vy_mean=0.0d0
  vz_mean=0.0d0

  V_rms=0.0d0
  Vy_rms=0.0d0
  Vz_rms=0.0d0
  Vyz_rms=0.0d0
  Vtot_rms=0.0d0


  do i=1,ncells
    if(active_cell(i)==1)then !Exclude ghost cells

    lenght = lenght + dx(i,1)  !Here we assume delta_x to be constant 

    endif
  end do


  do i=1,ncells
    if(active_cell(i)==1)then

    vx_mean = vx_mean + q(i,ivx)*dx(i,1)/lenght
    vy_mean = vy_mean + q(i,ivy)*dx(i,1)/lenght   
    vz_mean = vz_mean + q(i,ivz)*dx(i,1)/lenght

    endif
  end do

  do i=1,ncells
    if(active_cell(i)==1)then

    V_rms = V_rms + (q(i,ivx)-vx_mean)**2*dx(i,1)/lenght !This is V_rms**2
    Vy_rms = Vy_rms + (q(i,ivy)-vy_mean)**2*dx(i,1)/lenght
    Vz_rms = Vz_rms + (q(i,ivz)-vz_mean)**2*dx(i,1)/lenght

    endif
  end do




  Vtot_rms = SQRT(V_rms+Vy_rms+Vz_rms) !Here V_rms is actually V_rms**2
  Vyz_rms = SQRT(Vy_rms+Vz_rms) !Here V_rms is actually V_rms**2
  V_rms = SQRT(V_rms) !Take sqrt to get right V_rms. We should probably define a new variable
  Vy_rms = SQRT(Vy_rms)
  Vz_rms = SQRT(Vz_rms)




  return





end subroutine compute_rms_velocity





subroutine add_driven_turb_kick

!!Update momentum to include turbulence kick


  use commons
  use parameters
  use units
  use precision

  implicit none

  integer :: i,i_turb,idust
  real(dp) ::  mom_x,mom_y,mom_z,rho,ax_kick,ay_kick,az_kick,mom_x_d,mom_y_d,mom_z_d,rho_d


  do i=1,ncells
    if(active_cell(i)==1) then

	  mom_x = u_prim(i,ivx)
    mom_y = u_prim(i,ivy)
    mom_z = u_prim(i,ivz)

	  rho = u_prim(i,irho)

    ax_kick = 0.0d0
    ay_kick = 0.0d0
    az_kick = 0.0d0

    if (turb_compressive) then
      do i_turb=1,nb_turb_modes_driven

    	 ax_kick = ax_kick + random_array_ax(i_turb)*sin(k_turb_driven(i_turb)*position(i,1) + random_array_phix(i_turb)) !Beware of the definition of k in your setup/nml
    	
      end do

      u_prim(i,ivx) = mom_x + rho*ax_kick*dt

      if (turb_dust) then

        do idust=1,ndust

          mom_x_d = u_prim(i,ivdx(idust))
          rho_d= u_prim(i,irhod(idust))

          u_prim(i,ivdx(idust)) = mom_x_d + rho_d*ax_kick*dt

        end do
      end if





    end if


    if (turb_solenoidal) then

      do i_turb=1,nb_turb_modes_driven

       ay_kick = ay_kick + random_array_ay(i_turb)*sin(k_turb_driven(i_turb)*position(i,1) + random_array_phiy(i_turb)) !Beware of the definition of k in your setup/nml
       az_kick = az_kick + random_array_az(i_turb)*sin(k_turb_driven(i_turb)*position(i,1) + random_array_phiz(i_turb)) !Beware of the definition of k in your setup/nml
      
      end do

      u_prim(i,ivy) = mom_y + rho*ay_kick*dt
      u_prim(i,ivz) = mom_z + rho*az_kick*dt

      if (turb_dust) then

        do idust=1,ndust

          mom_y_d = u_prim(i,ivdy(idust))
          mom_z_d = u_prim(i,ivdz(idust))
          rho_d= u_prim(i,irhod(idust))

          u_prim(i,ivdy(idust)) = mom_y_d + rho_d*ay_kick*dt
          u_prim(i,ivdz(idust)) = mom_z_d + rho_d*az_kick*dt


        end do
      end if


    end if 


    end if

  end do 



end subroutine add_driven_turb_kick



subroutine kick_phase_drift
!Allows the phase of each turb kick modes to slowly drift as time goes by

  use commons
  use parameters
  use units
  use precision
  use random

  implicit none

  integer :: i_turb
  integer ,dimension(1,1:IRandNumSize)    :: allseed
  integer,dimension(IRandNumSize) :: localseed=-1
  real(dp) ::  rand_nb1,rand_nb2,rand_nb3

  
   if (localseed(1)==-1) then


      call rans(1,iseed_phase_drift,allseed)
      localseed = allseed(1,1:IRandNumSize)
 
  end if


  do i_turb=1,nb_turb_modes_driven
    call ranf(localseed,rand_nb1) !Generate random numbers from a uniform distrib for phases
    random_array_phix(i_turb) = random_array_phix(i_turb) + 2.0d0*pi*rand_nb1*dt/turnover_time
    call ranf(localseed,rand_nb2)
    random_array_phiy(i_turb) = random_array_phiy(i_turb) + 2.0d0*pi*rand_nb2*dt/turnover_time
    call ranf(localseed,rand_nb3)
    random_array_phiz(i_turb) = random_array_phiz(i_turb)+ 2.0d0*pi*rand_nb3*dt/turnover_time


  end do

  !print *, 'rand_nb1', rand_nb1
  !print *, 'rand_nb2', rand_nb2
  !print *, 'rand_nb3', rand_nb3


end subroutine kick_phase_drift




subroutine adjust_yz_kick_intensity
!To prevent transversal velocities from endlessly building up, adjust corrector to match targeted Mach

  use commons
  use parameters
  use units
  use precision
  use random

  implicit none


 
  call compute_rms_velocity

  if (Vyz_rms/cs_0 > Mach_yz_target) then



      random_array_ay = random_array_ay*0.97  
      random_array_az = random_array_az*0.97 

  end if 

  if (V_rms/cs_0 > Mach_x_target) then



      random_array_ax = random_array_ax*0.999  



  end if
end subroutine adjust_yz_kick_intensity



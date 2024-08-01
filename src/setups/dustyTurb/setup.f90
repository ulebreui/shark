subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud
  real(dp) :: rmax,vol_tot
  real(dp) :: B_field

  integer :: i,idust,imax,ix,iy,icell,i_turb

  box_l = box_l/unit_l


  call allocate_init
  call gridinit(box_l,box_l)
  q=0.0d0
  iso_cs=1


! #if TURB>0
!   !!!Decaying turb!!!
!   if (decaying_turb_compressive) then
!     allocate(k_turb(1:nb_turb_modes))
!     allocate(vx_turb(1:nb_turb_modes))
!     allocate(phix_turb(1:nb_turb_modes))


!     !Make directories for each Mach needed. Define a Mach param in setup_commons

!     !print *, trim(decay_turb_random_path) // trim('/wavenumber_turb_modes.dat')

!     open(15,file=trim(decay_turb_random_path) // trim('/wavenumber_turb_modes.dat')) !Read python-generated files
!     open(16,file=trim(decay_turb_random_path) // trim('/vx_turb_modes.dat'))
!     open(19,file=trim(decay_turb_random_path) // trim('/phix_turb_modes.dat'))


!     do i_turb = 1,nb_turb_modes
!       read(15,*) k_turb(i_turb)
!       read(16,*) vx_turb(i_turb)
!       read(19,*) phix_turb(i_turb)
!     end do

!     close(15)
!     close(16)
!     close(19)

!   endif

!   if (decaying_turb_solenoidal) then
!     allocate(k_turb(1:nb_turb_modes))
!     allocate(vy_turb(1:nb_turb_modes))
!     allocate(vz_turb(1:nb_turb_modes))
!     allocate(phiy_turb(1:nb_turb_modes))
!     allocate(phiz_turb(1:nb_turb_modes))


!     open(15,file=trim(decay_turb_random_path) // trim('/wavenumber_turb_modes.dat'))
!     open(17,file=trim(decay_turb_random_path) // trim('/vy_turb_modes.dat'))
!     open(18,file=trim(decay_turb_random_path) // trim('/vz_turb_modes.dat'))
!     open(20,file=trim(decay_turb_random_path) // trim('/phiy_turb_modes.dat'))
!     open(21,file=trim(decay_turb_random_path) // trim('/phiz_turb_modes.dat'))

!     do i_turb = 1,nb_turb_modes

!       read(15,*) k_turb(i_turb)
!       read(17,*) vy_turb(i_turb)
!       read(18,*) vz_turb(i_turb)
!       read(20,*) phiy_turb(i_turb)
!       read(21,*) phiz_turb(i_turb)
!     end do

!     close(15)
!     close(17)
!     close(18)
!     close(20)
!     close(21)

!   endif
! #endif


!Generate initial conditions (velocity) for turb
#if TURB>0

    call random_acceleration_and_velocity(.true.)
#endif



#if NDUST>0

    call distribution_dust(.true.)

    do i=1,ncells
      do idust=1,ndust
        sdust(i,idust)    = smax/unit_l !if a single grain
      end do
    end do  






#endif

  do i = 1,ncells

      q(i,irho) = rho_0/unit_d 

      cs(i)=cs_0/unit_v

      q(i,iP)=q(i,irho)*cs(i)**2


#if TURB>0
       if (turb_compressive) then 

      !   do i_turb = 1,nb_turb_modes

      !     q(i,ivx) = q(i,ivx) + vx_turb(i_turb)/unit_v*(sin(2.0d0*pi*position(i,1)/(box_l)*k_turb(i_turb)+phix_turb(i_turb))) !Decaying turb (compressive modes) for the gas. Make sure Vrms/cs = Mach.

      !   end do 
      ! end if

      !if (driven_turb_compressive) then !Driven turb initial condition!

        do i_turb = 1,nb_turb_modes_driven


          q(i,ivx) = q(i,ivx) + random_array_vx(i_turb)/unit_v*(sin(position(i,1)*k_turb_driven(i_turb)+random_array_phix(i_turb))) !Driven turb (compressive modes) for the gas. Beware: k in the form 2pi/l 

        end do
      endif

#endif




#if NDUST>0



        do idust=1,ndust


            !OLD SETUP
            epsilondust(i,idust) = dust2gas_ratio(idust) !To remove
            !epsilondust(i,idust) = 1.0d0 
            


            q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho)



            !St_0(idust)=1d-2.0d0*10000-9999*1d-2.0d0*(idust-1)
            !St_0(idust)=1d-2+9999*1d-2.0d0*(idust-1)

            !sdust(i,idust)=St_0(idust)*box_l !TODO edit when adding coagulation


         !q(i,ivd(idust)) = 0.0d0
         !print *, "Vxd=", q(i,ivd(idust))*unit_v

        end do


#endif

#if MHD==1
if(beta_0>0) then

     q(i,iBx)=dsqrt(rho_0/unit_d)*cs(i)/sqrt(beta_0) !todo : display

endif
#endif
     


#if TURB>0

  ! if (decaying_turb_solenoidal) then

  !   do i_turb = 1,nb_turb_modes

  !     q(i,ivy) = q(i,ivy) + vy_turb(i_turb)/unit_v*(sin(2.0d0*pi*position(i,1)/(box_l)*k_turb(i_turb)+phiy_turb(i_turb))) !Solenoidal modes. Make sure sqrt(sum(vx_turb**2+vy_turb**2))/cs = Mach.

  !     q(i,ivz) = q(i,ivz) + vz_turb(i_turb)/unit_v*(cos(2.0d0*pi*position(i,1)/(box_l)*k_turb(i_turb)+phiz_turb(i_turb)))
    
  !   end do

  ! endif 


  if (turb_solenoidal) then
      do i_turb = 1,nb_turb_modes_driven

            q(i,ivy) = q(i,ivy) + random_array_vy(i_turb)/unit_v*(sin(position(i,1)*k_turb_driven(i_turb)+random_array_phiy(i_turb))) !Solenoidal modes.

            q(i,ivz) = q(i,ivz) + random_array_vz(i_turb)/unit_v*(cos(position(i,1)*k_turb_driven(i_turb)+random_array_phiz(i_turb)))
    

      end do

  endif

#endif




#if MHD==1
if(beta_0==0) then

     q(i,iBx)=0.0d0
     q(i,iBy)=0.0d0
     q(i,iBz)=0.0d0

endif

#endif 



  end do


!!!Correct initial velocity to meet desired Mach!!!
#if TURB>0

  print *,'rndom_vy', random_array_vy
  print *,'rndom_vz', random_array_vz


  call compute_rms_velocity  !Compute rms for initial velocity profile provided by random generation of random_vx_array
  print *,'vrms=',V_rms
  print *,'vy_rms=',Vy_rms
  print *,'vz_rms=',Vz_rms
  print *,'vtot_rms=',Vtot_rms

  call initial_velocity_correction  !Uses the rms_velocity recently calculated
  !print *, 'vx corrected',random_array_vx
  print *,'rndom_vy corrected', random_array_vy
  print *,'rndom_vz corrected', random_array_vz




    do i = 1,ncells

      if (turb_compressive) then 

        q(i,ivx) = 0.0d0 !!!Don't forget to set it to 0 before correcting!!!


        do i_turb = 1,nb_turb_modes_driven !Reset initial conditions with corrected random_array_vx

          q(i,ivx) = q(i,ivx) + random_array_vx(i_turb)/unit_v*(sin(position(i,1)*k_turb_driven(i_turb)+random_array_phix(i_turb))) !Driven turb (compressive modes) for the gas. Beware: k in the form 2pi/l 

        end do
      endif



      if (turb_solenoidal) then


        q(i,ivy) = 0.0d0
        q(i,ivz) = 0.0d0

        do i_turb = 1,nb_turb_modes_driven

          q(i,ivy) = q(i,ivy) + random_array_vy(i_turb)/unit_v*(sin(position(i,1)*k_turb_driven(i_turb)+random_array_phiy(i_turb))) !Solenoidal modes.
          q(i,ivz) = q(i,ivz) + random_array_vz(i_turb)/unit_v*(cos(position(i,1)*k_turb_driven(i_turb)+random_array_phiz(i_turb)))
    

        end do

      endif

    end do !i boucle



  call compute_rms_velocity  !Call it again to check if correction performed correctly (read it in output_0000)
  print *,"updated vrms = ", V_rms
  print *,"updated vyrms = ", Vy_rms
  print *,"updated vzrms = ", Vz_rms
  print *,"updated vtotrms = ", Vtot_rms


#endif


#if MHD==1

  do i=1,ncells
  
    q(i,iBy) = delta_B*q(i,iBx)/unit_B*(sin(position(i,1)*k_mag)) !Alfven perturbation
    q(i,iBz) = delta_B*q(i,iBx)/unit_B*(cos(position(i,1)*k_mag)) 

  end do
#endif



  tend = tend/unit_t


  call apply_boundaries
  call primtoc
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
  namelist/setup_params/box_l,rho_0,St_0,dust2gas_ratio,beta_0,cs_0,k_mag,delta_B,T_cloud
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
     if(outputing) print *, "Total momentum is", sum(u_prim(:,ivx)+u_prim(:,ivy)+u_prim(:,ivz))
     if(outputing) print *, "Total energy is", sum(u_prim(:,iP))

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
   implicit none
   !call output(1)
   return
end subroutine setup_inloop

 subroutine update_force_setup
   use parameters
   use commons
   use units
   implicit none
  
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
  real(dp):: pn,rhon
  !Re-calc distribution


  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do i=1,ncells
   if(active_cell(i)==1) then
     do idust=1,ndust
        tstop(i,idust) = St_0(idust)*box_l*rho_0/cs(i)/q(i,irho)
     end do
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine compute_tstop
#endif


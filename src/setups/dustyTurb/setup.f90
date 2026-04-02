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



!Generate initial conditions (velocity) for turb
#if TURB>0

    call random_acceleration_and_velocity(.true.)
#endif



  do i = 1,ncells

      q(i,irho) = rho_0/unit_d 

      cs(i)=cs_0/unit_v

      q(i,iP)=q(i,irho)*cs(i)**2

      q(i,ivx) = delta_vx/unit_v*(cos(position(i,1)*k_mag)) 


      q(i,ivy) = delta_vy/unit_v*(sin(position(i,1)*k_mag)) 
      q(i,ivz) = delta_vz/unit_v*(cos(position(i,1)*k_mag)) 

  end do




#if NDUST>1

    call distribution_dust(.true.) !creates sdust, epsilondust and mdust arrays

#endif



  do i=1,ncells
    do idust=1,ndust

#if NDUST>1

    St(i,idust) = dsqrt(pi/8) * rhograin * sdust(i,idust) / (q(i,irho) *  box_l)
    q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho)

#endif

#if NDUST==1
        sdust(i,idust)    = smax/unit_l !if a single grain
        mdust(i,idust)    = (4./3.*pi*smax**3*rhograin)/unit_m
        St(i,idust)       = St_0(idust)
        epsilondust(i,idust) = dust2gas_ratio(idust) !To remove
        q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho)


#if NDUSTPSCAL==1
        q(i,idust_pscal(idust,1)) = smax/unit_l
        epsilondust(i,idust) = dust2gas_ratio(idust) !To remove
        q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho)
#endif

#endif

        q(i,ivdx(idust)) = delta_vdx/unit_v*(cos(position(i,1)*k_mag)) 

        q(i,ivdy(idust)) = delta_vdy/unit_v*(sin(position(i,1)*k_mag)) 
        q(i,ivdz(idust)) = delta_vdz/unit_v*(cos(position(i,1)*k_mag)) !Alfven perturbation
            !OLD SETUP
        

#if DUST_PRESSURE==1
        q(i,iPd(idust))=q(i,irhod(idust))*(delta_dust_cs*cs(i))**2
#endif



      end do
    end do  




  do i = 1,ncells


#if TURB>0
       if (decaying_turb_compressive) then 

        do i_turb = 1,nb_turb_modes

          q(i,ivx) = q(i,ivx) + vx_turb(i_turb)/unit_v*(sin(2.0d0*pi*position(i,1)/(box_l)*k_turb(i_turb)+phix_turb(i_turb))) !Decaying turb (compressive modes) for the gas. Make sure Vrms/cs = Mach.

        end do 
      end if

      if (turb_compressive) then !Driven turb initial condition!

        do i_turb = 1,nb_turb_modes_driven


          q(i,ivx) = q(i,ivx) + random_array_vx(i_turb)/unit_v*(sin(position(i,1)*k_turb_driven(i_turb)+random_array_phix(i_turb))) !Driven turb (compressive modes) for the gas. Beware: k in the form 2pi/l 

        end do
      endif

#endif




#if MHD==1
if(beta_0>0) then

     q(i,iBx)=dsqrt(4*pi*rho_0/unit_d)*cs(i)/dsqrt(beta_0) !todo : display

endif
#endif
     


#if TURB>0

  if (decaying_turb_solenoidal) then

    do i_turb = 1,nb_turb_modes

      q(i,ivy) = q(i,ivy) + vy_turb(i_turb)/unit_v*(sin(2.0d0*pi*position(i,1)/(box_l)*k_turb(i_turb)+phiy_turb(i_turb))) !Solenoidal modes. Make sure sqrt(sum(vx_turb**2+vy_turb**2))/cs = Mach.

      q(i,ivz) = q(i,ivz) + vz_turb(i_turb)/unit_v*(cos(2.0d0*pi*position(i,1)/(box_l)*k_turb(i_turb)+phiz_turb(i_turb)))
    
    end do

  endif 


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


!!!Correct initial velocity to meet desired initial Mach!!!
#if TURB>0


  call compute_rms_velocity  !Compute rms for initial velocity profile provided by random generation of random_vx_array
  call initial_velocity_correction  !Uses the rms_velocity recently calculated

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
  namelist/setup_params/box_l,rho_0,St_0,dust2gas_ratio,beta_0,cs_0,k_mag,delta_B,T_cloud,delta_vdx,delta_vdy,delta_vdz,delta_vx,delta_vy,delta_vz
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

   if(nrestart>0) then 
    call restart_setup_quantities
  endif

  
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
  integer :: i,idust,ipscal
  real(dp):: pn,rhon
  !Re-calc distribution


  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do i=1,ncells
   if(active_cell(i)==1) then
     do idust=1,ndust
        tstop(i,idust) = St_0(idust)*box_l*rho_0/cs(i)/q(i,irho) !St_0 is a namelist input. Should match with the def of smax
        St(i,idust) = tstop(i,idust) * cs(i) / box_l

#if NDUSTPSCAL>0
    ! Dust growth via monodisperse approach
    ipscal = 1 !The first passive scalar is the grain size
    tstop(i,idust) = dsqrt(pi/8) * rhograin * q(i,idust_pscal(idust,ipscal)) / (q(i,irho) * cs(i))
    St(i,idust) = tstop(i,idust) * cs(i) / box_l
#endif
     end do
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine compute_tstop
#endif

#if NDUST>0
! Dust coagulation time
subroutine compute_tcoag
  
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust,ipscal
  real(dp):: mgrain_step,cs_modified,St_frag 
  !Re-calc distribution


  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do i=1,ncells
   if(active_cell(i)==1) then
     do idust=1,ndust

        ipscal = 1 !The first passive scalar is the grain size


        mgrain_step = 4./3. * pi * q(i,idust_pscal(idust,1))**3 * rhograin
        if (St(i,idust) <= 1 .and. St(i,idust) < 0.1) dv_ormel_step(i,idust) = dsqrt(1.97*alpha_turb)*cs(i)*dsqrt(St(i,idust)) !Regime II
        if (St(i,idust) <= 1 .and. St(i,idust) >= 0.1) dv_ormel_step(i,idust) = dsqrt(alpha_turb)*cs(i)*dsqrt(St(i,idust)) !Regime II
        if (St(i,idust) > 1) dv_ormel_step(i,idust) = dsqrt(alpha_turb)*cs(i)*dsqrt(2./(1.+St(i,idust))) !Regime III


        if (modified_Ormel .eqv. .true.) then !Modify soundspeed due to dust backreaction

          cs_modified = cs(i)/dsqrt(1+SUM((q(i,irhod(:))/q(i,irho))/(1+St(i,:))))

          if (St(i,idust) <= 1 .and. St(i,idust) < 0.1) dv_ormel_step(i,idust) = dsqrt(1.97*alpha_turb)*cs_modified*dsqrt(St(i,idust)) !Regime II         
          if (St(i,idust) <= 1 .and. St(i,idust) >= 0.1) dv_ormel_step(i,idust) = dsqrt(alpha_turb)*cs_modified*dsqrt(St(i,idust)) !Regime II
          if (St(i,idust) > 1) dv_ormel_step(i,idust) = dsqrt(alpha_turb)*cs_modified*dsqrt(2./(1.+St(i,idust))) !Regime III
          

        endif

        !t_coag -> 3*t_coag in equa diff of grain size
        tcoag(i,idust) =  3.0/(sqrt(8/3*pi)*pi*4*q(i,idust_pscal(idust,ipscal))**2*dv_ormel_step(i,idust)*q(i,irhod(idust))/mgrain_step)

        if (frag_step) then

          St_frag = vfrag**2/(alpha_turb*cs(i)**2)

          if (modified_Ormel .eqv. .true.) St_frag = vfrag**2/(alpha_turb*cs_modified**2)

          size_frag_Ormel(i,idust) = St_frag*q(i,irho)*cs(i)*(box_l/cs(i))/rhograin/sqrt(pi/8.) !Define the stepinski size from the Stokes provided in the namelist
          !print*, 'sfrag=',sfrag(idust,ix,iy)

        endif


     end do
  end if
end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine compute_tcoag
#endif


subroutine restart_setup_quantities
!!Retrieve seteup dependent quantities that are not in the uprim array!!

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,idust
  real(dp) :: xdp
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out

  do i=1,ncells
    do idust=1,ndust

#if NDUST==1     
      sdust(i,idust)    = smax / unit_l !A priori not needed because no growth and this line is already called in setup.
      mdust(i,idust)    = (4./3.*pi*smax**3*rhograin)/unit_m !sdust and mdust needed for chemical network.

#if NDUSTPSCAL==1 
! --> q(i,idust_pscal(idust,1)) has been recovered in the uprim vector of the restart output
!mdust is updated when calling dust_growth_stepinski. However, analytical charge (which needs mdust) is called before, thus:
      mdust(i,idust)    = (4./3.*pi*q(i,idust_pscal(idust,1))**3*rhograin)/unit_m
#endif
#endif

  end do
end do
!Note: cs(i) is computed at the end of predictor step, and used in Riemann solvers and courant.
!But courant called first --> we need to retrieve the soundspeed from restart output or recompute it. HERE: isothermal, thus using the initial value cs_0 is valid.
end subroutine restart_setup_quantities



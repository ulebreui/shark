subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud
  real(dp) :: rmax,vol_tot
  real(dp) :: B_field

  integer :: i,idust,imax,ix,iy,icell

  box_l = box_l/unit_l


  call allocate_init
  call gridinit(box_l,box_l)
  q=0.0d0
  iso_cs=1




#if NDUST>0

    call distribution_dust(.true.)

#endif

  do i = 1,ncells

      q(i,irho) = rho_0/unit_d !Spatially uniform gas density

      cs(i)=cs_0/unit_v

      q(i,iP)=q(i,irho)*cs(i)**2


      !q(i,irho)= delta_rho*rho_0/unit_d*sin(2.0d0*pi*position(i,1)/box_l*kx_wave) !Density perturbation for the dust

      q(i,ivx) = vx_0/unit_v*(sin(2.0d0*pi*position(i,1)/(box_l)*kx_wave)) !Wave-like velocity perturbation for the gas



#if NDUST>0



        do idust=1,ndust


            !OLD SETUP
            epsilondust(i,idust) = dust2gas_ratio(idust) !To remove
            !epsilondust(i,idust) = 1.0d0 
            


            q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho) !Spatially uniform dust density



            !sdust(i,idust)=St_0(idust)*box_l !TODO edit when adding coagulation

         

            q(i,ivdx(idust)) = vdx_0/unit_v*sin(2.0d0*pi*position(i,1)/box_l*kx_wave_d) !Velocity perturbation for the dust

            !q(i,irhod(idust))= delta_rho_d*epsilondust(i,idust)*q(i,irho)*sin(2.0d0*pi*position(i,1)/box_l*kx_wave_d) !Density perturbation for the dust


         !q(i,ivd(idust)) = 0.0d0
         !print *, "Vxd=", q(i,ivd(idust))*unit_v

        end do


#endif

#if MHD==1
if(beta_0>0) then

     q(i,iBx)=dsqrt(rho_0/unit_d)*cs(i)/sqrt(beta_0) !todo : display
     q(i,iBy)=delta_B*q(i,iBx)*sin(2.0d0*pi*position(i,1)/box_l*ky_wave_B) !Transversal magnetic perturbations
     q(i,iBz)=delta_B*q(i,iBx)*cos(2.0d0*pi*position(i,1)/box_l*kz_wave_B)

     !delta_Bnorm = q(i,iBy)**2.0d0+q(i,iBz)**2.0d0
    !print *, "Bx =", q(i,iBx)
    !print *, "By =", q(i,iBy)
    !print *, "Bz =", q(i,iBz)

     !print *, "Bnorm squared =", Bnorm2





    q(i,ivy) = vy_0/unit_v*(sin(2.0d0*pi*position(i,1)/(box_l)*ky_wave)) !Wave-like velocity perturbation for the gas

    q(i,ivz) = vz_0/unit_v*(cos(2.0d0*pi*position(i,1)/(box_l)*kz_wave)) !Wave-like velocity perturbation for the gas


#if NDUST>0
     do idust=1,ndust
          q(i,ivdy(idust)) = vdy_0/unit_v*sin(2.0d0*pi*position(i,1)/box_l*ky_wave_d) !Velocity perturbation for the dust
          q(i,ivdz(idust)) = vdz_0/unit_v*sin(2.0d0*pi*position(i,1)/box_l*kz_wave_d) !Velocity perturbation for the dust
     end do
#endif





endif !end of if (beta>0). Perhaps to be placed elsewhere in the future



if(beta_0==0) then

     q(i,iBx)=0.0d0
     q(i,iBy)=0.0d0
     q(i,iBz)=0.0d0
endif

#endif
  



  end do



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
  namelist/setup_params/box_l,rho_0,St_0,dust2gas_ratio,vdx_0,vdy_0,vdz_0,vx_0,vy_0,vz_0,ky_wave_B,kz_wave_B,kx_wave,ky_wave,kz_wave,kx_wave_d,ky_wave_d,kz_wave_d,beta_0,cs_0,delta_B,delta_rho,delta_rho_d
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
     if(outputing) print *, "Total momentum is", sum(u_prim(:,ivx))
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


subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud
  real(dp) :: rmax,vol_tot

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

      !q(i,irho) = rho_0/unit_d 
      q(i,irho) = rho_0/unit_d*(1 + delta_rho*sin(2.0d0*pi*position(i,1)/box_l*kx_wave)) !DUSTYWAVE
      q(i,ivx)  = vx_0/unit_v*sin(2.0d0*pi*position(i,1)/box_l*kx_wave) !DUSTYWAVE


      cs(i)=cs_0/unit_v
      q(i,iP) = cs(i)**2*q(i,irho)


#if NDUST>0



        do idust=1,ndust


            epsilondust(i,idust) = dust2gas_ratio(idust) !Allows the user to define by hand the dust2gas of each grain species (and thus overwrites that produced by call distribution_dust(.true.))
            


            !q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho)

            q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho) !DUSTYWAVE


         

            q(i,ivdx(idust)) = vx_0/unit_v*sin(2.0d0*pi*position(i,1)/box_l*kx_wave) !Velocity perturbation for the dust


        end do


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
  namelist/setup_params/box_l,rho_0,K,dust2gas_ratio,vx_0,kx_wave,cs_0,delta_rho,delta_rho_d,i_coupled_species
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
        !tstop(i,idust) = St_0(idust)*box_l*rho_0/cs(i)/q(i,irho)
        !tstop(i,idust) = q(i,irhod(idust))*q(i,irho)/K(idust)/(q(i,irhod(idust))+q(i,irho))
        tstop(i,idust) = q(i,irhod(idust))*q(i,irho)/K(idust)/(q(i,irho))

     end do
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine compute_tstop
#endif


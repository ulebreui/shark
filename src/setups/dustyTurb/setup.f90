subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,perturbation
  real(dp) :: rmax,vol_tot,xx,yy,vkep,vx_nak,vy_nak
  real(dp) :: B_field
   real(dp), dimension(1:2*ndust+2) :: v_sol

  real (dp) :: A_nak, B_nak
  real(dp), dimension(1:ndust) :: ux_nak, uy_nak
  integer :: ind_ux,ind_uy,info


  integer :: i,idust,jdust,imax,ix,iy,icell
  call allocate_init



  tcross     = box_l/cs0/Mach

  allocate(f_hat(1:numb_modes))
  allocate(f_hat_old(1:numb_modes))
  allocate(phase_forcing(1:numb_modes))

  call distribution_dust(.true.)



  call gridinit(box_l,box_l_y)
  q = 0.0d0
  iso_cs = 1
  do i =1,ncells
      xx=position(i,1)-half*box_l  ! Boxlen already in pc
      yy=position(i,2)-half*box_l_y

        q(i,irho)  = rho_init

        q(i,ivx)   = 0.0d0
        q(i,ivy)   = 0.0d0
        q(i,iP)    = q(i,irho)*cs0**2.0
#if NDUST>0
     do idust=1,ndust

        q(i,irhod(idust))    = dust2gas_species(idust)*rho_init!+ perturbation
        epsilondust(i,idust) = dust2gas_species(idust)

        ! ts = St tcross = rhog sg/rho/cs
        ! sg = St tcross*rho*cs0/rhog

        sdust(i,idust)       = Stokes_species(idust)*rho_init*cs0/rhograin*tcross

        q(i,ivdx(idust))      = 0.0d0

        q(i,ivdy(idust))     = 0.0d0
     end do
#endif
  end do
  cs=cs0


call primtoc
call apply_boundaries

call update_force_setup

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
  namelist/setup_params/box_l,box_l_y,rho_init,Stokes_species,dust2gas_species
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
   integer :: i,idust
   real(dp) :: u, v,Ohmdt, AA, BB,d,old_ek,new_ek
   
   return

end subroutine setup_inloop


 subroutine update_force_setup
   use parameters
   use commons
   use units
   implicit none

   integer :: i,idust,k_mode,ix,iy,icell
   real(dp) :: xx,yy,intensity
   ! NOT THE RIGHT way

   
!    do k_mode = 1, numb_modes
!       call random_number(intensity)
!       f_hat(k_mode)       =  2.0d0*(intensity-1.0d0) * Mach*cs0/ turnover_time
!       call random_number(intensity)      
!       phase_forcing(k_mode)= 2.0*pi*intensity
!     end do
!     force      = 0.0d0
!     force_dust = 0.0d0
!     do iy = first_active_y,last_active_y
!         do ix = first_active,last_active
!         xx=position(icell(ix,iy),1)-half*box_l  ! Boxlen already in pc
! #if NY>1
!         yy=position(icell(ix,iy),2)-half*box_l_y
! #endif    
!         do k_mode=1,numb_modes
!             force(icell(ix,iy),1)= force(icell(ix,iy),1) + f_hat(k_mode)* cos(phase_forcing(k_mode)*xx/box_l)  
! #if NY>1
!             force(icell(ix,iy),2)= force(icell(ix,iy),2) + f_hat(k_mode) *sin(phase_forcing(k_mode)*yy/box_l)
! #endif        
!       end do
!    end do
!   end do

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
        tstop(i,idust) = rhograin*sdust(i,idust)/rho_init/cs0
     end do
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine compute_tstop
#endif


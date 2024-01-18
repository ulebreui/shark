subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud,perturbation,cs0
  real(dp) :: rmax,vol_tot,xx,yy,vkep,vx_nak,vy_nak
  real(dp) :: B_field
  real(dp), dimension(1:2*ndust+2) :: v_sol
  real (dp) :: A_nak, B_nak
  real(dp), dimension(1:ndust) :: ux_nak, uy_nak
  integer :: ind_ux,ind_uy,info


  integer :: i,idust,jdust,imax,ix,iy,icell
  call allocate_init



  eta_stream = 0.05*HoverR
  box_l      = box_l  *eta_stream*rad0
  box_l_y    = box_l_y*eta_stream*rad0

  vkep    = rad0*Omega_shear
  cs0     = (hoverR)*vkep

  ice_mantle = 0.0d0 ! No ice for SI test
  rhograin   = rho_init/theta_dust

  smin    = Stokes_min*rho_init*cs0/rhograin/omega_shear
  smax    = Stokes_max*rho_init*cs0/rhograin/omega_shear
  scut    = Stokes_cut*rho_init*cs0/rhograin/omega_shear
  scutmin = Stokes_min*rho_init*cs0/rhograin/omega_shear


  vfrag = vfrag * cs0 ! Quantify vfrag in terms of cs
  call distribution_dust(.true.)

  if(stokes_distrib) then

    do idust=1,ndust
        dust2gas_species(idust) = epsilondust(1,idust)
        stokes_species(idust)   = sdust(1,idust)/(rho_init*cs0/rhograin/omega_shear)
    end do

  endif

  A_nak = 0.0d0
  B_nak = 1.0d0

  do idust=1,ndust

    A_nak = A_nak + dust2gas_species(idust)*Stokes_species(idust)/(1.0d0+Stokes_species(idust)**2)
    B_nak = B_nak + dust2gas_species(idust)/(1.0d0+Stokes_species(idust)**2)

  end do

  vx_nak = 2.0d0*eta_stream*vkep*A_nak/(A_nak+B_nak)
  vy_nak = -eta_stream*vkep*B_nak/(A_nak+B_nak)

  do idust=1,ndust

    ux_nak(idust) = (vx_nak+2.0d0*Stokes_species(idust)*vy_nak)/(1.0d0+Stokes_species(idust)**2)
    uy_nak(idust) = (vy_nak-half*Stokes_species(idust)*vx_nak)/(1.0d0+Stokes_species(idust)**2)

  end do

  call gridinit(box_l,box_l_y)
  q = 0.0d0
  iso_cs = 1
  do i =1,ncells
        xx=position(i,1)-half*box_l  ! Boxlen already in pc
        yy=position(i,2)-half*box_l_y

        q(i,irho)  = rho_init
        call get_rhoturb(2d-2*cs0,perturbation)

        q(i,ivx)    = vx_nak + perturbation
        call get_rhoturb(2d-2*cs0,perturbation)

        q(i,ivy)   = perturbation
        call get_rhoturb(2d-2*cs0,perturbation)
        q(i,ivz)   =  q_shear*Omega_shear*xx+vy_nak + perturbation
        q(i,iP)    = q(i,irho)*cs0**2.0
#if NDUST>0
     do idust=1,ndust

        q(i,irhod(idust))    = dust2gas_species(idust)*rho_init!+ perturbation
        epsilondust(i,idust) = dust2gas_species(idust)

        if(.not. stokes_distrib) sdust(i,idust)       = Stokes_species(idust)*rho_init*cs0/rhograin/omega_shear
        call get_rhoturb(2d-2*cs0,perturbation)

        q(i,ivdx(idust))      = perturbation+ux_nak(idust)
        call get_rhoturb(2d-2*cs0,perturbation)

        q(i,ivdy(idust))     =  perturbation
        call get_rhoturb(2d-2*cs0,perturbation)
        q(i,ivdz(idust))     = perturbation+ q_shear*Omega_shear*xx + uy_nak(idust)
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
  namelist/setup_params/box_l,box_l_y,rho_init,omega_shear,q_shear,eta_stream,HoverR,Stokes_species,dust2gas_species, theta_dust,Stokes_min,  Stokes_max,  Stokes_cut,stokes_distrib 
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

    do i=1,ncells
    if(active_cell(i))then
    !Crank nicholson scheme
        u = u_prim(i,ivz)/u_prim(i,irho) 
        v = u_prim(i,ivx)/u_prim(i,irho) 
        d = u_prim(i,irho) 
        Ohmdt = Omega_shear*dt
        AA    = 2.0d0*q_shear*Omega_shear*(position(i,1)-half*box_l)-2.0d0*rad0*eta_stream 
        u_prim(i,ivz) = d*( u*(1.-Ohmdt**2) + 2.*v*Ohmdt + AA*Ohmdt**2 ) / (1.+Ohmdt**2) 
        u_prim(i,ivx)  = d*( v*(1.-Ohmdt**2) - 2.*u*Ohmdt + AA*Ohmdt ) / (1.+Ohmdt**2) 
#if NDUST>0
        Ohmdt=Omega_shear*dt
        AA = 2.0d0*q_shear*Omega_shear*(position(i,1)-half*box_l)
        do idust=1,ndust
          u = u_prim(i,ivdz(idust))/u_prim(i,irhod(idust)) 
          v = u_prim(i,ivdx(idust))/u_prim(i,irhod(idust)) 
          d = u_prim(i,irhod(idust)) 
          u_prim(i,ivdz(idust)) = d*( u*(1.-Ohmdt**2) + 2.*v*Ohmdt + AA*Ohmdt**2 ) / (1.+Ohmdt**2) 
          u_prim(i,ivdx(idust))  = d*( v*(1.-Ohmdt**2) - 2.*u*Ohmdt + AA*Ohmdt ) / (1.+Ohmdt**2) 
        end do
#endif
    endif
   end do
end subroutine setup_inloop


 subroutine update_force_setup
   use parameters
   use commons
   use units
   implicit none

   integer :: i,idust

   !return
! #if NDUST>0
!     call compute_tstop
! #endif


   do i=1,ncells
    if(active_cell(i))then
        force(i,1)  = 2.0d0*q_shear*Omega_shear**2.*(position(i,1)-half*box_l) -2.0d0*Omega_shear*q(i,ivz) -  2.0d0*rad0*Omega_shear*eta_stream
        force(i,2)  = 0.0d0
        force(i,3)  = 2.0d0*Omega_shear*q(i,ivx)
#if NDUST>0
        do idust=1,ndust
         force_dust(i,1,idust)  = 2.0d0*q_shear*Omega_shear**2.*(position(i,1)-half*box_l) -2.0d0*Omega_shear*q(i,ivdz(idust))
         force_dust(i,2,idust)  = 0.0d0
         force_dust(i,3,idust)  = 2.0d0*Omega_shear*q(i,ivdx(idust))

         if(.not. drag) then
          force(i,1)  = force(i,1)  - q(i,irhod(idust))/q(i,irho)*(q(i,ivx) - q(i,ivdx(idust)))/tstop(i,idust)
          force(i,2)  = force(i,2)  - q(i,irhod(idust))/q(i,irho)*(q(i,ivy) - q(i,ivdy(idust)))/tstop(i,idust)
          force(i,3)  = force(i,3)  - q(i,irhod(idust))/q(i,irho)*(q(i,ivz) - q(i,ivdz(idust)))/tstop(i,idust)
          force_dust(i,1,idust)  = force_dust(i,1,idust)  + (q(i,ivx) - q(i,ivdx(idust)))/tstop(i,idust)
          force_dust(i,2,idust)  = force_dust(i,2,idust)  + (q(i,ivy) - q(i,ivdy(idust)))/tstop(i,idust)
          force_dust(i,3,idust)  = force_dust(i,3,idust)  + (q(i,ivz) - q(i,ivdz(idust)))/tstop(i,idust)
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
  real(dp):: pn,rhon
  !Re-calc distribution


  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do i=1,ncells
   if(active_cell(i)) then
     do idust=1,ndust
        tstop(i,idust) = rhograin*sdust(i,idust)/rho_init/(rad0*Omega_shear*hoverr)
     end do
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine compute_tstop
#endif


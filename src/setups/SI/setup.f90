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

  !real(dp) :: ux_nak, uy_nak

  integer :: ind_ux,ind_uy,info


  integer :: i, idust,jdust,imax,ix,iy,icell,ixx,iyy
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
  call distribution_dust

  if(stokes_distrib) then

    do idust=1,ndust
        dust2gas_species(idust) = epsilondust(1,idust)
        stokes_species(idust)   = sdust(idust)/(rho_init*cs0/rhograin/omega_shear)
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
        ix=ixx(i)
        iy=iyy(i)
        xx=position(ix,iy,1)-half*box_l  ! Boxlen already in pc
        yy=position(ix,iy,2)-half*box_l_y

        q(irho,ix,iy)  = rho_init
        call get_rhoturb(2d-2*cs0,perturbation)

        q(ivx,ix,iy)    = vx_nak + perturbation
        call get_rhoturb(2d-2*cs0,perturbation)

        q(ivy,ix,iy)   = perturbation
        call get_rhoturb(2d-2*cs0,perturbation)
        q(ivz,ix,iy)   =  q_shear*Omega_shear*xx+vy_nak + perturbation
        q(iP,ix,iy)    =  q(irho,ix,iy)*cs0**2.0
#if NDUST>0
     do idust=1,ndust

        q(irhod(idust),ix,iy)    = dust2gas_species(idust)*rho_init!+ perturbation
        epsilondust(i,idust)     = dust2gas_species(idust)

        if(.not. stokes_distrib) sdust(idust)       = Stokes_species(idust)*rho_init*cs0/rhograin/omega_shear
        call get_rhoturb(2d-2*cs0,perturbation)

        q(ivdx(idust),ix,iy)      = perturbation+ux_nak(idust)
        call get_rhoturb(2d-2*cs0,perturbation)

        q(ivdy(idust),ix,iy)     =  perturbation
        call get_rhoturb(2d-2*cs0,perturbation)
        q(ivdz(idust),ix,iy)     = perturbation+ q_shear*Omega_shear*xx + uy_nak(idust)
     end do
#endif
  end do
  cs=cs0


call primtoc
call apply_boundaries
call compute_tstop
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
   integer :: icount,iout,i,idust,iy,ix
   real(dp):: mtot,mdtot
   logical :: outputing,verbose

        !We make an output at a certain frequency rate of for specific values of the density. This can be tuned at will
     if(icount.eq.freq_out) then
        print *, "time =",time,' tend = ', tend
        verbose=.true.
        outputing=.true.
        icount=0
     endif

     if(time.eq.0) outputing = .true.   
     if(outputing) call output(iout)
     if(outputing) iout=iout+1
     if(outputing) print *, "Outputing data "

     if(outputing) then
        mtot=0.0d0
        mdtot=0.0d0
        do iy = first_active_y,last_active_y
          do ix = first_active,last_active
          mtot=mtot+u_prim(irho,ix,iy)
          do idust=1,ndust
            mdtot=mdtot+u_prim(irhod(idust),ix,iy)
          end do
          end do
        end do
     end if

     if(outputing) print *, " Total mass is", mtot
     if(outputing) print *, " Total dust mass is", mdtot

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
   integer :: idust,ix,iy
   real(dp) :: Ohmdt
   real(dp) ::  u, v, AA, BB,d
   
   if(force_kick) return

   Ohmdt = Omega_shear*dt
   !$omp parallel do default(shared) private(idust, ix,iy,u, v, Ohmdt, AA, BB,d)
    do iy = first_active_y,last_active_y
      do ix = first_active,last_active
          !Crank nicholson scheme
          u = u_prim(ivz,ix,iy)/u_prim(irho,ix,iy) 
          v = u_prim(ivx,ix,iy)/u_prim(irho,ix,iy) 
          d = u_prim(irho,ix,iy) 
          AA    = 2.0d0*q_shear*Omega_shear*(position(ix,iy,1)-half*box_l)-2.0d0*rad0*eta_stream 
          u_prim(ivz,ix,iy)  = d*( u*(1.-Ohmdt**2) + 2.*v*Ohmdt + AA*Ohmdt**2 ) / (1.+Ohmdt**2) 
          u_prim(ivx,ix,iy)  = d*( v*(1.-Ohmdt**2) - 2.*u*Ohmdt + AA*Ohmdt ) / (1.+Ohmdt**2) 
#if NDUST>0
          AA    = 2.0d0*q_shear*Omega_shear*(position(ix,iy,1)-half*box_l)
          do idust=1,ndust
            u = u_prim(ivdz(idust),ix,iy)/u_prim(irhod(idust),ix,iy) 
            v = u_prim(ivdx(idust),ix,iy)/u_prim(irhod(idust),ix,iy) 
            d = u_prim(irhod(idust),ix,iy) 
            u_prim(ivdz(idust),ix,iy)  = d*( u*(1.-Ohmdt**2) + 2.*v*Ohmdt + AA*Ohmdt**2 ) / (1.+Ohmdt**2) 
            u_prim(ivdx(idust),ix,iy)  = d*( v*(1.-Ohmdt**2) - 2.*u*Ohmdt + AA*Ohmdt )    / (1.+Ohmdt**2) 
          end do
#endif
     end do
   end do
end subroutine setup_inloop


 subroutine update_force_setup
   use parameters
   use commons
   use units
   use OMP_LIB
   implicit none

   integer ::idust,ix,iy

   if(.not.force_kick) return
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy)
    do iy = first_active_y,last_active_y
      do ix = first_active,last_active
        force_x(ix,iy)  = 2.0d0*q_shear*Omega_shear**2.*(position(ix,iy,1)-half*box_l) -2.0d0*Omega_shear*q(ivz,ix,iy) -  2.0d0*rad0*Omega_shear*eta_stream
        force_y(ix,iy)  = 0.0d0
        force_z(ix,iy)  = 2.0d0*Omega_shear*q(ivx,ix,iy)
#if NDUST>0
        do idust=1,ndust
         force_dust_x(idust,ix,iy)  = 2.0d0*q_shear*Omega_shear**2.*(position(ix,iy,1)-half*box_l) -2.0d0*Omega_shear*q(ivdz(idust),ix,iy)
         force_dust_y(idust,ix,iy)  = 0.0d0
         force_dust_z(idust,ix,iy)  = 2.0d0*Omega_shear*q(ivdx(idust),ix,iy)

         if(.not. drag) then
          force_x(ix,iy)  = force_x(ix,iy)  - q(irhod(idust),ix,iy)/q(irho,ix,iy)*(q(ivx,ix,iy) - q(ivdx(idust),ix,iy))/tstop(idust,ix,iy)
          force_y(ix,iy)  = force_y(ix,iy)  - q(irhod(idust),ix,iy)/q(irho,ix,iy)*(q(ivy,ix,iy) - q(ivdy(idust),ix,iy))/tstop(idust,ix,iy)
          force_z(ix,iy)  = force_z(ix,iy)  - q(irhod(idust),ix,iy)/q(irho,ix,iy)*(q(ivz,ix,iy) - q(ivdz(idust),ix,iy))/tstop(idust,ix,iy)
          force_dust_x(idust,ix,iy)  = force_dust_x(idust,ix,iy)  + (q(ivx,ix,iy) - q(ivdx(idust),ix,iy))/tstop(idust,ix,iy)
          force_dust_y(idust,ix,iy)  = force_dust_y(idust,ix,iy)  + (q(ivy,ix,iy) - q(ivdy(idust),ix,iy))/tstop(idust,ix,iy)
          force_dust_z(idust,ix,iy)  = force_dust_z(idust,ix,iy)  + (q(ivz,ix,iy) - q(ivdz(idust),ix,iy))/tstop(idust,ix,iy)
         endif
        end do
#endif
    end do
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
  real(dp):: pn,rhon
  !Re-calc distribution


   integer :: idust,ix,iy
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy)
   do iy = first_active_y,last_active_y
    do ix = first_active,last_active
     do idust=1,ndust
        tstop(idust,ix,iy) = rhograin*sdust(idust)/rho_init/(rad0*Omega_shear*hoverr)
      end do
     end do
  end do


end subroutine compute_tstop
#endif


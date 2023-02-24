subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud,perturbation,cs0
  real(dp) :: rmax,vol_tot,xx,yy
  real(dp) :: B_field

  integer :: i,idust,imax,ix,iy,icell

  print *, 'tend = ', tend
  call allocate_init
  box_l   = box_l*eta_stream*rad0
  box_l_y = box_l_y*eta_stream*rad0
  cs0     = ((hoverR)*rad0*Omega_shear)
  call gridinit(box_l,box_l_y)
  q=0.0d0
  do iy = 1,ny_max
    do ix = 1,nx_max
      xx=position(icell(ix,iy),1)-half*box_l  ! Boxlen already in pc
      yy=position(icell(ix,iy),2)-half*box_l_y

      q(icell(ix,iy),irho)  = rho_init*(1.0d0+ 1d-2*cos(2.*pi*6*xx/box_l+2.*pi*6*yy/box_l_y)) !+ perturbation
      call get_rhoturb(1d-2*cs0,perturbation)

      q(icell(ix,iy),iv)    = perturbation
      call get_rhoturb(1d-2*cs0,perturbation)

      q(icell(ix,iy),ivy)   = perturbation
#if IVZ==1
      call get_rhoturb(1d-2*cs0,perturbation)
      q(icell(ix,iy),ivz)   = q_shear*Omega_shear*xx + perturbation
#endif
      q(icell(ix,iy),iP)    = q(icell(ix,iy),irho)*cs0**2.0
  enddo
  enddo
#if NDUST>0
  call distribution_dust(.true.)
   do iy = 1,ny_max
    do ix = 1,nx_max
     do idust=1,ndust
        q(icell(ix,iy),irhod(idust))   = dust2gas*rho_init
        epsilondust(icell(ix,iy),idust)= dust2gas
        sdust(icell(ix,iy),idust)      = scut*rho_init*cs0/(rhograin/unit_d)/sqrt(pi*gamma/8.0d0)/omega_shear ! ts =0.1 at rho0
        q(icell(ix,iy),ivd(idust))     = 0.d0
        q(icell(ix,iy),ivdy(idust))    = 0.d0
        xx=position(icell(ix,iy),1)-half*box_l  ! Boxlen already in pc
#if IVZ==1        
        q(icell(ix,iy),ivdz(idust))   = q_shear*Omega_shear*xx!omega_shear*xx
#endif
     end do
  end do
  end do
#endif

  call primtoc
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
  namelist/setup_params/box_l,box_l_y,rho_init,omega_shear,q_shear,eta_stream,HoverR
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
     if(outputing) print *, "Total momentum is", sum(u_prim(:,iv)+u_prim(:,ivy)+u_prim(:,ivz))
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
   !return
   do i=1,ncells
    if(active_cell(i)==1)then
    !Crank nicholson scheme
        u = q(i,ivz)
        v = q(i,iv)
        d = q(i,irho) 
        old_ek= half*d*(v**2+u**2)
        Ohmdt = Omega_shear*dt
        AA = 2.0d0*q_shear*Omega_shear*(position(i,1)-half*box_l)-2.0d0*rad0*Omega_shear*eta_stream 
        u_prim(i,ivz) = d*( u*(1.-Ohmdt**2) + 2.*v*Ohmdt + AA*Ohmdt**2 ) / (1.+Ohmdt**2) 
        u_prim(i,iv)  = d*( v*(1.-Ohmdt**2) - 2.*u*Ohmdt + AA*Ohmdt ) / (1.+Ohmdt**2) 
        new_ek=half*(u_prim(i,ivz)**2+u_prim(i,iv)**2)/d
        u_prim(i,iP) = u_prim(i,iP) + new_ek-old_ek
#if NDUST>0
        Ohmdt=Omega_shear*dt
        AA = 2.0d0*q_shear*Omega_shear*(position(i,1)-half*box_l)
        do idust=1,ndust
          u = q(i,ivdz(idust))
          v = q(i,ivd(idust))
          d = q(i,irhod(idust)) 
          u_prim(i,ivdz(idust)) = d*( u*(1.-Ohmdt**2) + 2.*v*Ohmdt + AA*Ohmdt**2 ) / (1.+Ohmdt**2) 
          u_prim(i,ivd(idust))  = d*( v*(1.-Ohmdt**2) - 2.*u*Ohmdt + AA*Ohmdt ) / (1.+Ohmdt**2) 
        end do
#endif
    endif
   end do
  call ctoprim
end subroutine setup_inloop


 subroutine update_force_setup
   use parameters
   use commons
   use units
   implicit none

   integer :: i,idust


   return ! remove return for explicit force and add it to setup_inloop



   do i=1,ncells
    if(active_cell(i)==1)then
        force(i,1)  =  2.0d0*q_shear*Omega_shear**2.*(position(i,1)-half*box_l) -2.0d0*Omega_shear*q(i,ivz)
        force(i,1)  = force(i,1)-2.0d0*rad0*Omega_shear**2*eta_stream 
        force(i,2)  = 0.0d0
        force(i,3)  = 2.0d0*Omega_shear*q(i,iv)
#if NDUST>0
        do idust=1,ndust
         force_dust(i,1,idust)  = 2.0d0*q_shear*Omega_shear**2.*(position(i,1)-half*box_l) -2.0d0*Omega_shear*q(i,ivdz(idust))
         force_dust(i,2,idust)  = 0.0d0
         force_dust(i,3,idust)  = 2.0d0*Omega_shear*q(i,ivd(idust))
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
   if(active_cell(i)==1) then
     do idust=1,ndust
        tstop(i,idust) = sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(i,idust)/(q(i,irho)*cs(i))
     end do
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine compute_tstop
#endif

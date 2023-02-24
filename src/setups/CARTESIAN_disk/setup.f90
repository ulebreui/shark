subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud,perturbation,cs0,rr
  real(dp) :: rmax,vol_tot,xx,yy,Omega,H,cutoff
  real(dp) :: B_field

  integer :: i,idust,imax,ix,iy,icell

  call allocate_init

  ! Dimensionless units
  sigma_R0   = sigma_R0/unit_dcol
  R0_disk    = R0_disk*au/unit_l
  R_out_disk = R_out_disk*au/unit_l
  Rplanet    = Rplanet*au/unit_l

  smooth_r   = smooth_r*au/unit_l
  Mstar      = Mstar *Msun/unit_m
  Mplanet    = Mplanet *Msun/unit_m

  box_l   = box_l   * R_out_disk
  box_l_y = box_l_y * R_out_disk

  tend = tend /(sqrt(Mstar/R0_disk**3.0))

  print *, 'tend = ', tend, '1 au orbits'

  call gridinit(box_l,box_l_y)

  q=0.0d0
  iso_cs =1
  do iy = 1,ny_max
    do ix = 1,nx_max
      xx    = position(icell(ix,iy),1)-half*box_l  ! Boxlen already in pc
      yy    = position(icell(ix,iy),2)-half*box_l_y
      rr    = sqrt(xx**2+yy**2+(smooth_r*dx(icell(ix,iy),1))**2)
      H     = HoverR*rr
      Omega = sqrt(Mstar/rr**3) ! /!\ Units must have G = 1
      cs0   = H * Omega
      cutoff= 1.0d0
      if(rr>R_out_disk) cutoff = 1e-3
      q(icell(ix,iy),irho)  = sigma_R0*(R0_disk/rr)*cutoff
      q(icell(ix,iy),iv)    = -omega*yy
      q(icell(ix,iy),ivy)   = omega*xx
      q(icell(ix,iy),iP)    = q(icell(ix,iy),irho)*cs0**2.0
      cs(icell(ix,iy))      = cs0
#if NDUST>0
     do idust=1,ndust
        q(icell(ix,iy),irhod(idust))   = dust2gas*q(icell(ix,iy),irho)
        epsilondust(icell(ix,iy),idust)= dust2gas
        sdust(icell(ix,iy),idust)      = scut/unit_l ! The 2H term comes from the vertical integration 
        q(icell(ix,iy),ivd(idust))     = -omega*yy
        q(icell(ix,iy),ivdy(idust))    = omega*xx
     end do
#endif
  enddo
  enddo
#if NDUST>0  
  call distribution_dust(.true.)
#endif

  call primtoc
  call update_force_setup
#if NDUST>0  
  call compute_tstop
#endif
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
  namelist/setup_params/buffer_size,box_l,box_l_y, sigma_R0, R0_disk,Mstar , R_out_disk  , HoverR ,smooth_r,Mstar,Mplanet

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
     if(outputing) print *, "Total momentum is", sum(u_prim(:,iv)+u_prim(:,ivy))
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
   use OMP_LIB

   implicit none
   integer :: i,icell,idust
   real(dp) :: rr, H, Omega,cs0,xx,yy,Trel,rho_ana,v_ana,vy_ana,P_ana,rate,cutoff

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE( i,icell,idust, rr, H, Omega,cs0,xx,yy,Trel,rho_ana,v_ana,vy_ana,P_ana,rate,cutoff)
  !$OMP DO
!We reinitialise at the center
  do i= 1,ncells
    if(active_cell(i)==1) then
      xx    = position(i,1)-half*box_l  ! Boxlen already in pc
      yy    = position(i,2)-half*box_l_y
      rr    = sqrt(xx**2+yy**2+(smooth_r*dx(i,1))**2)

      if(rr<R_out_disk*buffer_size .or. rr > R_out_disk) then
        H     = HoverR*rr
        Omega = sqrt(Mstar/rr**3) ! /!\ Units must have G = 1
        !Relaxation timescale
        Trel  = 1./Omega
        rate  = 0.0d0
        if(rr<R0_disk .or. rr > R_out_disk) rate=0.0d0
        cs0   = H * Omega
        cutoff=1.0d0
        if(rr>R_out_disk) cutoff = 1e-3

        rho_ana = sigma_R0*(R0_disk/rr)*cutoff
        v_ana   = -omega*rho_ana

        u_prim(i,irho)  = rho_ana+(u_prim(i,irho)-rho_ana)*rate
        u_prim(i,iv)    =  -yy*v_ana +(u_prim(i,iv)+yy*v_ana)*rate
        u_prim(i,ivy)   =  xx*v_ana +(u_prim(i,ivy)-xx*v_ana)*rate
        u_prim(i,iP)    = (rho_ana+(u_prim(i,irho)-rho_ana)*rate)*cs0**2.0
#if NDUST>0
        do idust=1,ndust
            u_prim(i,irhod(idust))   = dust2gas*rho_ana+(u_prim(i,irhod(idust))-rho_ana*dust2gas)*rate
            u_prim(i,ivd(idust))     = -yy*v_ana*dust2gas +(u_prim(i,ivd(idust))+dust2gas*yy*v_ana)*rate
            u_prim(i,ivdy(idust))    = xx*v_ana*dust2gas +(u_prim(i,ivdy(idust))-dust2gas*xx*v_ana)*rate
        end do
#endif
    endif  
    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
call ctoprim
end subroutine setup_inloop


 subroutine update_force_setup
   use parameters
   use commons
   use units
   use OMP_LIB

   implicit none

   integer :: i,idust
   real(dp) :: xx, yy, rr, xx_p, yy_p, rr_p, theta_p
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,xx,yy,rr, xx_p, yy_p, rr_p,theta_p)
  !$OMP DO
   do i=1,ncells
    if(active_cell(i)==1)then

        rr = sqrt(xx**2+yy**2+(smooth_r*dx(i,1))**2)
        force(i,1)  = -Mstar*xx/rr**3
        force(i,2)  = -Mstar*yy/rr**3
#if NDUST>0
        do idust=1,ndust
         force_dust(i,1,idust)  = -Mstar*xx/rr**3
         force_dust(i,2,idust)  = -Mstar*yy/rr**3
        end do
#endif
    if(Mplanet>0.0d0) then
        theta_p= 2.0d0*pi*mod(time/(2.0d0*pi/sqrt(Mstar/Rplanet**3)),1.0d0) ! Orbital period
        xx_p = cos(theta_p)*Rplanet
        yy_p = sin(theta_p)*Rplanet
        rr_p = sqrt(xx_p**2+yy_p**2+(smooth_r*dx(i,1))**2)

        force(i,1)  = force(i,1) -Mplanet*xx_p/rr_p**3
        force(i,2)  = force(i,1)-Mplanet*yy_p/rr_p**3
#if NDUST>0
        do idust=1,ndust
         force_dust(i,1,idust)  = force_dust(i,1,idust) -Mplanet*xx_p/rr_p**3
         force_dust(i,2,idust)  = force_dust(i,2,idust) -Mplanet*yy_p/rr_p**3
        end do
#endif
    endif
    endif

   end do
  !$OMP END DO  
  !$OMP END PARALLEL  
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
  real(dp):: xx, yy, H,rr
  !Re-calc distribution


  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,xx,yy,H)
  !$OMP DO
  do i=1,ncells
   if(active_cell(i)==1) then
        xx=position(i,1)-half*box_l  ! Boxlen already in pc
        yy=position(i,2)-half*box_l_y
        rr    = sqrt(xx**2+yy**2+smooth_r**2)
        H     = HoverR*rr
        do idust=1,ndust
            tstop(i,idust) = sqrt(2.0d0*pi*H)*sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(i,idust)/(q(i,irho)*cs(i))
        end do
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine compute_tstop
#endif


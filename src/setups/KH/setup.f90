subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud,perturbation,cs_eos
  real(dp) :: rmax,vol_tot,xx,yy
  real(dp) :: B_field

  integer :: i,idust,imax,ix,iy,icell

  tend = tend * 3.15576e13/unit_t ! Convert tend from Myr to time unit
  print *, 'tend = ', tend
  call allocate_init
  call gridinit(box_l,box_l_y)
  q=0.0d0
  do iy = 1,ny_max
    do ix = 1,nx_max
      xx=position(icell(ix,iy),1)-half*box_l  ! Boxlen already in pc
      yy=position(icell(ix,iy),2)-half*box_l_y

      q(icell(ix,iy),irho)  = rho_dense/unit_d
      q(icell(ix,iy),iv)    = cs_eos(T_dense)*Mach/unit_v
      !call get_vxturb(vy0*cs_eos(T_dense)*Mach,perturbation)
      q(icell(ix,iy),ivy)   = (vy0*cs_eos(T_dense)*Mach*cos(xx*2.0d0*kx*acos(-1.0d0)/box_l))/unit_v*half
      q(icell(ix,iy),iP)    = rho_dense*cs_eos(T_dense)**2/unit_P  
      if(((yy/box_l_y)<0.d0)) then
          q(icell(ix,iy),irho) = rho_diffuse/unit_d
          q(icell(ix,iy),iv)   = -(vy0*cs_eos(T_dense)*Mach*cos(xx*2.0d0*kx*acos(-1.0d0)/box_l))/unit_v*half
          q(icell(ix,iy),iP)   = rho_dense*cs_eos(T_dense)**2/unit_P !rho_diffuse*cs_eos(T_diffuse)**2/unit_P 
      endif
    end do
  end do
#if NDUST>0
  call distribution_dust(.true.)
  do i=1,ncells
     do idust=1,ndust
        if(ndust==1)  epsilondust(i,idust) = 0.01 
        q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho)
        if(ndust==1)  sdust(i,idust) = scut/unit_l ! 1 micron
     end do
  end do
#endif
  call primtoc
end subroutine setup



subroutine get_vxturb(pert,delvx)

  use random
  use precision
  implicit none

  integer :: i
  integer::iseed=0         
  integer ,dimension(1,1:IRandNumSize)    :: allseed
  integer,dimension(IRandNumSize) :: localseed=-1
  real(dp) ::  pert,delvx,randno

  if (localseed(1)==-1) then
     call rans(1,iseed,allseed)
     localseed = allseed(1,1:IRandNumSize)
  end if

  call ranf(localseed,randno)

  delvx = randno*pert


end subroutine get_vxturb

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
  namelist/setup_params/box_l,box_l_y,rho_dense,rho_diffuse,P_dense,P_diffuse,v_dense,v_diffuse,vy0,kx,Mach,T_dense,T_diffuse
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
     if(outputing) print *, "Total mass is", sum(uold(:,irho))
     if(outputing) print *, "Total momentum is", sum(uold(:,iv))
     if(outputing) print *, "Total energy is", sum(uold(:,iP))

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

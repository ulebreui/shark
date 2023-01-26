subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud,perturbation
  real(dp) :: rmax,vol_tot,xx,yy
  real(dp) :: B_field

  integer :: i,idust,imax,ix,iy,icell


  call allocate_init
  call gridinit(box_l)
  q=0.0d0
  do iy = 1,ny_max
    do ix = 1,nx_max
      !print*, icell(ix,iy)
      !print*, active_cell(icell(ix,iy)),ix,iy
      xx=position(icell(ix,iy),1)-half*box_l
      yy=position(icell(ix,iy),2)-half*box_l
      q(icell(ix,iy),irho)  = rho_dense
      q(icell(ix,iy),iv)    = v_dense
      call get_vxturb(vy0,perturbation)
      !print *, perturbation
      q(icell(ix,iy),ivy)   = vy0*cos(xx*2.0d0*kx*acos(-1.0d0)/box_l) + perturbation
      q(icell(ix,iy),iP)    = P_dense   
      if((abs(yy)>box_l/4.0d0)) then
          q(icell(ix,iy),irho) = rho_diffuse
          q(icell(ix,iy),iv)   = v_diffuse
          q(icell(ix,iy),iP)   = P_diffuse
      endif
    end do
  end do
#if NDUST>0
  call distribution_dust(.true.)
  do i=1,ncells
     do idust=1,ndust
        epsilondust(i,idust) = 0.01 
        q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho)
        sdust(i,idust) = 1.0d0
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
  namelist/setup_params/box_l,rho_dense,rho_diffuse,P_dense,P_diffuse,v_dense,v_diffuse,vy0,kx
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

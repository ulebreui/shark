subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud
  real(dp) :: rmax,vol_tot
  real(dp) :: B_field

  integer :: i,idust,imax,ix,iy,icell


  call allocate_init
  call gridinit(box_l,box_l_y)
  q=0.0d0
  do iy = 1,ny_max
    do ix = 1,nx_max
      !print*, icell(ix,iy)
      !print*, active_cell(icell(ix,iy)),ix,iy
      q(icell(ix,iy),irho) = rho_l
      q(icell(ix,iy),ivx)   = v_l
      q(icell(ix,iy),iP)   = P_l
      ! if (ix-first_active+1==nx_max/2) then
      !     q(icell(ix,iy),irho) = half*(rho_l+rho_r)
      !     q(icell(ix,iy),ivx)   = half*(v_l+v_r)
      !     q(icell(ix,iy),iP)   = half*(P_l+P_r)
      ! endif
      if(direction_shock==0) then
      if(position(icell(ix,iy),1)>half) then
          q(icell(ix,iy),irho) = rho_r
          q(icell(ix,iy),ivx)   = v_r
          q(icell(ix,iy),iP)   = P_r
      endif
      else
      if(position(icell(ix,iy),2)>half) then
          q(icell(ix,iy),irho) = rho_r
          q(icell(ix,iy),ivx)   = v_r
          q(icell(ix,iy),iP)   = P_r
      endif
      endif
    end do
  end do
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
  namelist/setup_params/box_l,box_l_y,rho_l,rho_r,P_l,P_r,v_l,v_r,direction_shock
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
        tstop(i,idust) = sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(i,idust)/(q(i,irho)*cs(i))
     end do
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine compute_tstop
#endif


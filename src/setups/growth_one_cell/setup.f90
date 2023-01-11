subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud
  real(dp) :: rmax

  integer :: i,idust
  
  call allocate_init
  if(kernel_type==1)list_times=(/0.0d0,10.0d0,100.0d0,1000.0d0,1d4,1d5,1d6/)
  !if(kernel_type==2)list_times=(/0,2,4,6,8,8,8/)
  if(kernel_type==2)list_times=(/0,1,2,3,3,3,3/)
  if(kernel_type==1)list_times=(/0.0d0,1d4,2d4,3d4,3d4,3d4,3d4/)
  
  
  ilist=2
  tend=maxval(list_times)
 
  !initialisation of the hydro variables
  q=0
  uold=0
  unew=0
#if NDUST>0
  call distribution_dust(.true.)
  do i=1,ncells
        q(i,irho)= 1.0d0
        q(i,iv)= 0.0d0
     do idust=1,ndust
        q(i,irhod(idust))= epsilondust(i,idust)
     end do
  end do
#endif
  call primtoc
end subroutine setup



subroutine write_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  print *,'time', time 
  write(ilun,*) time
end subroutine write_setup_info

subroutine read_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun

end subroutine read_setup_info

subroutine read_setup_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/setup_params/list_times
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "Setup namelist reading  !"
   read(13,setup_params,IOSTAT=io)
   rewind(13)
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

   dt=(list_times(ilist)-time)
   
   ilist=ilist+1
   if(time>=tend) then
      dt=0.0
      continue_sim=.false.
   endif
 end subroutine flag_continue

 subroutine check_output(icount,iout,outputing,verbose)
   use parameters
   use commons
   use units
   implicit none
   integer :: icount,iout
   logical :: outputing,verbose

   
   outputing=.true.
   verbose=.true.
   call output(iout)
   iout=iout+1

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
end subroutine setup_inloop

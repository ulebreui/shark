subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud
  real(dp) :: rmax

  integer :: i,idust
  
  call allocate_init

  if(dt_strat.eq.1) then
     tend=maxval(list_times)
     list_times=(list_times)*(365.25*3600*24*1000.)/unit_t
  endif
  !initialisation of the hydro variables
  q=0
  uold=0
  unew=0
#if NDUST>0
  call distribution_dust(.true.)
  do i=1,ncells
        q(i,irho)= nh_cell/unit_nh
        q(i,iv)= 0.0d0
     do idust=1,ndust
        q(i,irhod(idust))= q(i,irho)*epsilondust(i,idust)
     end do
  end do
#endif
  call primtoc
  tend=tend*(365.25*3600*24*1000.)/unit_t
  if(charging) call charge
  print *, tend
end subroutine setup



subroutine write_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  write(ilun,*) time*unit_t/(365.*24.*3600.*1000.)   
end subroutine write_setup_info

subroutine read_setup_params(ilun,nmlfile)
  use parameters
  use commons
  use units

  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/setup_params/tend,nh_cell,T0_cloud,list_times,dt_fix,dt_strat,isotherm
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "Setup namelist reading  !"
   read(13,setup_params,IOSTAT=io)
   rewind(13)
   print *, 'Final time', tend
   
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"

   
 end subroutine read_setup_params
 subroutine flag_continue(continue_sim)
   use parameters
   use commons
   use units
   implicit none
   logical :: continue_sim
   !print*,'toto'

   !Here you can add flags to kill the simulation
   if(dt_strat.eq.0) then
      dt=tend
   else if (dt_strat.eq.1) then
      if(list_times(ilist+1)>tend) then
         list_times(ilist+1)=tend
      else if(list_times(ilist+1).eq.0.0) then
         list_times(ilist+1)=tend
      endif
      dt=(list_times(ilist+1)-list_times(ilist))!*(365.25*3600*24*10000.)/unit_t

      ilist=ilist+1
   else
      if(dt_fix*(365.25*3600*24*1000.)/unit_t>tend) then
         print *, tend,dt_fix*(365.25*3600*24*1000.)/unit_t
         print *, ' Warning dt > tend'
         stop
      endif
      dt=dt_fix*(365.25*3600*24*1000.)/unit_t
   endif
   !print*,'toto2'

   if(time>=tend) then
      dt=0.0
      print* , 'toto', tend
      !call output(2)
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
   !iout=2
   iout=iout+1

 end subroutine check_output

 subroutine setup_preloop
   use parameters
   use commons
   use units
   implicit none
   !call output(1)
end subroutine setup_preloop

 subroutine setup_inloop
   use parameters
   use commons
   use units
   implicit none
   !call output(1)
end subroutine setup_inloop

subroutine read_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  real(dp):: info1
  
  read(ilun,*) info1
  time=info1*(365.*24.*3600.*1000.)/unit_t
  close(ilun)
end subroutine read_setup_info

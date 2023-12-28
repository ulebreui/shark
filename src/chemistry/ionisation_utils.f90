subroutine read_ionisation_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/ionisation_params/charging,charging_all_the_time,B_0_lee,B_threshold
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "ionisation_params namelist reading  !"
   read(13,ionisation_params,IOSTAT=io)
   rewind(13)
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   if (io/=0) then
      write(*,*) 'Invalid line in ionisation namelist'
      stop
   end if
   
 end subroutine read_ionisation_params
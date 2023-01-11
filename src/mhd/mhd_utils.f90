subroutine read_mhd_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/mhd_params/B_0_lee,B_threshold
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "MHD namelist reading  !"
   read(13,mhd_params,IOSTAT=io)
   rewind(13)
   print *, 'B_0', B_0_lee
   print *, 'B_threshold', B_threshold
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   if (io/=0) then
      write(*,*) 'Invalid line in MHD namelist'
      stop
   end if
   
 end subroutine read_mhd_params

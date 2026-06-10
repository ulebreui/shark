subroutine read_ionisation_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/ionisation_params/charging,charging_all_the_time,B_0_lee,B_threshold,electrons,ions,dusty_nonideal_MHD,x,analytical_charging,&
  & analytical_charging_Wurster,analytical_charging_Shu,dusty_nonideal_MHD_no_electron,hyper_diffusion,apply_Lorentz_force_explicit,apply_Lorentz_force_implicit,only_Hall_effect,&
  & ni_coeff,call_electric_field,Hall_effect,hyper_diffusion_with_electrons, write_electric_field, write_Hall_factors,ideal_MHD
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





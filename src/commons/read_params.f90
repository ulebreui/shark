! Reads the params for the simulation in the .nml file specified when exectuting shark (e.g. ./shark my_super_setup.nml)

subroutine read_params
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile,infile
  integer :: io
  logical::nml_ok
  namelist/grid_params/tend,rin,freq_out,fargo,static,CFL,nrestart,gamma,force_kick,slope_type,viscosity,smallr,smallc,smallP

   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   CALL getarg(1,infile)
   nmlfile=TRIM(infile)
   INQUIRE(file=infile,exist=nml_ok)

   if(.not. nml_ok)then
      write(*,*)'File '//TRIM(infile)//' does not exist'
      stop
   end if
   print *, "Namelist reading : ", nmlfile
   
   !read(*,*) nmlfile
   !Reading the namelist
   open(13,file=nmlfile)
   print *, "Grid namelist reading  !"
   read(13,grid_params,IOSTAT=io)
   rewind(13)
   if (io/=0) then
      write(*,*) 'Invalid line in the grid namelist'
      stop
   end if

   if(static) then
      print*, 'Static gas'
   else
      print *, 'Dynamics activated'
   endif

#if NDUST>0
   print *, "Dust namelist reading  !"
   call read_dust_params(13,nmlfile)
#endif
#if TURB>0
   print *, "Turb namelist reading  !"
   call read_turb_params(13,nmlfile)
#endif
   print *, "Gravity namelist reading  !"
   call read_gravity_params(13,nmlfile)
   print *, "Ionisation namelist reading  !"
   call read_ionisation_params(13,nmlfile)
   print *, "Setup namelist reading  !"
   call read_setup_params(13,nmlfile)
   close(13)
   print *, "Namelist reading ok !"
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"

   
end subroutine read_params



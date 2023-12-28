subroutine read_turb_params(ilun,nmlfile)

  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  integer :: n
   integer,allocatable :: seed(:)
  namelist/turb_params/turb_turnover,v_rms,turb_correl
  print *, "########################################################################################################################################"
  print *, "########################################################################################################################################"
  print *, "turb namelist reading  !"
  read(13,turb_params,IOSTAT=io)
  rewind(13)
  if (io/=0) then
     write(*,*) 'Invalid line in the turb namelist'
     stop
  end if

   call random_seed(size=n)
   allocate(seed(n))
   seed = 123456789    ! putting arbitrary seed to all elements
   call random_seed(put=seed)
   deallocate(seed)

end subroutine read_turb_params
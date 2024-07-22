subroutine read_turb_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/turb_params/new_seed,nb_turb_modes_driven,Mach_nb_ini,Mach_yz_target,Mach_x_target,driven_turb,turb_compressive,turb_solenoidal,k_min,k_max,corrector,corrector_sol,iseed,turnover_time,phase_drift,turb_dust
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "turb_params namelist reading  !"
   read(13,turb_params,IOSTAT=io)
   rewind(13)
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   if (io/=0) then
      write(*,*) 'Invalid line in turb namelist'
      stop
   end if
   
 end subroutine read_turb_params




 ! Generates evenly spaced numbers from `from` to `to` (inclusive).
!
! Inputs:
! -------
!
! from, to : the lower and upper boundaries of the numbers to generate
!
! Outputs:
! -------
!
! array : Array of evenly spaced numbers
!
subroutine linspace(from, to, array)

  use parameters
  use commons
  use units
  use precision

  implicit none
    real(dp), intent(in) :: from, to
    real(dp), dimension(1:nb_turb_modes_driven), intent(out) :: array
    real(dp) :: range
    integer :: n, i
    n = size(array)
    range = to - from

    if (n == 0) return

    if (n == 1) then
        array(1) = from
        return
    end if


    do i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do
end subroutine linspace
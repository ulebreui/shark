subroutine read_output(iout)
  use parameters
  use commons
  use units
  implicit none
  integer  :: i,idust,jdust,iout
  real(dp) :: sum_dust,barotrop,B_field,read_var
  character(LEN = 5) :: nchar
  character(len=80)  :: path, makedirectory
 !Careful restart is not implemented.
end subroutine read_output


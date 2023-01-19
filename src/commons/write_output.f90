! Writes the simulation output. The routine is a bit ugly and will be rewritten soon
subroutine output(iout)
  use parameters
  use commons
  use units
  implicit none
  integer :: i,idust,jdust,iout
  real(dp):: sum_dust,barotrop,B_field
  character(LEN = 5) :: nchar
  character(len=80)  :: path, makedirectory

  path='output_'

  call title(iout,nchar)
  makedirectory = 'mkdir ' // trim(path) // trim(nchar)
  call system(makedirectory)
  ! Generic info file
  open(14,file=trim(path) // trim(nchar)//trim('/info.dat'))
  write(14,*) ncells-nghost*2
  call write_setup_info(14)
  close(14)
  open(14,file=trim(path) // trim(nchar)//trim('/x.dat'))
  do i = 1,ncells
   if(active_cell(i)==1) write(14,*) position(i,1)
  end do
  close(14)

  open(14,file=trim(path) // trim(nchar)//trim('/rho.dat'))
  do i = 1,ncells
   if(active_cell(i)==1) write(14,*) q(i,irho)*unit_d
  end do
  close(14)
  open(14,file=trim(path) // trim(nchar)//trim('/v.dat'))
  do i = 1,ncells
   if(active_cell(i)==1) write(14,*) q(i,iv)*unit_v
  end do
  close(14)
   if(ndim==2) then
   open(14,file=trim(path) // trim(nchar)//trim('/y.dat'))
   do i = 1,ncells
      if(active_cell(i)==1) write(14,*) position(i,2)
   end do
   close(14)
   open(14,file=trim(path) // trim(nchar)//trim('/vy.dat'))
   do i = 1,ncells
      if(active_cell(i)==1) write(14,*) q(i,ivy)*unit_v
   end do
   close(14)
   endif
  open(14,file=trim(path) // trim(nchar)//trim('/P.dat'))
  do i = 1,ncells
   if(active_cell(i)==1) write(14,*) q(i,iP)*unit_P
  end do
  close(14)

 end subroutine output

!===========================================================================
subroutine title(n,nchar) ! Shamelessely stolen from RAMSES (Teyssier, 2002)
!===========================================================================
  implicit none
  integer::n
  character(LEN=5)::nchar

  character(LEN=1)::nchar1
  character(LEN=2)::nchar2
  character(LEN=3)::nchar3
  character(LEN=4)::nchar4
  character(LEN=5)::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title


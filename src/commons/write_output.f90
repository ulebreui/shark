! Writes the simulation output. The routine is a bit ugly and will be rewritten soon
subroutine output(iout)
  use parameters
  use commons
  use units
  implicit none
  integer :: i,idust,jdust,iout,ilun=50
  real(dp):: sum_dust,barotrop,B_field
  real(dp), dimension(1: nx*ny) :: xdp
  character(LEN = 5) :: nchar
  character(len=80)  :: path, makedirectory,format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)
  makedirectory = 'mkdir ' // trim(path) // trim(nchar)
  call system(makedirectory)
  ! Generic info file
  open(ilun,file=trim(path) // trim(nchar)//trim('/info.dat'))
  write(ilun,*) nvar
  write(ilun,*) ncells-nghost*2
  call write_setup_info(ilun)
  close(ilun)
  open(ilun,file=trim(path) // trim(nchar)//trim('/x'), form=format_out,access='stream')
  do i = 1,ncells
   if(active_cell(i)==1) write(ilun) position(i,1)
  end do
  !xdp
  close(ilun)

  open(ilun,file=trim(path) // trim(nchar)//trim('/rho'), form=format_out,access='stream')
  do i = 1,ncells
   if(active_cell(i)==1) write(ilun) q(i,irho)*unit_d
  end do
  close(ilun)  
  open(ilun,file=trim(path) // trim(nchar)//trim('/v'), form=format_out,access='stream')
  do i = 1,ncells
   if(active_cell(i)==1) write(ilun) q(i,iv)*unit_v
  end do
  close(ilun)
   if(ndim==2) then
   open(ilun,file=trim(path) // trim(nchar)//trim('/y'), form=format_out,access='stream')
   do i = 1,ncells
      if(active_cell(i)==1) write(ilun) position(i,2)
   end do
   close(ilun)
   open(ilun,file=trim(path) // trim(nchar)//trim('/vy'), form=format_out,access='stream')
   do i = 1,ncells
      if(active_cell(i)==1) write(ilun) q(i,ivy)*unit_v
   end do
   close(ilun)
   endif
  open(ilun,file=trim(path) // trim(nchar)//trim('/P'), form=format_out,access='stream')
  do i = 1,ncells
   if(active_cell(i)==1) write(ilun) q(i,iP)*unit_P
  end do
  close(ilun)

#if NDUST>0
  open(ilun,file=trim(path) // trim(nchar)//trim('/rhod'), form=format_out,access='stream')
  do idust=1,ndust
   do i = 1,ncells
      if(active_cell(i)==1) write(ilun) q(i,irhod(idust))*unit_d
   end do
  end do
  close(ilun)
  open(ilun,file=trim(path) // trim(nchar)//trim('/vd'), form=format_out,access='stream')
  do idust=1,ndust
   do i = 1,ncells
      if(active_cell(i)==1) write(ilun) q(i,ivd(idust))*unit_v
   end do
  end do
  close(ilun)
  open(ilun,file=trim(path) // trim(nchar)//trim('/vdy'), form=format_out,access='stream')
  do idust=1,ndust
   do i = 1,ncells
      if(active_cell(i)==1) write(ilun) q(i,ivdy(idust))*unit_v
   end do
  end do
  close(ilun)
#endif
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


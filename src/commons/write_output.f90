! Writes the simulation output. The routine is a bit ugly and will be rewritten soon
subroutine output(iout)
  use parameters
  use commons
  use units
  implicit none
  integer  :: i,idust,jdust,iout,ilun=50,ipscal,ixx,iyy,ix,iy
  real(dp) :: sum_dust,barotrop,B_field,rhod_tot
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
  write(ilun,'("nvar        =",I11)')nvar
  write(ilun,'("ndust        =",I11)')NDUST
  write(ilun,'("ndustpscal   =",I11)')NDUSTPSCAL
  write(ilun,'("NX        =",I11)')NX
  write(ilun,'("NY        =",I11)')NY
  write(ilun,'("GEOM        =",I11)')GEOM
  write(ilun,'("GRIDSPACE       =",I11)')GRIDSPACE
  write(ilun,'("time        =",E23.15)')time
  write(ilun,'("unit_t        =",E23.15)')unit_t
  write(ilun,'("unit_d        =",E23.15)')unit_d
  write(ilun,'("unit_l        =",E23.15)')unit_l
  write(ilun,'("unit_v        =",E23.15)')unit_v
  write(ilun,'("unit_p        =",E23.15)')unit_p
  !call write_setup_info(ilun)
  close(ilun)
  
  open(ilun,file=trim(path) // trim(nchar)//trim('/x'), form=format_out,access='stream')
  do i = 1,ncells
#if GEOM<2
   if(active_cell(i)==1) then
      ix=ixx(i)
      iy=iyy(i)
      write(ilun) position(ix,iy,1) !
   endif 
#endif
#if GEOM==2
   if(active_cell(i)==1) write(ilun) radii_c(i) !
#endif
#if GEOM==4
   if(active_cell(i)==1) write(ilun) radii_c(i) !
#endif
  end do
  !xdp
  close(ilun)

  open(ilun,file=trim(path) // trim(nchar)//trim('/rho'), form=format_out,access='stream')

  do iy = first_active_y,last_active_y
      do ix = first_active,last_active
#if GEOM<2
         write(ilun) q(ix,iy,irho)*unit_d
#endif
#if GEOM==2
         write(ilun) q(ix,iy,irho)*unit_dcol
#endif
#if GEOM>2
         write(ilun) q(ix,iy,irho)*unit_d
#endif
      end do 
  end do
  close(ilun)  

  open(ilun,file=trim(path) // trim(nchar)//trim('/v'), form=format_out,access='stream')

  do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         write(ilun) q(ix,iy,ivx)*unit_v !
       end do
  end do
  close(ilun)

   if(ndim==2) then

   open(ilun,file=trim(path) // trim(nchar)//trim('/y'), form=format_out,access='stream')
   do i = 1,ncells
      if(active_cell(i)==1) then
         ix=ixx(i)
         iy=iyy(i)
         write(ilun) position(ix,iy,2)
      endif
   end do
   close(ilun)
   endif

   open(ilun,file=trim(path) // trim(nchar)//trim('/vy'), form=format_out,access='stream')
   do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         write(ilun) q(ix,iy,ivy)*unit_v
      end do
   end do
   close(ilun)
   open(ilun,file=trim(path) // trim(nchar)//trim('/vz'), form=format_out,access='stream')

   do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         write(ilun) q(ix,iy,ivz)*unit_v
      end do
   end do

   close(ilun)
  open(ilun,file=trim(path) // trim(nchar)//trim('/P'), form=format_out,access='stream')
   do iy = first_active_y,last_active_y
      do ix = first_active,last_active

#if GEOM<2
         write(ilun) q(ix,iy,iP)*unit_P
#endif
#if GEOM==2
         write(ilun) q(ix,iy,iP)*unit_P*unit_l
#endif
#if GEOM>2
         write(ilun) q(ix,iy,iP)*unit_P
#endif
      end do
  end do
  close(ilun)

#if NDUST>0
  open(ilun,file=trim(path) // trim(nchar)//trim('/rhod_tot'), form=format_out,access='stream')

   do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         rhod_tot=0.d0
         do idust=1,ndust
#if GEOM<2
            rhod_tot = rhod_tot + q(ix,iy,irhod(idust))*unit_d
#endif
#if GEOM==2
            rhod_tot  = rhod_tot +  q(ix,iy,irhod(idust))*unit_dcol
#endif
#if GEOM>2
            rhod_tot  = rhod_tot + q(ix,iy,irhod(idust))*unit_d
#endif
         end do
         write(ilun) rhod_tot
      end do 
   end do

  close(ilun)
  open(ilun,file=trim(path) // trim(nchar)//trim('/rhod'), form=format_out,access='stream')
  do idust=1,ndust
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
#if GEOM<2
            write(ilun) q(ix,iy,irhod(idust))*unit_d
#endif
#if GEOM==2
            write(ilun) q(ix,iy,irhod(idust))*unit_dcol
#endif
#if GEOM>2
            write(ilun) q(ix,iy,irhod(idust))*unit_d
#endif
         end do
      end do
  end do
  close(ilun)

#if NDUSTPSCAL>0
 open(ilun,file=trim(path) // trim(nchar)//trim('/dustpscal'), form=format_out,access='stream')

 do idust=1,ndust
   do ipscal=1,ndustpscal
     do idust=1,ndust
         do iy = first_active_y,last_active_y
            write(ilun) q(ix,iy,idust_pscal(idust,ipscal))
      end do
   end do
 end do
 close(ilun)

#endif

  open(ilun,file=trim(path) // trim(nchar)//trim('/sd'), form=format_out,access='stream')
  do idust=1,ndust
   do i = 1,ncells
      if(active_cell(i)==1) write(ilun) sdust(i,idust)*unit_l 
   end do
  end do
  close(ilun)

  open(ilun,file=trim(path) // trim(nchar)//trim('/vd'), form=format_out,access='stream')
  do idust=1,ndust
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
            write(ilun) q(ix,iy,ivdx(idust))*unit_v
         end do
      end do
  end do
  close(ilun)

  open(ilun,file=trim(path) // trim(nchar)//trim('/vdy'), form=format_out,access='stream')

  do idust=1,ndust
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
            write(ilun) q(ix,iy,ivdy(idust))*unit_v
         end do
      end do
  end do

  close(ilun)
  open(ilun,file=trim(path) // trim(nchar)//trim('/vdz'), form=format_out,access='stream')
  do idust=1,ndust
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
            write(ilun) q(ix,iy,ivdy(idust))*unit_v
         end do
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


! subroutine write_variable(name_var,xdp)

!   use parameters
!   use commons
!   use units
!   implicit none
!   integer  :: i,ilun=50
!   real(dp), dimension(1: ncells) :: xdp
!   open(ilun,file=trim(path) // trim(nchar)//trim(name_var), form=format_out,access='stream')
!   do i = 1,ncells
!    if(active_cell(i)==1) write(ilun) xdp(i)
!   end do
!   !xdp
!   close(ilun)
! end subroutine write_variable

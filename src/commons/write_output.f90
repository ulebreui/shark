! Writes the simulation output. The routine is a bit ugly and will be rewritten soon
subroutine output(iout)
  use parameters
  use commons
  use units
  implicit none
  integer  :: i,idust,jdust,iout,ilun=50,ipscal,ixx,iyy,ix,iy
  real(dp) :: sum_dust,barotrop,B_field,rhod_tot,unit_dloc,unit_Ploc
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

  write(ilun,'("nvar        =",I11)')      nvar
  write(ilun,'("ndust        =",I11)')     NDUST
  write(ilun,'("ndustpscal   =",I11)')     NDUSTPSCAL
  write(ilun,'("NX        =",I11)')        NX
  write(ilun,'("NY        =",I11)')        NY
  write(ilun,'("GEOM        =",I11)')      GEOM
  write(ilun,'("GRIDSPACE       =",I11)')  GRIDSPACE
  write(ilun,'("time        =",E23.15)')   time
  write(ilun,'("unit_t        =",E23.15)') unit_t
  write(ilun,'("unit_d        =",E23.15)') unit_d
  write(ilun,'("unit_l        =",E23.15)') unit_l
  write(ilun,'("unit_v        =",E23.15)') unit_v
  write(ilun,'("unit_p        =",E23.15)') unit_p
  !call write_setup_info(ilun)
  close(ilun)
  
  unit_dloc = unit_d
  unit_Ploc = unit_P

#if GEOM==2
  unit_dloc = unit_dcol
  unit_Ploc =unit_P*unit_l
#endif
 
  call write_backup(iout)
  call write_x(iout)
  call write_y(iout)

  call write_rho(iout,unit_dloc)
  call write_v(iout,unit_v)
  call write_vy(iout,unit_v)
  call write_vz(iout,unit_v)
  call write_P(iout,unit_Ploc)


#if NDUST>0
  call write_rhod(iout,unit_dloc)
  call write_vdx(iout,unit_v)
  call write_vdy(iout,unit_v)
  call write_vdz(iout,unit_v)
  call write_sdust(iout,unit_l)
  call write_St(iout)


#if NDUSTPSCAL>0
  call write_rho_pscal(iout)
#endif
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

subroutine write_x(iout)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout
  character(LEN = 5) :: nchar
  character(len=80)  :: path,format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)
  open(ilun,file=trim(path) // trim(nchar)//trim('/x'), form=format_out,access='stream')

   do iy = first_active_y,last_active_y
      do ix = first_active,last_active
#if GEOM<2
         write(ilun) position(ix,iy,1) !
#endif
#if GEOM==2
         write(ilun) radii(ix,iy) !
#endif
#if GEOM==4
      write(ilun) radii(ix,iy) !
#endif
  end do
  end do
  close(ilun)

end subroutine write_x


subroutine write_y(iout)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout
  character(LEN = 5) :: nchar
  character(len=80)  :: path,format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)
  open(ilun,file=trim(path) // trim(nchar)//trim('/y'), form=format_out,access='stream')

    do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         write(ilun) position(ix,iy,2)
      end do
  end do
  close(ilun)

end subroutine write_y


subroutine write_rho(iout,unit_loc)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out
  real(dp)           :: unit_loc

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/rho'), form=format_out,access='stream')

  do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         write(ilun) q(irho,ix,iy)*unit_loc
      end do 
  end do

  close(ilun)
end subroutine write_rho


subroutine write_v(iout,unit_loc)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out
  real(dp) :: unit_loc

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/v'), form=format_out,access='stream')

  do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         write(ilun) q(ivx,ix,iy)*unit_loc
      end do 
  end do

  close(ilun)
end subroutine write_v

subroutine write_vy(iout,unit_loc)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out
  real(dp) :: unit_loc

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/vy'), form=format_out,access='stream')

  do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         write(ilun) q(ivx,ix,iy)*unit_loc
      end do 
  end do

  close(ilun)
end subroutine write_vy

subroutine write_vz(iout,unit_loc)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out
  real(dp) :: unit_loc

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/vz'), form=format_out,access='stream')

  do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         write(ilun) q(ivx,ix,iy)*unit_loc
      end do 
  end do

  close(ilun)
end subroutine write_vz


subroutine write_P(iout,unit_loc)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out
  real(dp) :: unit_loc
  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)
  open(ilun,file=trim(path) // trim(nchar)//trim('/P'), form=format_out,access='stream')
   do iy = first_active_y,last_active_y
      do ix = first_active,last_active

         write(ilun) q(iP,ix,iy)*unit_P*unit_loc

      end do
  end do
  close(ilun)

  end subroutine write_P

#if NDUST>0
subroutine write_rhod(iout,unit_loc)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout,idust
  real(dp) :: rhod_tot
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out
  real(dp) :: unit_loc

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/rhod'), form=format_out,access='stream')
  do idust=1,ndust
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
            write(ilun) q(irhod(idust),ix,iy)*unit_loc
         end do
      end do
  end do
  close(ilun)

  open(ilun,file=trim(path) // trim(nchar)//trim('/rhod_tot'), form=format_out,access='stream')

  do iy = first_active_y,last_active_y
   do ix = first_active,last_active
      rhod_tot=0.d0
         do idust=1,ndust
            rhod_tot  = rhod_tot +  q(irhod(idust),ix,iy)*unit_loc
         end do
         write(ilun) rhod_tot
      end do 
   end do

  close(ilun)
  end subroutine write_rhod

#if NDUSTPSCAL>0
subroutine write_rho_pscal(iout)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout,idust
  real(dp) :: rhod_tot
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

 open(ilun,file=trim(path) // trim(nchar)//trim('/dustpscal'), form=format_out,access='stream')

 do idust=1,ndust
   do ipscal=1,ndustpscal
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
            write(ilun) q(idust_pscal(idust,ipscal),ix,iy)
      end do
   end do
 end do
 close(ilun)
 end subroutine write_rho_pscal

#endif


subroutine write_sdust(iout,unit_loc)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout,idust
  real(dp) :: rhod_tot,unit_loc
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/sd'), form=format_out,access='stream')
   do idust=1,ndust
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
            write(ilun) sdust(idust)*unit_loc
      end do
   end do
 end do
  close(ilun)


end subroutine write_sdust

subroutine write_vdx(iout,unit_loc)


  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout,idust
  real(dp) :: unit_loc
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/vd'), form=format_out,access='stream')
  do idust=1,ndust
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
            write(ilun) q(ivdx(idust),ix,iy)*unit_loc
         end do
      end do
  end do
  close(ilun)
end subroutine write_vdx


subroutine write_vdy(iout,unit_loc)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout,idust
  real(dp) :: unit_loc

  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/vdy'), form=format_out,access='stream')
  do idust=1,ndust
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
            write(ilun) q(ivdy(idust),ix,iy)*unit_loc
         end do
      end do
  end do
  close(ilun)
end subroutine write_vdy

subroutine write_vdz(iout,unit_loc)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout,idust
  real(dp) :: unit_loc

  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)
  open(ilun,file=trim(path) // trim(nchar)//trim('/vdz'), form=format_out,access='stream')
  do idust=1,ndust
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
            write(ilun) q(ivdz(idust),ix,iy)*unit_loc
         end do
      end do
  end do
  close(ilun)
end subroutine write_vdz


subroutine write_St(iout)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout,idust
  real(dp) :: rhod_tot
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/St'), form=format_out,access='stream')
  do idust=1,ndust
      do iy = first_active_y,last_active_y
         do ix = first_active,last_active
            write(ilun) St(idust,ix,iy)
         end do
      end do
  end do
  close(ilun)
end subroutine write_St




#endif


subroutine write_backup(iout)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout,ivar
  real(dp) :: rhod_tot
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/uprim'), form=format_out,access='stream')
  do ivar=1,nvar
      do iy = 1,ny_max
         do ix = 1,nx_max
            write(ilun) u_prim(ivar,ix,iy)
         end do
      end do
  end do
  close(ilun)

  open(ilun,file=trim(path) // trim(nchar)//trim('/time'), form=format_out,access='stream')
  write(ilun) time
  close(ilun)


  end subroutine write_backup


 subroutine restart(iout)

  use parameters
  use commons
  use units
  implicit none
  integer  :: i,ilun,ix,iy,iout,ivar
  real(dp) :: xdp
  character(LEN = 5) :: nchar
  character(len=80)  :: path, format_out

  path='output_'
  format_out=trim("unformatted")
  ilun=20
  call title(iout,nchar)

  open(ilun,file=trim(path) // trim(nchar)//trim('/uprim'), form=format_out,access='stream')
  do ivar=1,nvar
      do iy = 1,ny_max
         do ix = 1,nx_max
            read(ilun) xdp
            u_prim(ivar,ix,iy) = xdp
         end do
      end do
  end do
  close(ilun)

  open(ilun,file=trim(path) // trim(nchar)//trim('/time'), form=format_out,access='stream')
  read(ilun) xdp
  close(ilun)
  time = xdp

  call ctoprim
  call apply_boundaries
  end subroutine restart

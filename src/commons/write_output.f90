! Writes the simulation output. The routine is a bit ugly and will be rewritten soon
subroutine output(iout)
  use parameters
  use commons
  use units
  implicit none
  integer :: i,idust,jdust,iout
  real(dp):: sum_dust,barotrop,B_field
  character(LEN = 5) :: nchar
  character(len=80) :: path, makedirectory

  path='output_'

  call title(iout,nchar)
  makedirectory = 'mkdir ' // trim(path) // trim(nchar)
  call system(makedirectory)
  ! Generic info file
  open(14,file=trim(path) // trim(nchar)//trim('/info.dat'))
  write(14,*) ncells-nghost*2
  call write_setup_info(14)
  ! Gas dump
  open(15,file=trim(path) // trim(nchar)//trim('/grid.dat'))
  open(16,file=trim(path) // trim(nchar)//trim('/rho.dat'))
  open(17,file=trim(path) // trim(nchar)//trim('/v.dat'))
  open(33,file=trim(path) // trim(nchar)//trim('/nh.dat'))
  open(34,file=trim(path) // trim(nchar)//trim('/temperature.dat'))
  open(35,file=trim(path) // trim(nchar)//trim('/B.dat'))

  do i =nghost+1,ncells-nghost
     write(15,*) r_c(i)*unit_l/au
     write(16,*) q(i,irho)*unit_d
     write(17,*) q(i,iv)*unit_v
     write(33,*) q(i,irho)*unit_d/(mu_gas*mH)
     write(34,*) barotrop(q(i,irho))
#if MHD==0     
!     write(35,*) B_field(q(i,irho))
     write(35,*) B_cell(i)
     
#endif     
  end do
  close(14)
  close(15)
  close(16)
  close(17)
  
#if NDUST>0
  !Dust info file
  open (30,file=trim(path) // trim(nchar)//trim('/info_dust.dat'))
  write(30,*) ndust
  write(30,*) smin
  write(30,*) smax
  write(30,*) scutmin
  write(30,*) scut  
  write(30,*) mrn

  !Dust dump
  open(18,file=trim(path) // trim(nchar)//trim('/epsilon.dat'))
  open(19,file=trim(path) // trim(nchar)//trim('/rhodust.dat'))
  open(20,file=trim(path) // trim(nchar)//trim('/sdust.dat'))
  open(21,file=trim(path) // trim(nchar)//trim('/mdust.dat'))
  open(22,file=trim(path) // trim(nchar)//trim('/epsdust.dat'))
  open(23,file=trim(path) // trim(nchar)//trim('/vdust.dat'))
  open(24,file=trim(path) // trim(nchar)//trim('/eta_a.dat'))
  open(25,file=trim(path) // trim(nchar)//trim('/eta_h.dat'))
  open(26,file=trim(path) // trim(nchar)//trim('/eta_o.dat'))
  open(27,file=trim(path) // trim(nchar)//trim('/ni.dat'))
  open(28,file=trim(path) // trim(nchar)//trim('/ne.dat'))
  open(29,file=trim(path) // trim(nchar)//trim('/zd.dat'))
  open(31,file=trim(path) // trim(nchar)//trim('/mminus.dat'))
  open(32,file=trim(path) // trim(nchar)//trim('/mplus.dat'))
  open(36,file=trim(path) // trim(nchar)//trim('/sigma_o.dat'))
  open(37,file=trim(path) // trim(nchar)//trim('/sigma_p.dat'))
  open(38,file=trim(path) // trim(nchar)//trim('/sigma_h.dat'))
  open(39,file=trim(path) // trim(nchar)//trim('/aminus.dat'))
  open(40,file=trim(path) // trim(nchar)//trim('/aplus.dat'))
  open(41,file=trim(path) // trim(nchar)//trim('/dust_drifts.dat'))
  open(42,file=trim(path) // trim(nchar)//trim('/gamma_d.dat'))

  do i =nghost+1,ncells-nghost
     write(24,*) eta_a(i)!/unit_t
     write(25,*) eta_h(i)!/unit_t
     write(26,*) eta_o(i)!/unit_t
     write(27,*) ni(i)
     write(28,*) ne(i)
     write(36,*) sigma_o(i)
     write(37,*) sigma_p(i)
     write(38,*) sigma_h(i)

     sum_dust=0.0
     do idust=1,ndust
        sum_dust=sum_dust + epsilondust(i,idust)
     end do
     write(18,*) sum_dust
  end do
  do idust=1,ndust
     do i =nghost+1,ncells-nghost
        write(19,*) q(i,irhod(idust))*unit_d
        write(20,*) sdust(i,idust)*unit_l
        write(21,*) mdust(i,idust)*unit_m
        write(22,*) epsilondust(i,idust)
        write(23,*) q(i,ivd(idust))*unit_v
        write(29,*) zd(i,idust)
        write(42,*) gamma_d(i,idust)
     end do
  end do
  do idust=1,ndust
     write(31,*) mminus(idust)*unit_m
     write(32,*) mplus(idust)*unit_m
     write(39,*) aminus(idust)*unit_l
     write(40,*) aplus(idust)*unit_l
     
  end do
  !Write the drift velocities
  do i =nghost+1,ncells-nghost
     do idust=1,ndust
        do jdust=1,ndust
           write(41,*) vdrift_turb(i,idust,jdust)*unit_v, vdrift_hydro(i,idust,jdust)*unit_v,vdrift_ad(i,idust,jdust)*unit_v, vdrift_brow(i,idust,jdust)*unit_v
        end do
     end do
  end do
  
  close(18)
  close(19)
  close(20)
  close(21)
  close(22)
  close(23)
  close(24)
  close(25)
  close(26)
  close(27)
  close(28)
  close(29)
  close(30)
  close(31)
  close(32)
  close(33)
  close(34)
  close(35)
  close(36)
  close(37)
  close(38)
  close(39)
  close(40)
  close(41)
  close(42)
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


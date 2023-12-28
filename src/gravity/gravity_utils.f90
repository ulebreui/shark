!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This routine computes the cumulative mass in the cells Mc which is needed for the
!gravity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mtot
  use parameters
  use commons
  use units
  implicit none
  integer :: i,idust
  real(dp), dimension(1:ncells) :: Mint,dmint
  real(dp):: fourpi

  fourpi = 4.0d0*pi
  Mint   = 0.0d0
  dmint  = 0.0d0
  M_tot  = 0.0d0

#if NDUST>0  
  M_tot_dust =0.0d0
#endif  

#if NY==1
  do i=first_active,last_active
     Mint(i)  = q(i,irho)*vol(i)*fourpi
     M_tot    = M_tot+Mint(i)
     dmint(i) = q(i,irho)*dvol(i)*fourpi
#if NDUST>0     
     do idust=1,ndust
        Mint(i)    = Mint(i)+q(i,irhod(idust))*vol(i)*fourpi
        M_tot_dust = M_tot_dust+q(i,irhod(idust))*vol(i)*fourpi
        dmint(i)   = dmint(i)+q(i,irhod(idust))*dvol(i)*fourpi
     end do
#endif     
  end do
  !Compute mass within cells
  Mc(first_active-1) = M_cent

  do i=first_active,last_active
     Mc(i) = Mc(i-1)+Mint(i)
  end do
  Mc = Mc-dmint
#endif

end subroutine mtot


subroutine read_gravity_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/gravity_params/l_soft,self_gravity
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "Gravity namelist reading  !"
   read(13,gravity_params,IOSTAT=io)
   rewind(13)
   print *, 'Gravity is used'
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   if (io/=0) then
      write(*,*) 'Invalid line in gravity namelist'
      stop
   end if
   
 end subroutine read_gravity_params

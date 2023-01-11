subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: rho_cloud,r_cloud,mcloud
  real(dp) :: rmax,vol_tot,zeta_nh

  integer :: i,idust,imax
  zero_flux_in=.true. ! Mandatory for collapse calculations

  tend=0.0d0
  call allocate_init
  zeta_nh=(nh_max/(nh_min))**(1.0d0/(ncells-2*nghost-1))

  !initialisation of the hydro variables
  do i=nghost+1,ncells-nghost
     q(i,irho)=(nh_min)*zeta_nh**(i-1)
  end do

#if NDUST>0
  call distribution_dust(.true.)
  do i=1,ncells
     do idust=1,ndust
        q(i,irhod(idust))= epsilondust(i,idust)*q(i,irho)
     end do
  end do
#endif
  call primtoc


  call Mtot
  print *, 'Total mass is',M_tot*unit_m/2d33, 'Solar masses'
#if NDUST>0
  call charge 
#endif
  
end subroutine setup

subroutine gridinit(rmax)
  use parameters
  use commons
  use units
  implicit none

  real(dp):: rmax,zeta_r,ms
  integer :: i

  !Grid initialisation
  print *,' Logarithmic grid'
  zeta_r=(rmax/(rin*au))**(1.0d0/(ncells-2*nghost-1))
  r=0.0d0
  do i=1,ncells-nghost
     r(nghost+i)= (rin*au/unit_l)*zeta_r**(i-1)
  enddo
  do i=nghost+1,ncells-nghost
     dx(i)=r(i)-r(i-1)
  end do
  do i=1,nghost+1
     dx(i)=dx(nghost+1)
  end do

  !Cell center
  do i=nghost+1,ncells!-nghost
     r_c(i)     = ( (r(i)**3 + r(i-1)**3) / 2.)**(1./3.)
  end do

 !Cell volume
  do i=nghost+1,ncells-nghost
     vol_cell(i) = (r(i)**3.-r(i-1)**3.)/3.0d0
  end do
  !vol_cell(nghost+1)=r(nghost+1)**3/3.0d0

  !Cell volume from center to the right 
  do i=nghost+1,ncells-nghost
     dvol(i)  =(r(i)**3-(r_c(i)**3) )/3.0d0
   end do

  do i=nghost+1,ncells-nghost
     Surf_p(i) = r(i)**2.
  end do
  Surf_p(nghost)=0.0d0
  !if(sink)Surf_p(nghost)=((rin*au)/unit_l)**2.
  !Surf_p(ncells-nghost)= r(ncells-nghost)**2.

  !Needed for 2nd order
  
  do i=nghost+1,ncells-1!-nghost
     dx_c(i)      = r_c(i+1)- r_c(i-1)
     dx_r(i)      = r_c(i+1)- r_c(i)
     dx_l(i)      = r_c(i)- r_c(i-1)
     dx_r_cell(i) = r_c(i+1)- r(i)
     dx_l_cell(i) = r(i)- r_c(i)
  end do
  do i=1,nghost
     dx_l(i)      = r_c(nghost+1)
     dx_l_cell(i) = r_c(nghost+1)
     dx_r(i)      = dx_r(nghost+1)
     dx_r_cell(i) = dx_r_cell(nghost+1)
     !dx_r(ncells+1-i)=dx_r(ncells-nghost)
     !dx_r_cell(ncells+1-i)=dx_r(ncells-nghost)*0.5d0
     !dx_l_cell(ncells+1-i)=dx_l(ncells-nghost)*0.5d0
  end do
  
  
  !Not dx_l etc at the boundary dont really matter, they are not really used
  print *, 'r'
  print *, r
  print *, 'rc'
  print *, r_c
  print *, 'dx'
  print *, dx
  print *, 'vol_cell'
  print *, vol_cell
  print *, 'dvol'
  print *, dvol
  print *, 'dxc'
  print *, dx_c 
  print *, 'dxr'
  print *, dx_r
  print *, 'dxl'
  print *, dx_l
  print *, 'dxrcell'
  print *, dx_r_cell
  print *, 'dxlcell'
  print *, dx_l_cell
  !stop
end subroutine gridinit


subroutine write_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  write(ilun,*) time*unit_t/(365.*24.*3600.*1000.)
  write(ilun,*) alpha_cloud
end subroutine write_setup_info

subroutine read_setup_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/setup_params/M_cloud,alpha_cloud,T0_cloud,gamma,grid_type,ramses_baro,rho_max_sim,more_outputs,nh_min,nh_max
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "Setup namelist reading  !"
   read(13,setup_params,IOSTAT=io)
   rewind(13)
   if (io/=0) then
      write(*,*) 'Invalid line in the setup namelist'
      stop
   end if
   
   print *,"Mass of the cloud: "     ,  M_cloud
   print *,"Virial parameter: "      ,  alpha_cloud
   print *,"Isothermal temperature: ",  T0_cloud
   print *,"Adiabatic index: "       ,  gamma
   print*, 'Simulation will end at density', rho_max_sim
   if(ramses_baro)print *,"Barotrop as in RAMSES"
   if(.not.ramses_baro)print *,"Full barotrop"
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"

   
 end subroutine read_setup_params

 subroutine flag_continue(continue_sim)
   use parameters
   use commons
   use units
   implicit none
   logical :: continue_sim
   !Here you can add flags to kill the simulation
    continue_sim=.false.

 end subroutine flag_continue

 subroutine check_output(icount,iout,outputing,verbose)
   use parameters
   use commons
   use units
   implicit none
   integer :: icount,iout
   logical :: outputing,verbose

 end subroutine check_output

 subroutine setup_preloop
   use parameters
   use commons
   use units
   implicit none

end subroutine setup_preloop

subroutine setup_inloop
   use parameters
   use commons
   use units
   implicit none
   !call output(1)
   return
end subroutine setup_inloop

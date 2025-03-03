subroutine setup
  use parameters
  use commons
  use units
  implicit none

  real(dp) :: vol_tot,xx,yy,vkep,vx_nak,vy_nak


  real(dp) :: rho_cloud,H_cloud,mcloud,Rstar_cu,omega_cloud,cs_cloud
  real(dp) :: rmax,rr,y
  
  integer :: i, idust,jdust,imax,ix,iy,icell,ixx,iyy
  call allocate_init

  box_l=box_l*1.5d13/unit_l
  box_l_y=box_l_y*1.5d13/unit_l
  call gridinit(box_l,box_l_y)

  Rstar_cu    = Rstar*1.5d13/unit_l
  H_cloud     = Rstar_cu * HoverR
  mcloud      = M_cloud*Msun/unit_m
  omega_cloud = sqrt(mcloud/Rstar_cu**3.)
  cs_cloud    = H_cloud*omega_cloud*unit_l/unit_t
  cs          = cs_cloud/unit_v
  rho_cloud   = 6d-13/unit_d
  q = 0.0d0
  iso_cs = 1
  do i =1,ncells
        ix=ixx(i)
        iy=iyy(i)
        xx=position(ix,iy,1)-half*box_l  ! Boxlen already in pc
        yy=position(ix,iy,2)-half*box_l_y
        rr=sqrt(yy**2.+Rstar_cu**2)
        q(irho,ix,iy) = max(rho_cloud*exp(-yy**2./(2.0d0*H_cloud**2)),1d-25/unit_d)
        q(ivx,ix,iy)  = 0.0d0
        q(ivy,ix,iy)  = 0.0d0
        q(iP,ix,iy)   = q(irho,ix,iy)*cs(ix,iy)**2
  end do

  call distribution_dust
  sdust(1)=1d-5/unit_l
  sdust(2)=2.78d-5/unit_l
  sdust(3)=7.74d-5/unit_l
  sdust(4)=2.15d-4/unit_l
  sdust(5)=5.99d-4/unit_l
  sdust(6)=1.67d-3/unit_l
  sdust(7)=4.64d-3/unit_l
  sdust(8)=1.29d-2/unit_l
  sdust(9)=3.59d-2/unit_l
  sdust(10)=0.1/unit_l
  
  epsilondust(1,1)=3.99d-5
  epsilondust(1,2)=6.65d-5
  epsilondust(1,3)=1.11d-4
  epsilondust(1,4)=1.85d-4
  epsilondust(1,5)=3.09d-4
  epsilondust(1,6)=5.15d-4
  epsilondust(1,7)=8.59d-4
  epsilondust(1,8)=1.43d-3
  epsilondust(1,9)=2.39d-3
  epsilondust(1,10)=3.99d-3

  do i =1,ncells
        ix=ixx(i)
        iy=iyy(i)
        do idust=1,ndust
          q(irhod(idust),ix,iy)= epsilondust(1,idust)*q(irho,ix,iy)
        end do
  end do
call primtoc
call apply_boundaries
call compute_tstop
call update_force_setup
  tend =tend*2.0d0*3.14159265d0/omega_cloud

end subroutine setup



subroutine get_rhoturb(pert,del)

  use random
  use precision
  implicit none

  integer :: i
  integer::iseed=0         
  integer ,dimension(1,1:IRandNumSize)    :: allseed
  integer,dimension(IRandNumSize) :: localseed=-1
  real(dp) ::  pert,del,randno

  if (localseed(1)==-1) then
     call rans(1,iseed,allseed)
     localseed = allseed(1,1:IRandNumSize)
  end if

  call ranf(localseed,randno)

  del = randno*pert - 0.5d0*pert


end subroutine get_rhoturb 

subroutine write_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  write(ilun,*) time
end subroutine write_setup_info


subroutine read_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  real(dp):: info1
  
  read(ilun,*) info1
  time=info1
  close(ilun)
end subroutine read_setup_info


subroutine read_setup_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/setup_params/box_l,box_l_y
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "Setup namelist reading  !"
   read(13,setup_params,IOSTAT=io)
   rewind(13)
   if (io/=0) then
      write(*,*) 'Invalid line in the setup namelist'
      stop
   end if
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
   if(time>=tend) continue_sim=.false.

 end subroutine flag_continue

 subroutine check_output(icount,iout,outputing,verbose)
   use parameters
   use commons
   use units
   implicit none
   integer :: icount,iout,i,idust,iy,ix
   real(dp):: mtot,mdtot
   logical :: outputing,verbose

        !We make an output at a certain frequency rate of for specific values of the density. This can be tuned at will
     if(icount.eq.freq_out) then
        print *, "time =",time,' tend = ', tend
        verbose=.true.
        outputing=.true.
        icount=0
     endif

     if(time.eq.0) outputing = .true.   
     if(outputing) call output(iout)
     if(outputing) iout=iout+1
     if(outputing) print *, "Outputing data "

     if(outputing) then
        mtot=0.0d0
        mdtot=0.0d0
        do iy = first_active_y,last_active_y
          do ix = first_active,last_active
          mtot=mtot+u_prim(irho,ix,iy)
          do idust=1,ndust
            mdtot=mdtot+u_prim(irhod(idust),ix,iy)
          end do
          end do
        end do
     end if

     if(outputing) print *, " Total mass is", mtot
     if(outputing) print *, " Total dust mass is", mdtot

     outputing=.false.

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
   use OMP_LIB
   implicit none
   integer :: idust,ix,iy
   real(dp) :: Ohmdt
   real(dp) ::  u, v, AA, BB,d
   

end subroutine setup_inloop


 subroutine update_force_setup
   use parameters
   use commons
   use units
   use OMP_LIB
   implicit none

   integer ::idust,ix,iy
   real(dp):: xx,yy,rr,Rstar_cu
    Rstar_cu    = Rstar*1.5d13/unit_l

   if(.not.force_kick) return
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy)
    do iy = first_active_y,last_active_y
       do ix = first_active,last_active
          
        xx=position(ix,iy,1)-half*box_l  ! Boxlen already in pc
        yy=position(ix,iy,2)-half*box_l_y
        rr=sqrt(yy**2.+Rstar_cu**2)
        force_x(ix,iy)  = 0.0d0
        force_y(ix,iy)  = - m_cloud *Msun/unit_m*yy/rr**3.
        force_z(ix,iy)  = 0.0d0
#if NDUST>0
        do idust=1,ndust
         force_dust_x(idust,ix,iy)  = 0.0d0
         force_dust_y(idust,ix,iy)  =  - m_cloud*Msun/unit_m *yy/rr**3.
         force_dust_z(idust,ix,iy)  = 0.0d0

         if(.not. drag) then
          force_x(ix,iy)  = force_x(ix,iy)  - q(irhod(idust),ix,iy)/q(irho,ix,iy)*(q(ivx,ix,iy) - q(ivdx(idust),ix,iy))/tstop(idust,ix,iy)
          force_y(ix,iy)  = force_y(ix,iy)  - q(irhod(idust),ix,iy)/q(irho,ix,iy)*(q(ivy,ix,iy) - q(ivdy(idust),ix,iy))/tstop(idust,ix,iy)
          force_z(ix,iy)  = force_z(ix,iy)  - q(irhod(idust),ix,iy)/q(irho,ix,iy)*(q(ivz,ix,iy) - q(ivdz(idust),ix,iy))/tstop(idust,ix,iy)
          force_dust_x(idust,ix,iy)  = force_dust_x(idust,ix,iy)  + (q(ivx,ix,iy) - q(ivdx(idust),ix,iy))/tstop(idust,ix,iy)
          force_dust_y(idust,ix,iy)  = force_dust_y(idust,ix,iy)  + (q(ivy,ix,iy) - q(ivdy(idust),ix,iy))/tstop(idust,ix,iy)
          force_dust_z(idust,ix,iy)  = force_dust_z(idust,ix,iy)  + (q(ivz,ix,iy) - q(ivdz(idust),ix,iy))/tstop(idust,ix,iy)
         endif
        end do
#endif
    end do
   end do
end subroutine update_force_setup

#if NDUST>0
! Dust stopping time
subroutine compute_tstop
  
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  real(dp):: pn,rhon
  !Re-calc distribution


   integer :: idust,ix,iy
   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy)
   do iy = first_active_y,last_active_y
    do ix = first_active,last_active
     do idust=1,ndust
        tstop(idust,ix,iy) = sqrt(pi*gamma/8.0d0)*rhograin/unit_d*sdust(idust)/q(irho,ix,iy)/cs(ix,iy)
      end do
     end do
  end do


end subroutine compute_tstop
#endif


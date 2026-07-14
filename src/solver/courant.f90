!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine computes the timestep according to the CFL conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine courant
  use parameters
  use commons
  use units
  implicit none

  integer :: i,idust
  real(dp) :: vmax,dxx,force_max,ca,c_fast,cw,magnetosonic_fast,vv,fratio,D_max

  ! real(dp), dimension(:)  , allocatable :: dt_cw
  ! real(dp), dimension(:)  , allocatable :: dt_diffusion






  if(static)then
     return 
  endif


  

  ! allocate(dt_cw(1:ncells))
  ! dt_cw=0.0d0
  ! allocate(dt_diffusion(1:ncells))
  ! dt_diffusion=0.0d0

  


  dt=2d44


  ca=0.0d0

  do i = 1,ncells
   if(active_cell(i)==1) then 

  
   !Cas 1D   
#if NY==1
   dxx = dx(i,1)
#else
   dxx  = min(dx(i,1),dx(i,2))
#if GEOM==2
   dxx  = min(dx(i,1),radii_c(i)*dx(i,2))
#endif
#endif

   vv   = abs(q(i,ivx)) + abs(q(i,ivy))+ abs(q(i,ivz))
   vmax = cs(i)+vv

#if MHD==1


      ca   = dsqrt((q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/(4*pi*q(i,irho))) !Alfven speed of the dust fluid
      c_fast = dsqrt((cs(i))**2 + ca**2) !Safer to use this one, is the maxmimum speed (get it with theta = pi/2)
      vv   =  abs(q(i,ivx)) + abs(q(i,ivy)) + abs(q(i,ivz)) 

      !magnetosonic_fast = dsqrt(half*(cs(i)**2+(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irho) + dsqrt((cs(i)**2+(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irho))**2-4*cs(i)**2*q(i,iBx)**2/q(i,irho)))) !You may need 4pi factors
      vmax=  max(vmax,c_fast+vv)

#endif

#if NDUST==1     
   idust=1
#if MHD==1
if (ideal_MHD .eqv. .false.) then

      ca   = dsqrt((q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/(4*pi*q(i,irhod(idust)))) !Alfven speed of the dust fluid
      vv   =  abs(q(i,ivdx(idust))) + abs(q(i,ivdy(idust))) + abs(q(i,ivdz(idust))) 
      vmax =  max(vmax,ca+vv)
#if DUST_PRESSURE==1

      c_fast = dsqrt((delta_dust_cs*cs(i))**2 + ca**2) !Safer to use this one

      !magnetosonic_fast = dsqrt(half*(c_fast**2 + dsqrt(c_fast**4-4*(delta_dust_cs*cs(i))**2*ca**2))) !In 1D along B: reduces to a simple soundwave

      vmax =  max(vmax,c_fast+vv)

#endif
endif
#endif


#endif
   dt = min(dt,CFL*dxx/abs(vmax))


#if GRAVITY==1   
      dt = min(dt,CFL*dxx/sqrt(Mc(i)/sqrt(radii_c(i)**2.+(l_soft/unit_l)**2.)))
#endif

if(force_kick) then

   force_max= sqrt(force(i,1)**2+force(i,2)**2+force(i,3)**2)

#if NDUST>0     
   do idust=1,ndust

      force_max= max(force_max,sqrt(force_dust(i,1,idust)**2+force_dust(i,2,idust)**2+force_dust(i,3,idust)**2))

   end do
#endif 
endif


#if MHD==1

if (dusty_nonideal_MHD_no_electron) then !!Adapt timestep to hyper_diffusion in induction equation and Lorentz force (source term)



    if (hyper_diffusion) then
      D_max = max(abs(eta_eff_ohm(i)),abs(eta_eff_Hall_y(i)),abs(eta_eff_Hall_z(i))) !Is necessarily in cgs because resistivities cannot be rendered dimensionless

      dt = min(dt,0.5d0*dxx**2/D_max)

 

   endif

   if (Hall_effect) then

      cw = abs(eta_eff_Hall_y(i))*pi/(2*dxx) + dsqrt((abs(eta_eff_Hall_y(i))*pi/(2*dxx))**2 + ca**2) 

      vmax = max(vmax,cw+vv)

      dt = min(dt,CFL*dxx/abs(vmax))

    endif


endif

if (dusty_nonideal_MHD) then !!Adapt timestep to to magnetocompressive modes, hyper_diffusion in induction equation and Lorentz force (source term)

#if NDUST>1

   do idust=1,ndust

      vv   =  abs(q(i,ivdx(idust))) + abs(q(i,ivdy(idust))) + abs(q(i,ivdz(idust))) 
      vmax =  max(vmax,c_ms_d(i)+vv)

      ! dt_cms(i) = CFL*(dxx/abs(c_ms_d(i)+vv))


      if (Hall_effect) then

         cw = abs(clight**2/(4*pi)*eta_H(i))*pi/(2*dxx) + dsqrt(((clight**2/(4*pi)*eta_H(i))*pi/(2*dxx))**2 + c_ms_d(i)**2)!Which ca should we use?
         vmax = max(vmax,cw+vv)

         ! dt_cw(i) = CFL*(dxx/abs(cw+vv))

      endif

   enddo

   dt = min(dt,CFL*dxx/abs(vmax))

#endif
   

    if (hyper_diffusion_with_electrons) then
      D_max = max(abs(clight**2/(4*pi)*eta_o(i)),abs(2*clight**2/(4*pi)*eta_a(i))) !Is necessarily in cgs because resistivities cannot be rendered dimensionless

      ! print *,'eta_a(i)',eta_a(i)
      !!!D_max can be very high because of eta_AD being high too. This is due to the absence of the dust contribution in the computation of the resistivities!!!
      !!!We have to apply a cap to work with reasonable timestep!!!

      dt = min(dt,0.4d0*dxx**2/D_max)

      ! dt_diffusion(i) = 0.4d0*dxx**2/D_max







 

   endif


endif

#if NDUST>0
if(dusty_nonideal_MHD_no_electron .or. dusty_nonideal_MHD) then


   if (apply_Lorentz_force_explicit) then
      do idust=1,ndust
         if (FLor_x_d(i,idust) /= 0.0d0 .or. FLor_y_d(i,idust) /= 0.0d0 .or. FLor_z_d(i,idust) /= 0.0d0) then
            dt = min(dt,CFL*dsqrt(dxx/dsqrt(FLor_x_d(i,idust)**2+FLor_y_d(i,idust)**2+FLor_z_d(i,idust)**2)/q(i,irhod(idust))))

         endif
      end do
   endif
end if
#endif

if(dusty_nonideal_MHD_no_electron .or. dusty_nonideal_MHD) then

   if (apply_Lorentz_force_explicit) then
      if (FLor_x(i) /= 0.0d0 .or. FLor_y(i) /= 0.0d0 .or. FLor_z(i) /= 0.0d0) then
         dt=min(dt,CFL*dsqrt(dxx/dsqrt(FLor_x(i)**2+FLor_y(i)**2+FLor_z(i)**2)/q(i,irho)))

      endif
   endif
endif

#endif

   if(vv.ne.0.0d0) then
      fratio = max(force_max*dxx/vv**2,1d-3)
      dt = min(dt,CFL*dxx/vv*(sqrt(1.0d0+2.0d0*CFL*fratio)-1.0d0)/fratio)      !print*, 'dt dfratio',dt

   endif

   endif
  end do


 ! print*, 'dt final',dt

 ! CFL_deactivate_diffusion = .false.

 ! if (dt > minval(dt_diffusion)) CFL_deactivate_diffusion = .true.



! deallocate(dt_cw)
! deallocate(dt_diffusion)


#if NY>1
  !print *, 'time = ', time, 'dt = ', dt
#endif 
end subroutine courant








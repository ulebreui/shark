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
  real(dp) :: vmax,dxx,force_max,ca,magnetosonic_fast,vv,fratio,D_max


  if(static)then
     return 
  endif


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
#if NDUST==0
      magnetosonic_fast = dsqrt(half*(cs(i)**2+(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irho) + dsqrt((cs(i)**2+(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irho))**2-4*cs(i)**2*q(i,iBx)**2/q(i,irho))))
      vmax=  max(vmax,magnetosonic_fast+abs(q(i,ivx)))
#endif
#endif

#if NDUST>0     
   do idust=1,ndust
#if MHD==1
      ca   = dsqrt((q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irhod(idust))) !TODO , to modify when accounting for a dust distribution
#endif
      vv   =  abs(q(i,ivdx(idust))) + abs(q(i,ivdy(idust))) + abs(q(i,ivdz(idust)))
      vmax =  max(vmax,ca+vv)
   enddo

#endif
   !print(vmax)
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

#if MHD==1
if (dusty_nonideal_MHD) then !!Adapt timestep to hyper_diffusion in induction equation

   call effective_diffusion_coef_induction
   D_max = max(eta_eff_yy(i),eta_eff_yz(i),eta_eff_zy(i),eta_eff_zz(i)) !Is necessarily in cgs because resistivities cannot be rendered dimensionless
   print *, D_max
   dt = min(dt,CFL*dxx**2/D_max)

end if

#endif

   if(vv.ne.0.0d0) then
      fratio = max(force_max*dxx/vv**2,1d-3)
      dt = min(dt,CFL*dxx/vv*(sqrt(1.0d0+2.0d0*CFL*fratio)-1.0d0)/fratio)
   endif
endif
   endif
  end do
#if NY>1
  !print *, 'time = ', time, 'dt = ', dt
#endif 
end subroutine courant








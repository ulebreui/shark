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



#if NDUST>0     
   do idust=1,ndust
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








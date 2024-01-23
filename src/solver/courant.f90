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
  real(dp) :: vmax,dxx,force_max,ca,magnetosonic_fast
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
   vmax = cs(i)+abs(q(i,ivx))
#if NY==1
   vmax = max(vmax,abs(q(i,ivy)))!cs propagates only in x here because 1D
#endif
#if NY>1
   vmax = max(vmax,abs(q(i,ivy))+cs(i))
#endif
   vmax = max(vmax,abs(q(i,ivz)))!cs propagates only in x here because 1D

#if MHD==1
#if NDUST==0
      magnetosonic_fast = dsqrt(half*(cs(i)**2+(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irho) + dsqrt((cs(i)**2+(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irho))**2-4*cs(i)**2*q(i,iBx)**2/q(i,irho))))
      vmax=  max(vmax,magnetosonic_fast+abs(q(i,ivx)))


#endif
#endif

#if NDUST>0     
   do idust=1,ndust
#if MHD==1
      ca = dsqrt((q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irhod(idust))) !TODO , to modify when accounting for a dust distribution
#endif
      vmax=  max(vmax,ca+abs(q(i,ivdx(idust))))
      vmax=  max(vmax,abs(q(i,ivdy(idust))))
      vmax=  max(vmax,abs(q(i,ivdz(idust))))

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

   if(force_max>0.d0)dt = min(dt,CFL*sqrt(dxx/force_max))
endif
   endif
  end do
#if NY>1
  !print *, 'time = ', time, 'dt = ', dt
#endif 
end subroutine courant

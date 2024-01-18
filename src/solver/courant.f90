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
  use OMP_LIB

  implicit none

  integer  :: i,idust
  real(dp) :: vmax,dxx,force_max,ca,magnetosonic_fast,vmax_x,vmax_y
  if(static)then
     return
  endif
  ca=0.0d0
  dtcells=1d10

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,vmax,dxx,force_max,ca,magnetosonic_fast,vmax_x,vmax_y)
  !$OMP DO
  do i = 1,ncells
   if(active_cell(i)) then   

   vmax_x = cs(i)+abs(q(i,ivx))
#if NY>1
   vmax_y = abs(q(i,ivy))+cs(i)
#endif

#if MHD==1
#if NDUST==0
      magnetosonic_fast   = dsqrt(half*(cs(i)**2+(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irho) + dsqrt((cs(i)**2+(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irho))**2-4*cs(i)**2*q(i,iBx)**2/q(i,irho))))
      vmax_x              = max(vmax_x,magnetosonic_fast+abs(q(i,ivx)))
#endif
#endif

#if NDUST>0     
   do idust=1,ndust
#if MHD==1
      ca     = dsqrt((q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)/q(i,irhod(idust))) !TODO , to modify when accounting for a dust distribution
#endif
      vmax_x =  max(vmax_x,ca+abs(q(i,ivdx(idust))))
#if NY>1
      vmax_y =  max(vmax_y,abs(q(i,ivdy(idust))))
#endif
   enddo

#endif

      dtcells(i) = min(dtcells(i) , CFL*dx(i,1)/abs(vmax_x))
#if NY>1
#if GEOM==0
      dtcells(i) = min(dtcells(i) , CFL*dx(i,2)/abs(vmax_y))
#endif
#if GEOM==2
      dtcells(i) = min(dtcells(i) , CFL*radii_c(i)*dx(i,2)/abs(vmax_y))
#endif
#endif

if(force_kick) then
#if NY==1
   dxx  = dx(i,1)
#else
   dxx  = min(dx(i,1),dx(i,2))
#if GEOM==2
   dxx  = min(dx(i,1),radii_c(i)*dx(i,2))
#endif
#endif

   force_max = sqrt(force(i,1)**2+force(i,2)**2+force(i,3)**2)

#if NDUST>0     
   do idust=1,ndust

      force_max = max(force_max,sqrt(force_dust(i,1,idust)**2+force_dust(i,2,idust)**2+force_dust(i,3,idust)**2))

   end do
#endif 

   if(force_max>0.d0) dtcells(i)  = min(dtcells(i) ,CFL*sqrt(dxx/force_max))
endif


#if GRAVITY==1   
      dtcells(i)  = min(dtcells(i) ,CFL*dx(i,1)/sqrt(Mc(i)/sqrt(radii_c(i)**2.+(l_soft/unit_l)**2.)))
#endif
   endif


  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  dt = MINVAL(dtcells)

#if NY>1
  !print *, 'time = ', time, 'dt = ', dt
#endif 
end subroutine courant

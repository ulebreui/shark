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

  integer :: i,idust,ix,iy,icell
  real(dp) :: vmax,dxx,force_max,ca,magnetosonic_fast,vv,fratio,D_max


  if(static)then
     return 
  endif


  dt=2d44
  ca=0.0d0

  do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         i = icell(ix,iy)
         !Cas 1D   
#if NY==1
         dxx = dx(ix,iy,1)
#else
         dxx  = min(dx(ix,iy,1),dx(ix,iy,2))
#if GEOM==2
         dxx  = min(dx(ix,iy,1),radii_c(i)*dx(ix,iy,2))
#endif
#endif

         vv   = abs(q(ix,iy,ivx)) + abs(q(ix,iy,ivy))+ abs(q(ix,iy,ivz))
         vmax = cs(i)+vv



#if NDUST>0     
         do idust=1,ndust
            vv   =  abs(q(ix,iy,ivdx(idust))) + abs(q(ix,iy,ivdy(idust))) + abs(q(ix,iy,ivdz(idust)))
            vmax =  max(vmax,ca+vv)
         enddo
#endif
         !print(vmax)
         dt = min(dt,CFL*dxx/abs(vmax))
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
      end do
  end do
end subroutine courant








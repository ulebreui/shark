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
  real(dp) :: vmax,dxx,force_max
  if(static)then
     return
  endif
  dt=2d44


  do i = 1,ncells
   if(active_cell(i)==1) then   
   !Cas 1D   
#if NY==1
   vmax= abs(q(i,ivx))
   dxx = dx(i,1)
#endif


!Cas 2D/2.5D
#if NY>1
   dxx  = min(dx(i,1),dx(i,2))
   !vmax = sqrt(q(i,ivx)**2+q(i,ivy)**2+q(i,ivz)**2)
   vmax = sqrt(q(i,ivx)**2+q(i,ivy)**2)

#endif
   dt = min(dt,CFL*dxx/(cs(i)+abs(vmax))) 

#if NDUST>0     
   do idust=1,ndust
#if NY==1
      vmax= abs(q(i,ivdx(idust)))
#endif
#if NY>1
      vmax=  sqrt(q(i,ivdx(idust))**2+q(i,ivdy(idust))**2)
      !vmax=  sqrt(q(i,ivdx(idust))**2+q(i,ivdy(idust))**2+q(i,ivdz(idust))**2)

#endif
      dt = min(dt,CFL*dxx/abs(vmax))
   enddo
#endif
! if(force_kick) then
! #if NY==1
!    force_max = abs(force(i,1))
! #endif
! #if NY>1
!    force_max= sqrt(force(i,1)**2+force(i,2)**2+force(i,3)**2)
! #endif

! #if NDUST>0     
!    do idust=1,ndust
! #if NY==1
!       force_max= max(force_max,sqrt(force_dust(i,1,idust)**2))
! #endif
! #if NY>1
!       force_max= max(force_max,sqrt(force_dust(i,1,idust)**2+force_dust(i,2,idust)**2+force_dust(i,3,idust)**2))

! #endif
!    end do
! #endif 

!    if(force_max>0.d0)dt = min(dt,CFL*sqrt(dxx/force_max))
! endif
   endif
  end do

!print *, 'time = ', time, 'dt = ', dt

end subroutine courant

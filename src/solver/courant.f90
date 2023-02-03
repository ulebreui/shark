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
  real(dp) :: vmax,dxx
  if(static)then
     return
  endif
  dt=2d44

  do i = 1,ncells
   if(active_cell(i)==1) then
      if(ndim==1)vmax = abs(q(i,iv))
      if(ndim==2)vmax = sqrt(q(i,iv)**2+q(i,ivy)**2)
      dxx=dx(i,1)
      if(ndim==2)dxx = min(dx(i,1),dx(i,2))
      dt = min(dt,CFL*dxx/(sqrt(gamma*q(i,iP)/q(i,irho))+abs(vmax)))    
#if NDUST>0     
     do idust=1,ndust
         if(ndim==1)vmax = abs(q(i,ivd(idust)))
         if(ndim==2)vmax = sqrt(q(i,ivd(idust))**2+q(i,ivdy(idust))**2)
         dxx=dx(i,1)
         if(ndim==2)dxx = min(dx(i,1),dx(i,2))
         dt = min(dt,CFL*dxx/abs(vmax))
     end do
#endif     
   endif
  end do

  print *, 'time = ', time, 'dt = ', dt

end subroutine courant

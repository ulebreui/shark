! This routine performs the predictor operation ! 
subroutine predictor
  use parameters
  use commons
  use units
  use OMP_LIB
  implicit none
  integer :: i,ix,iy,icell,ivv,iuu,iymin,iymax,idust
  integer :: ixx,iyy
  integer :: ivar,idim,ix0,iy0
  real(dp) :: slope_left,slope_right,slope_lim
  real(dp), dimension(:,:,:),allocatable :: dq
  real(dp), dimension(:,:)  ,allocatable :: qpred
  if(static) return

  iymin=1
  iymax=1
  if(ndim==2) then
   iymin=2
   iymax=ny_max-1
  endif

  ! Initialise to zero
  allocate(dq(1:ncells,1:nvar,1:ndim))
  allocate(qpred(1:ncells,1:nvar))

  dq    = 0.0d0
  qp    = 0.0d0
  qm    = 0.0d0

  call apply_boundaries(1,uold,ncells,nvar)

  ! Computes primitive variables
  call ctoprim
  qpred = q
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,ivar,ix,iy,idim,ix0,iy0,slope_left,slope_right,ivv,iuu,idust)


  !$OMP DO
  do i=1,ncells
    if(active_cell_predictor(i)==1) then
      ix=ixx(i)
      iy=iyy(i)
      do ivar = 1,nvar
        slope_left  = 2.0d0*(q(icell(ix,iy),ivar) - q(icell(ix-1,iy),ivar))/(dx(icell(ix,iy),1)+dx(icell(ix-1,iy),1))
        slope_right = 2.0d0*(q(icell(ix+1,iy),ivar) - q(icell(ix,iy),ivar))/(dx(icell(ix+1,iy),1)+dx(icell(ix,iy),1))
        dq(icell(ix,iy),ivar,1) = slope_left
        if(abs(slope_right)<abs(slope_left)) dq(icell(ix,iy),ivar,1) = slope_right
        if(slope_right*slope_left<0.0d0) dq(icell(ix,iy),ivar,1) = 0.0d0
#if NY>1
        slope_left  = 2.0d0*(q(icell(ix,iy),ivar) - q(icell(ix,iy-1),ivar))/(dx(icell(ix,iy),2)+dx(icell(ix,iy-1),2))
        slope_right = 2.0d0*(q(icell(ix,iy+1),ivar) - q(icell(ix,iy),ivar))/(dx(icell(ix,iy+1),2)+dx(icell(ix,iy),2))
        dq(icell(ix,iy),ivar,2) = slope_left
        if(abs(slope_right)<abs(slope_left)) dq(icell(ix,iy),ivar,2) = slope_right
        if(slope_right*slope_left<0.0d0) dq(icell(ix,iy),ivar,2) = 0.0d0
#endif            
      end do

      !First, we take care of the 1D terms
      qpred(i,irho) = qpred(i,irho) + half*dt*(-dq(i,irho,1)*q(i,iv)-dq(i,iv,1)*q(i,irho))
      qpred(i,iv)   = qpred(i,iv)   + half*dt*(-dq(i,iv,1)*q(i,iv)-dq(i,iP,1)/q(i,irho))
      qpred(i,iP)   = qpred(i,iP)   + half*dt*(-q(i,iv)*dq(i,iP,1)-gamma*q(i,iP)*dq(i,iv,1))

      !We now add the terms that are specific to 2D problems
#if NY>1
        qpred(i,irho) = qpred(i,irho) + half*dt*(-dq(i,irho,2)*q(i,ivy)-dq(i,ivy,2)*q(i,irho))
        qpred(i,iv)   = qpred(i,iv)   + half*dt*(-dq(i,iv,2)*q(i,ivy))
        qpred(i,ivy)  = qpred(i,ivy)  + half*dt*(-dq(i,ivy,2)*q(i,ivy)-dq(i,ivy,1)*q(i,iv)-dq(i,iP,2)/q(i,irho))
        qpred(i,iP)   = qpred(i,iP)   + half*dt*(-q(i,ivy)*dq(i,iP,2)-gamma*q(i,iP)*dq(i,ivy,2))
#endif        
      !Dust terms
#if NDUST>0
      do idust=1,ndust
        qpred(i,irhod(idust)) = qpred(i,irhod(idust)) + half*dt*(-dq(i,irhod(idust),1)*q(i,ivd(idust))-dq(i,ivd(idust),1)*q(i,irhod(idust)))
        qpred(i,ivd(idust))   = qpred(i,ivd(idust))   + half*dt*(-dq(i,ivd(idust),1)*q(i,ivd(idust)))
#if NY>1
        qpred(i,irhod(idust)) = qpred(i,irhod(idust)) + half*dt*(-dq(i,irhod(idust),2)*q(i,ivdy(idust))-dq(i,ivdy(idust),2)*q(i,irhod(idust)))
        qpred(i,ivd(idust))   = qpred(i,ivd(idust))   + half*dt*(-dq(i,ivd(idust),2)*q(i,ivdy(idust)))
        qpred(i,ivdy(idust))  = qpred(i,ivdy(idust))  + half*dt*(-dq(i,ivdy(idust),2)*q(i,ivdy(idust))-dq(i,ivdy(idust),1)*q(i,ivd(idust)))
#endif          
      end do
#endif
    end if
  end do
  !$OMP END DO

  !$OMP BARRIER

  !$OMP DO
      do i=1,ncells
        do ivar = 1,nvar
          qm(i,ivar,1)=qpred(i,ivar)-half*dq(i,ivar,1)*dx(i,1)
          qp(i,ivar,1)=qpred(i,ivar)+half*dq(i,ivar,1)*dx(i,1)
#if NY>1    
          qm(i,ivar,2)=qpred(i,ivar)-half*dq(i,ivar,2)*dx(i,2)        
          qp(i,ivar,2)=qpred(i,ivar)+half*dq(i,ivar,2)*dx(i,2)
#endif          
      end do
    end do
  !$OMP END DO
  !$OMP END PARALLEL
  deallocate(dq)
  deallocate(qpred)
end subroutine predictor

! This is the Riemmann solver : llf + Godunov scheme
subroutine add_delta_u
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust,ivar,ix,iy,il,ily,icell,idim,ivn,ivt
  integer :: ixx,iyy

  real(dp), dimension(:,:)  , allocatable  :: delta_U
  real(dp), dimension(:,:,:), allocatable  :: flux
  real(dp), dimension(1:nvar) :: flux_left,qleft,flux_right,qright,uleft,uright,lambda_llf,ustar,qstar,fstar,qstarleft,qstarright

  real(dp) :: S_left,S_right,hllc_l,hllc_r,r_o,u_o,P_o,e_o
  
  if(static) return

  !Initialisation of the flux,lambda_llf and delta U: they must be allocatable for 2D simus with h-res


  allocate(delta_U(1:ncells,1:nvar))
  allocate(flux(1:ncells,1:nvar,1:ndim))
  flux       = 0.0d0
  delta_U    = 0.0d0

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&

  !$OMP PRIVATE(i,idust,ivar,ix,iy,il,ily,idim,ivn,ivt,lambda_llf,flux_left,flux_right,uleft,uright,qleft,qright,S_right,S_left,hllc_r,hllc_l,qstarleft,qstarright,r_o,u_o,P_o,e_o,qstar)

  !$OMP DO
  do i = 1, ncells
    if(active_cell_predictor(i)==1) then
      ix=ixx(i)
      iy=iyy(i)
      do idim = 1,ndim 
        il = icell(ix+1,iy)
        ivn = iv
        ivt = ivy
        if(idim==2) then
          il = icell(ix,iy+1)
          ivn = ivy
          ivt = iv
        endif 

        flux_left  = 0.0d0
        flux_right = 0.0d0
        uleft      = 0.0d0
        uright     = 0.0d0
        qleft      = 0.0d0
        qright     = 0.0d0
        lambda_llf = 0.0d0
        S_right    = 0.0d0
        S_left     = 0.0d0
#if SOLVER==2
          qstar      = 0.0d0
          qstarleft  = 0.0d0
          qstarright = 0.0d0
#endif
        
        !Density
        qright(irho)  = qm(il,irho,idim) !rho
        qleft(irho)   = qp(i,irho,idim)

        uright(irho)  = qm(il,irho,idim)
        uleft(irho)   = qp(i,irho,idim)

        flux_right(irho) = qm(il,irho,idim) * qm(il,ivn,idim) ! rho u
        flux_left(irho)  = qp(i,irho,idim)  * qp(i,ivn,idim)

        !Momentum
        qright(ivn)   = qm(il,ivn,idim) ! v
        qleft(ivn)    = qp(i,ivn,idim)

        uright(ivn)   = qm(il,irho,idim) * qm(il,ivn,idim)  ! rho u
        uleft(ivn)    = qp(i,irho,idim)  * qp(i,ivn,idim)  
        
        flux_right(ivn) = qm(il,irho,idim) * qm(il,ivn,idim) * qm(il,ivn,idim) + qm(il,iP,idim) ! rho u u + P
        flux_left(ivn)  = qp(i,irho,idim)  * qp(i,ivn,idim)  * qp(i,ivn,idim)  + qp(i,iP,idim)


#if NY>1
          !Transverse momentum
          qright(ivt)     = qm(il,ivt,idim) ! v
          qleft(ivt)      = qp(i,ivt,idim)

          uright(ivt)     = qm(il,irho,idim) * qm(il,ivt,idim) ! rho v
          uleft(ivt)      = qp(i,irho,idim)  * qp(i,ivt,idim) 
        
          flux_right(ivt) = qm(il,irho,idim) * qm(il,ivn,idim) * qm(il,ivt,idim) ! rho u v
          flux_left(ivt)  = qp(i,irho,idim)  * qp(i,ivn,idim)  * qp(i,ivt,idim) 
#endif

        !Energy

        qright(iP)    = qm(il,iP,idim) ! P
        qleft(iP)     = qp(i,iP,idim)

        uright(iP)    = qm(il,iP,idim)/(gamma-1.d0) + half * qm(il,irho,idim) * qm(il,ivn,idim) * qm(il,ivn,idim)! E
        uleft(iP)     = qp(i,iP,idim) /(gamma-1.d0) + half * qp(i,irho,idim)  * qp(i,ivn,idim)  *  qp(i,ivn,idim)
#if NY>1
        uright(iP)   = uright(iP)  + half * qm(il,irho,idim) * qm(il,ivt,idim) * qm(il,ivt,idim)
        uleft(iP)    = uleft(iP)   + half * qp(i,irho,idim)  * qp(i,ivt,idim)  * qp(i,ivt,idim)
#endif
        flux_right(iP) = (uright(iP)+qright(iP)) * qright(ivn) ! (E+P) v
        flux_left(iP)  = (uleft(iP) +qleft(iP))  * qleft(ivn)

#if NDUST>0
        do idust=1,ndust
          ivn = ivd(idust)
          ivt = ivdy(idust)
          if(idim==2) then
            ivn = ivdy(idust)
            ivt = ivd(idust)
          endif 
          !Dust density
          qright(irhod(idust))  = qm(il,irhod(idust),idim)
          qleft(irhod(idust))   = qp(i,irhod(idust),idim)

          uright(irhod(idust))  = qm(il,irhod(idust),idim)
          uleft(irhod(idust))   = qp(i,irhod(idust),idim)

          flux_right(irhod(idust))=qm(il,irhod(idust),idim) * qm(il,ivn,idim)
          flux_left(irhod(idust)) =qp(i,irhod(idust),idim)  * qp(i,ivn,idim)

          !Dust momentum
          qright(ivn)    = qm(il,ivn,idim)
          qleft(ivn)     = qp(i,ivn,idim)

          uright(ivn)    = qm(il,irhod(idust),idim) * qm(il,ivn,idim)
          uleft(ivn)     = qp(i,irhod(idust),idim)  * qp(i,ivn,idim)

          flux_right(ivn)= qm(il,irhod(idust),idim) * qm(il,ivn,idim)**2
          flux_left(ivn) = qp(i,irhod(idust),idim)  * qp(i,ivn,idim)**2
#if NY>1
          !Dust transverse momentum
            qright(ivt)    = qm(il,ivt,idim)
            qleft(ivt)     = qp(i,ivt,idim)

            uright(ivt)    = qm(il,ivt,idim)* qm(il,irhod(idust),idim)
            uleft(ivt)     = qp(i,ivt,idim) * qp(i,irhod(idust),idim)

            flux_right(ivt)= qm(il,irhod(idust),idim)*qm(il,ivn,idim)*qm(il,ivt,idim)
            flux_left(ivt) = qp(i,irhod(idust),idim)*qp(i,ivn,idim)*qp(i,ivt,idim)
#endif
        end do
#endif

!LLF
#if SOLVER==0 
        ivn = iv
        ivt = ivy
        if(idim==2) then
          ivn = ivy
          ivt = iv
        endif 
          lambda_llf = max(abs(qleft(ivn))+sqrt(gamma*qleft(iP)/qleft(irho)),abs(qright(ivn))+sqrt(gamma*qright(iP)/qright(irho)))
#if NDUST>0
          do idust=1,ndust
            ivn = ivd(idust)
            if(idim==2) then
              ivn = ivdy(idust)
            endif 
            lambda_llf(irhod(idust)) = max(abs(qleft(ivn)),abs(qright(ivn)))
            lambda_llf(ivd(idust))   = max(abs(qleft(ivn)),abs(qright(ivn)))
            if(ndim==2)lambda_llf(ivdy(idust)) = max(abs(qleft(ivn)),abs(qright(ivn)))
          end do
#endif  
        do ivar=1,nvar
              flux(i,ivar,idim)  = half  * (flux_left(ivar)+flux_right(ivar))-half*lambda_llf(ivar)* (uright(ivar)-uleft(ivar))
        end do
#endif

!HLL
#if SOLVER==1
          ivn = iv
          ivt = ivy
          if(idim==2) then
            ivn = ivy
            ivt = iv
          endif 
          S_left  = min(min(qleft(ivn),qright(ivn))-max(sqrt(gamma*qleft(iP)/qleft(irho)),sqrt(gamma*qright(iP)/qright(irho))),0.0d0)
          S_right = max(max(qleft(ivn),qright(ivn))+max(sqrt(gamma*qleft(iP)/qleft(irho)),sqrt(gamma*qright(iP)/qright(irho))),0.0d0)
#if NDUST>0
          do idust=1,ndust
            ivn = ivd(idust)
            if(idim==2) then
              ivn = ivdy(idust)
            endif 
            lambda_llf(irhod(idust)) = max(abs(qleft(ivn)),abs(qright(ivn)))
            lambda_llf(ivd(idust))   = max(abs(qleft(ivn)),abs(qright(ivn)))
            if(ndim==2)lambda_llf(ivdy(idust)) = max(abs(qleft(ivn)),abs(qright(ivn)))            
          end do
#endif  
          do ivar = 1,iP
              flux(i,ivar,idim)  = (S_right*flux_left(ivar)-S_left*flux_right(ivar)+S_right*S_left*(uright(ivar)-uleft(ivar)))/(S_right-S_left)
          end do

#if NDUST>0
          do idust=1,ndust
            flux(i,irhod(idust),idim)  = half  * (flux_left(irhod(idust))+flux_right(irhod(idust)))-half*lambda_llf(irhod(idust))* (uright(irhod(idust))-uleft(irhod(idust)))
            flux(i,ivd(idust),idim)    = half  * (flux_left(ivd(idust))+flux_right(ivd(idust)))-half*lambda_llf(ivd(idust))* (uright(ivd(idust))-uleft(ivd(idust)))
            if(ndim==2)flux(i,ivdy(idust),idim)   = half  * (flux_left(ivdy(idust))+flux_right(ivdy(idust)))-half*lambda_llf(ivdy(idust))* (uright(ivdy(idust))-uleft(ivdy(idust)))
        enddo
#endif            
#endif

!HLLC
#if SOLVER==2
          ivn = iv
          ivt = ivy
          if(idim==2) then
            ivn = ivy
            ivt = iv
          endif 
          S_left  = min(min(qleft(ivn),qright(ivn))-max(sqrt(gamma*qleft(iP)/qleft(irho)),sqrt(gamma*qright(iP)/qright(irho))),0.0d0)
          S_right = max(max(qleft(ivn),qright(ivn))+max(sqrt(gamma*qleft(iP)/qleft(irho)),sqrt(gamma*qright(iP)/qright(irho))),0.0d0) 

          ! Compute lagrangian sound speed
          hllc_l = qleft(irho)*(qleft(ivn)-S_left)
          hllc_r = qright(irho)*(S_right-qright(ivn))

          qstar(ivn) = (hllc_r*qright(ivn)   +hllc_l*qleft(ivn)   +  (qleft(iP)-qright(iP)))/(hllc_r+hllc_l)
          qstar(iP)  = (hllc_r*qleft(iP)+hllc_l*qright(iP)+hllc_l*hllc_r*(qleft(ivn)-qright(ivn)))/(hllc_r+hllc_l)
          ! Left star region variables
          qstarleft(irho)=qleft(irho)*(S_left-qleft(ivn))/(S_left-qstar(ivn))
          qstarleft(iP)=((S_left-qleft(ivn))*uleft(iP)-qleft(iP)*qleft(ivn)+qstar(iP)*qstar(ivn))/(S_left-qstar(ivn))
          ! Right star region variables
          qstarright(irho)=qright(irho)*(S_right-qright(ivn))/(S_right-qstar(ivn))
          qstarright(iP)=((S_right-qright(ivn))*uright(iP)-qright(iP)*qright(ivn)+qstar(iP)*qstar(ivn))/(S_right-qstar(ivn))

        ! Sample the solution at x/t=0
        if(S_left>0.0d0)then
          r_o=qleft(irho)
          u_o=qleft(ivn)
          P_o=qleft(iP)
          e_o=uleft(iP)
        else if(qstar(ivn)>0.0d0)then
          r_o=qstarleft(irho)
          u_o=qstar(ivn)
          P_o=qstar(iP)
          e_o=qstarleft(iP)
        else if (S_right>0d0)then
          r_o=qstarright(irho)
          u_o=qstar(ivn)
          P_o=qstar(iP)
          e_o=qstarright(iP)
        else
          r_o=qright(irho)
          u_o=qright(ivn)
          P_o=qright(iP)
          e_o=uright(iP)
          !eo=er
        end if
        flux(i,irho,idim) = r_o*u_o
        flux(i,ivn,idim)  = r_o*u_o*u_o+P_o
        flux(i,iP,idim)   = (e_o+P_o)*u_o
        if(ndim==2) then
          if(qstar(ivn)>0.0d0) then
            flux(i,ivt,idim)  = r_o*u_o*qleft(ivt)
          else
            flux(i,ivt,idim)  = r_o*u_o*qright(ivt)
          endif
        endif
#if NDUST>0
          do idust=1,ndust
            ivn = ivd(idust)
            if(idim==2) then
              ivn = ivdy(idust)
            endif 
            lambda_llf(irhod(idust)) = max(abs(qleft(ivn)),abs(qright(ivn)))
            lambda_llf(ivd(idust))   = max(abs(qleft(ivn)),abs(qright(ivn)))
            if(ndim==2)lambda_llf(ivdy(idust)) = max(abs(qleft(ivn)),abs(qright(ivn)))            
          end do
          do idust=1,ndust
            flux(i,irhod(idust),idim)  = half  * (flux_left(irhod(idust))+flux_right(irhod(idust)))- half*lambda_llf(irhod(idust))* (uright(irhod(idust))-uleft(irhod(idust)))
            flux(i,ivd(idust),idim)    = half  * (flux_left(ivd(idust))+flux_right(ivd(idust)))    - half*lambda_llf(ivd(idust))* (uright(ivd(idust))-uleft(ivd(idust)))
            if(ndim==2)flux(i,ivdy(idust),idim)   = half  * (flux_left(ivdy(idust))+flux_right(ivdy(idust)))-half*lambda_llf(ivdy(idust))* (uright(ivdy(idust))-uleft(ivdy(idust)))
        enddo
#endif 
#endif

      end do
    end if
  end do
  !$OMP END DO

  !$OMP BARRIER

  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
        ix = ixx(i)
        iy = iyy(i)
        do ivar = 1, nvar
            il = icell(ix-1,iy)
            delta_U(i,ivar)=(flux(il,ivar,1)*surf(il,1)-flux(i,ivar,1)*surf(i,1))/vol(i)*dt
            if(ndim==2) then
              il = icell(ix,iy-1)
              delta_U(i,ivar)=delta_U(i,ivar)+(flux(il,ivar,2)*surf(il,2)-flux(i,ivar,2)*surf(i,2))/vol(i)*dt
            endif
        end do
      endif
    end do
  !end do
  !$OMP END DO
  !$OMP END PARALLEL
  !Update state vector 
  unew=unew+delta_U


  deallocate(delta_U)
  deallocate(flux)
end subroutine add_delta_u

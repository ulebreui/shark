! This routine performs the predictor operation ! 
subroutine predictor
  use parameters
  use commons
  use units
  use OMP_LIB
  implicit none
  integer :: i,ix,iy,icell,ivv,iuu,iymin,iymax,idust
  integer :: ivar,idim,ix0,iy0
  real(dp) :: slope_left,slope_right,minmod
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
  !$OMP PRIVATE(i,ivar,ix,iy,idim,ix0,iy0,slope_left,slope_right,minmod,ivv,iuu,idust)
  !$OMP DO
    do ix = 2, nx_max-1
    do iy = iymin, iymax
      do ivar = 1,nvar
        do idim=1,ndim
          ix0=0
          iy0=0
          if(idim==1) ix0=1
          if(idim==2) iy0=1
          slope_left  = 2.0d0*(q(icell(ix,iy),ivar) - q(icell(ix-ix0,iy-iy0),ivar))/(dx(icell(ix,iy),idim)+dx(icell(ix-ix0,iy-iy0),idim))
          slope_right = 2.0d0*(q(icell(ix+ix0,iy+iy0),ivar) - q(icell(ix,iy),ivar))/(dx(icell(ix+ix0,iy+iy0),idim)+dx(icell(ix,iy),idim))
          minmod = slope_left
          if(abs(slope_right)<abs(slope_left)) minmod = slope_right
          if(slope_right*slope_left<0.0d0) minmod = 0.0d0
          dq(icell(ix,iy),ivar,idim)=minmod  
        end do
      end do 
    end do
  end do
  !$OMP END DO
  
  !$OMP BARRIER 

  !$OMP DO
  do ix = 2, nx_max-1
    do iy = iymin,iymax
      i=icell(ix,iy)
      !print*, ix,iy
      !First, we take care of the 1D terms
      qpred(i,irho) = qpred(i,irho) + half*dt*(-dq(i,irho,1)*q(i,iv)-dq(i,iv,1)*q(i,irho))
      qpred(i,iv)   = qpred(i,iv)   + half*dt*(-dq(i,iv,1)*q(i,iv)-dq(i,iP,1)/q(i,irho))
      qpred(i,iP)   = qpred(i,iP)   + half*dt*(-q(i,iv)*dq(i,iP,1)-gamma*q(i,iP)*dq(i,iv,1))

      !We now add the terms that are specific to 2D problems
      if(ndim==2) then
        qpred(i,irho) = qpred(i,irho) + half*dt*(-dq(i,irho,2)*q(i,ivy)-dq(i,ivy,2)*q(i,irho))
        qpred(i,iv)   = qpred(i,iv)   + half*dt*(-dq(i,iv,2)*q(i,ivy))
        qpred(i,ivy)  = qpred(i,ivy)  + half*dt*(-dq(i,ivy,2)*q(i,ivy)-dq(i,ivy,1)*q(i,iv)-dq(i,iP,2)/q(i,irho))
        qpred(i,iP)   = qpred(i,iP)   + half*dt*(-q(i,ivy)*dq(i,iP,2)-gamma*q(i,iP)*dq(i,ivy,2))
      endif
      !Dust terms
#if NDUST>0
      do idust=1,ndust
        qpred(i,irhod(idust)) = qpred(i,irhod(idust)) + half*dt*(-dq(i,irhod(idust),1)*q(i,ivd(idust))-dq(i,ivd(idust),1)*q(i,irhod(idust)))
        qpred(i,ivd(idust))   = qpred(i,ivd(idust))   + half*dt*(-dq(i,ivd(idust),1)*q(i,ivd(idust)))
        if(ndim==2) then
          qpred(i,irhod(idust)) = qpred(i,irhod(idust)) + half*dt*(-dq(i,irhod(idust),2)*q(i,ivdy(idust))-dq(i,ivdy(idust),2)*q(i,irhod(idust)))
          qpred(i,ivd(idust))   = qpred(i,ivd(idust))   + half*dt*(-dq(i,ivd(idust),2)*q(i,ivdy(idust)))
          qpred(i,ivdy(idust))  = qpred(i,ivdy(idust))  + half*dt*(-dq(i,ivdy(idust),2)*q(i,ivdy(idust))-dq(i,ivdy(idust),1)*q(i,ivd(idust)))
        endif
      end do
#endif
    end do
  end do
  !$OMP END DO

  !$OMP BARRIER

  !$OMP DO
      do ix=1,nx_max
        do iy=1,ny_max
          do idim =1,ndim
           do ivar =1,nvar
            i = icell(ix,iy)
            qm(i,ivar,idim)=qpred(i,ivar)-half*dq(i,ivar,idim)*dx(i,idim)
            qp(i,ivar,idim)=qpred(i,ivar)+half*dq(i,ivar,idim)*dx(i,idim)
        end do
      end do
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
  real(dp), dimension(:,:)  , allocatable  :: lambda_llf
  real(dp), dimension(:,:)  , allocatable  :: delta_U
  real(dp), dimension(:,:,:), allocatable  :: flux
  real(dp) :: eleft, eright, vleft, vright, vtleft, vtright, rholeft,rhoright,Pleft,Pright
  if(static) return

  !Initialisation of the flux,lambda_llf and delta U: they must be allocatable for 2D simus with h-res
  allocate(lambda_llf(1:ncells,1:ndim))
  allocate(delta_U(1:ncells,1:nvar))
  allocate(flux(1:ncells,1:nvar,1:ndim))
  flux       = 0.0d0
  delta_U    = 0.0d0
  lambda_llf = 0.0d0
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&

  !$OMP PRIVATE(i,idust,ivar,ix,iy,il,ily,idim,ivn,ivt,eleft, eright, vleft, vright, vtleft, vtright, rholeft,rhoright, Pleft,Pright)

  !$OMP DO
  do ix=2,nx_max-1 
    do iy=2,ny_max-1
      il   = icell(ix+1,iy)
      i  = icell(ix,iy)
      do idim = 1,ndim
        ivn = iv
        if(idim==2)    i   = icell(ix,iy+1)
        if(idim==2)    ivn = ivy
        lambda_llf(i,idim) = max(abs(qp(i,ivn,idim))+sqrt(gamma*qp(i,iP,idim)/qp(i,irho,idim)),abs(qm(il,ivn,idim))+sqrt(gamma*qm(il,iP,idim)/qm(il,irho,idim)))
#if NDUST>0
        do idust=1,ndust
          ivn = ivd(idust)
          if(idim==2)    ivn = ivdy(idust)
          lambda_llf(i,idim) = max(lambda_llf(i,idim),max(abs(qp(i,ivn,idim)),abs(qm(il,ivn,idim))))
        end do
#endif      
      end do
    end do
  end do
  !$OMP END DO

  !$OMP BARRIER

  !$OMP DO
  do ix = 2,nx_max-1 
    do iy = 2,ny_max-1
      do idim = 1,ndim 
        il = icell(ix+1,iy)
        i = icell(ix,iy)
        ivn = iv
        ivt = ivy
        if(idim==2) then
          il = icell(ix,iy+1)
          ivn = ivy
          ivt = iv
        endif 
        rhoright  = qm(il,irho,idim)
        rholeft   = qp(i,irho,idim)
        vright    = qm(il,ivn,idim)
        vleft     = qp(i,ivn,idim)
        Pright    = qm(il,iP,idim)
        Pleft     = qp(i,iP,idim)
        eleft    = Pleft /(gamma-1.d0) + half*rholeft*vleft**2.
        eright   = Pright/(gamma-1.d0) + half*rhoright*vright**2.
        if(ndim==2) then
          vtright   = qm(il,ivt,idim)
          vtleft  = qp(i,ivt,idim)
          eleft    = eleft   + half*rholeft *vtleft**2.
          eright   = eright  + half*rhoright*vtright**2.
        endif
        flux(i,irho,idim)  = half  * (rholeft*vleft+rhoright*vright) &
        & - half*lambda_llf(il,idim)* (rhoright-rholeft)
        flux(i,ivn,idim)   = half  * (rholeft*vleft**2+Pleft+rhoright*vright**2+Pright)&
        & - half*lambda_llf(il,idim)* (rhoright*vright-rholeft*vleft)
        if(ndim==2) then
          flux(i,ivt,idim)   = half  * (rholeft*vleft*vtleft+rhoright*vright*vtright)&
          & - half*lambda_llf(il,idim)* (rhoright*vtright-rholeft*vtleft)
        endif
        flux(i,iP,idim)    = half*((eleft+Pleft)*vleft+(eright+Pright)*vright)&
        & - half*lambda_llf(il,idim)*(eright-eleft)

        !Dust fluxes
#if NDUST>0
        do idust=1,ndust
          ivn = ivd(idust)
          ivt = ivdy(idust)
          if(idim==2) then
            il = icell(ix,iy+1)
            ivn = ivdy(idust)
            ivt = ivd(idust)
          endif 
          rhoright  = qm(il,irhod(idust),idim)
          rholeft   = qp(i,irhod(idust),idim)
          vright    = qm(il,ivn,idim)
          vleft     = qp(i,ivn,idim)
          if(ndim==2) then
            vtright   = qm(il,ivt,idim)
            vtleft    = qp(i,ivt,idim)
          endif
          flux(i,irhod(idust),idim)  = half  * (rholeft*vleft+rhoright*vright) &
          & - half*lambda_llf(il,idim)* (rhoright-rholeft)
          flux(i,ivn,idim)   = half  * (rholeft*vleft**2+rhoright*vright**2)&
          & - half*lambda_llf(il,idim)* (rhoright*vright-rholeft*vleft)
          if(ndim==2) then
            flux(i,ivt,idim)   = half  * (rholeft*vleft*vtleft+rhoright*vright*vtright)&
            & - half*lambda_llf(il,idim)* (rhoright*vtright-rholeft*vtleft)
          endif
        end do
#endif
      end do
    end do
  end do
  !$OMP END DO

  !$OMP BARRIER

  !$OMP DO
  do ix = first_active,last_active
    do iy = first_active_y,last_active_y
        do ivar = 1, nvar
            i  = icell(ix,iy)
            il = icell(ix-1,iy)
            delta_U(i,ivar)=(flux(il,ivar,1)*surf(il,1)-flux(i,ivar,1)*surf(i,1))/vol(i)*dt
            if(ndim==2) then
              i  = icell(ix,iy)
              il = icell(ix,iy-1)
              delta_U(i,ivar)=delta_U(i,ivar)+(flux(il,ivar,2)*surf(il,2)-flux(i,ivar,2)*surf(i,2))/vol(i)*dt
            endif
        end do
      end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  !Update state vector 
  unew=unew+delta_U

  deallocate(lambda_llf)
  deallocate(delta_U)
  deallocate(flux)
end subroutine add_delta_u

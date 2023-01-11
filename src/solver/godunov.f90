! This routine performs the predictor operation
subroutine predictor
  use parameters
  use commons
  use units
  use OMP_LIB
  implicit none
  integer :: i,idust,ivar
  real(dp):: barotrop,cs_eos
  real(dp), dimension(1:ncells,1:nvar) :: dq,dq1
  real(dp), dimension(1:ncells,1:nvar) :: qpred

  real(dp), dimension(1:ncells) ::dqq,dPP,PP,dqq1,dPP1
  if(static) return
  call apply_boundaries(1,uold,ncells,nvar)

  ! Initialise to zero
  dq    = 0.0d0
  dpp   = 0.0d0
  dq1   = 0.0d0
  qr    = 0.0d0
  ql    = 0.0d0
  csr   = 0.0d0
  csl   = 0.0d0
  qpred = 0.0d0

  ! Computes primitive variables

  call ctoprim

  ! Computes slopes for the variables (and for P when using barotropic EOS)
  do ivar=1,nvar
     CALL calc_slope(q(:,ivar),dqq,dqq1)
     do i =1,ncells
        dq(i,ivar)=dqq(i)
        dq1(i,ivar)=dqq1(i)
     end do
  end do
  do i=1,ncells
     PP(i)=q(i,irho)*cs_eos(barotrop(q(i,irho)))**2.
  end do
  call calc_slope(PP,dPP,dPP1)

  !Predictor step
  call predicting_gas(qpred,dq,dPP) 
#if NDUST>0  
  call predicting_dust(qpred,dq)
#endif

  ! Computes left and right states
  qr=q
  ql=q
  do i =1,ncells
     do ivar=1,nvar
        qr(i,ivar)=qr(i,ivar)+dq1(i,ivar)*dx_r_cell(i)
        ql(i,ivar)=ql(i,ivar)-dq1(i,ivar)*dx_l_cell(i)
     end do
  end do
  qr=qr+qpred
  ql=ql+qpred

  do i=1,ncells
     csl(i) = cs_eos(barotrop(ql(i,irho)))
     csr(i) = cs_eos(barotrop(qr(i,irho)))
     cs(i)  = cs_eos(barotrop(q(i,irho)))
  end do

#if SHPERE==1 
    ! Sphere setups assumes zero gradient at the boundary
    ql(inner_bound,:) = qr(first_active,:)
    qr(inner_bound,:) = qr(first_active,:)
    ql(outer_bound,:) = ql(last_active,:)
    qr(outer_bound,:) = ql(last_active,:)
    csl(inner_bound)  = csr(first_active)
    csr(inner_bound)  = csr(first_active)
    csl(outer_bound)  = csl(last_active)
    csr(outer_bound)  = csl(last_active)
#endif

#if GRAVITY==1  
  call mtot
#endif  
end subroutine predictor

! This routine computes the slope (with a slope limiter) for the primitive variables
subroutine calc_slope(qq,dq,dqq)
  use parameters
  use commons
  use units
  use OMP_LIB
  implicit none
  integer :: i,idust,ivar,ifirst,ilast
  real(dp) :: T,barotrop,cs_eos
  real(dp), dimension(1:ncells),intent(inout)  :: dq,dqq
  real(dp), dimension(1:ncells),intent(in)     :: qq

  dq = 0.0d0
  dqq= 0.0d0


  ifirst = 2
  ilast  = ncells-1
#if SPHERE==1
    ifirst = first_active ! Impose zero gradient for collapse, so no need to compute below -> it must be 0
    ilast  = last_active  ! Impose zero gradient for collapse, so no need to compute above -> it must be 0
#endif
!MINMOD slope limiter
  do i = 2,ncells-1
     if(abs((qq(i)-qq(i-1))) < abs((qq(i+1)-qq(i)))) then
        dq(i) =(qq(i)-qq(i-1))
        dqq(i)=(qq(i)-qq(i-1))/dx_l(i)
     else if ((qq(i)-qq(i-1))*(qq(i+1)-qq(i))<0.0d0) then
        dq(i) = 0.0d0
        dqq(i)= 0.0d0
     else
        dq(i) = (qq(i+1)-qq(i))
        dqq(i)= (qq(i+1)-qq(i))/dx_r(i)
     endif
  end do

end subroutine calc_slope

! Predictor step for the gas

subroutine predicting_gas(qpred,dq,dpP)
  use parameters
  use commons
  use units
  use OMP_LIB
  implicit none
  integer :: i,ivar
  real(dp), dimension(1:ncells,1:nvar),intent(inout) :: qpred,dq
  real(dp), dimension(1:ncells)       ,intent(inout) :: dPP
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i)
  !$OMP DO
  do i=first_active,last_active
     qpred(i,irho)=qpred(i,irho)+dt/2.0d0*(-dq(i,irho)*q(i,iv)-dq(i,iv)*q(i,irho))/dx_c(i)
#if SPHERE==1
     !Spherical geometry source term
     qpred(i,irho)=qpred(i,irho)-dt*q(i,irho)*q(i,iv)/r_c(i)
#endif    
     qpred(i,iv)=qpred(i,iv)+dt/2.0d0*(-dq(i,iv)*q(i,iv)/dx_c(i)-dPP(i)/dx_c(i)/q(i,irho))
#if GRAVITY==1
     !Gravity source term
     qpred(i,iv)=qpred(i,iv)+dt/2.0d0*(-Mc(i)/(r_c(i)**2.+(l_soft/unit_l)**2.))
#endif
     !Force source term
     qpred(i,iv)=qpred(i,iv)+dt/2.0d0*force(i)

  end do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine predicting_gas

! This is the Riemmann solver : llf + Godunov scheme
subroutine add_delta_u
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust,ivar
  real(dp), dimension(1:ncells)        :: lambda_llf
  real(dp), dimension(1:ncells,1:nvar) :: flux,delta_U
#if GRAVITY>0  
  real(dp), dimension(1:ir_sink)       :: f_acc_mass,delta_mass
  real(dp) :: m_add,delta_md
#endif  
  if(static) return

  !Initialisation of the flux and delta U
  flux    = 0.0d0
  delta_U = 0.0d0

  ! Compute the max wave velocity
  do i =inner_bound,last_active
     lambda_llf(i)=max(abs(qr(i,iv))+abs(csr(i)),abs(ql(i+1,iv))+abs(csl(i+1)))
  end do

#if NDUST>0     
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do idust=1,ndust
     do i =inner_bound,last_active
        lambda_llf(i)=max(lambda_llf(i),max(abs(qr(i,ivd(idust))),abs(ql(i+1,ivd(idust)))))
     end do
  end do
  !$OMP END DO
 !$OMP END PARALLEL
#endif


  
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO  
  do i=inner_bound,last_active
     flux(i,irho) = (0.5d0*(ql(i+1,irho)*ql(i+1,iv)+qr(i,irho)*qr(i,iv))+0.5d0*(qr(i,irho)-ql(i+1,irho))*lambda_llf(i)) 
     flux(i,iv)   = (0.5d0*(ql(i+1,irho)*ql(i+1,iv)**2.+qr(i,irho)*qr(i,iv)**2.+ql(i+1,irho)*csl(i+1)**2.+qr(i,irho)*csr(i)**2.)+0.5d0*(qr(i,irho)*qr(i,iv)-ql(i+1,irho)*ql(i+1,iv))*lambda_llf(i))
  end do
  !$OMP END DO
  !$OMP END PARALLEL

#if NDUST>0  
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do idust=1,ndust
     do i=inner_bound,last_active
        flux(i,irhod(idust)) = (0.5d0*(ql(i+1,irhod(idust))*ql(i+1,ivd(idust))+qr(i,irhod(idust))*qr(i,ivd(idust)))+0.5d0*(qr(i,irhod(idust))-ql(i+1,irhod(idust)))*lambda_llf(i))
        flux(i,ivd(idust))   = (0.5d0*(ql(i+1,irhod(idust))*ql(i+1,ivd(idust))**2.+qr(i,irhod(idust))*qr(i,ivd(idust))**2.)+0.5d0*(qr(i,irhod(idust))*qr(i,ivd(idust))-ql(i+1,irhod(idust))*ql(i+1,ivd(idust)))*lambda_llf(i))
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
#endif

! Godunov scheme dU/dt + div F (U) = 0
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,ivar)
  !$OMP DO
  do ivar=1,nvar
     do i=first_active,last_active
        delta_U(i,ivar)=(flux(i,ivar)*surf_p(i)-flux(i-1,ivar)*surf_p(i-1))/vol_cell(i)*dt
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL


  !Update state vector 
  unew=unew-delta_U
  
end subroutine add_delta_u

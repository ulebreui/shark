! This routine performs the dust step
subroutine dust(verbose,outputing)
  use parameters
  use commons
  use units
  implicit none
  logical :: verbose,outputing

  !We re-compute the primitive variables to compute the dust source terms
  call apply_boundaries(1,uold,ncells,nvar) 
  call ctoprim
  call compute_B
  !We account for the gas-dust drag
  call dust_drag
  call apply_boundaries(1,uold,ncells,nvar) 
  !We accound for the dust evolution, ie growth/fragmentation   
  call ctoprim  
  if(growth)call dust_growth(verbose)      
  !call set_uold
  call apply_boundaries(1,uold,ncells,nvar) 

  !We compute the charge of the dust grains and the values of the resistivities
  if(charging.and.charging_all_the_time) call charge
  call apply_boundaries(1,uold,ncells,nvar)  ! Last application of boundaries
  
end subroutine dust


! This routine performs the dust predictor step
subroutine predicting_dust(qpred,dq)
  use parameters
  use commons
  use units
  use OMP_LIB
  implicit none
  integer :: i,idust,ind_rhod,ind_vd
  real(dp), dimension(1:ncells,1:nvar),intent(inout)  :: qpred,dq
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,ind_rhod,ind_vd)
  !$OMP DO
  do idust=1,ndust
     do i=first_active,last_active
        ind_rhod = irhod(idust)
        ind_vd   = ivd(idust)
        qpred(i,ind_rhod) = qpred(i,ind_rhod)+0.5d0*dt*(-dq(i,ind_rhod)*q(i,ind_vd)-dq(i,ind_vd)*q(i,ind_rhod))/dx_c(i)
#if SPHERE==1
        !Spherical geometry source term
        qpred(i,ind_rhod) = qpred(i,ind_rhod)-dt*q(i,ind_rhod)*q(i,ind_vd)/r_c(i)
#endif       
        qpred(i,ind_vd)   = qpred(i,ind_vd)+0.5d0*dt*(-dq(i,ind_vd)*q(i,ind_vd)/dx_c(i))
#if GRAVITY==1
        !Gravity source term
        qpred(i,ind_vd)   = qpred(i,ind_vd)+0.5d0*dt*(-Mc(i)/(r_c(i)**2.+(l_soft/unit_l)**2.))
#endif
        !Force source term
        qpred(i,ind_vd)   = qpred(i,ind_vd)+dt/2.0d0*force(i)
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine predicting_dust

! Dust drag is computed (implicitely)
subroutine dust_drag
  
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust
  real(dp),dimension(1:ndust):: alphak
  real(dp):: pn,rhon
  if(static) return
  !Re-calc distribution
  call distribution_dust(.false.)

  ! Here we apply the Krapp et al. implict scheme to compute the dust drag source terms

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(pn,rhon,alphak,i,idust)
  !$OMP DO
  do i=first_active,last_active
     pn   = uold(i,iv)
     rhon = uold(i,irho)
     do idust=1,ndust
        ! alpha = 1/ts
        alphak(idust) = 1.0d0/(sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(i,idust)/(q(i,irho)*cs(i)))
        if(coupled_dust)  alphak(idust)=1.0d0/(sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*(1d-25/unit_l)/(q(i,irho)*cs(i))) 

        pn   = pn+(alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*q(i,irhod(idust))*q(i,ivd(idust))
        rhon = rhon+(alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*q(i,irhod(idust))
     end do
     do idust = 1,ndust
        uold(i,ivd(idust)) = q(i,irhod(idust))*(q(i,ivd(idust))/(1.0d0+alphak(idust)*dt)+ (alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*pn/rhon)
     end do 
     uold(i,iv) = pn/rhon*uold(i,irho)
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Regularisation to avoid negative dust densities
  do idust=1,ndust
     uold(:,irhod(idust)) = max(uold(:,irho)*dust_ratio_min, uold(:,irhod(idust)))
  enddo
  
  
end subroutine dust_drag


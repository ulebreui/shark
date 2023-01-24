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
  !We account for the gas-dust drag
  call dust_drag
  call apply_boundaries(1,uold,ncells,nvar) 
  
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
  do i=1,ncells
   if(active_cell(i)==1) then
     ! Direction - x
     pn   = uold(i,iv)
     rhon = uold(i,irho)
     do idust=1,ndust
        ! alpha = 1/ts
        alphak(idust) = 1.0d0/(sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(i,idust)/(q(i,irho)*cs(i)))
        pn   = pn+(alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*q(i,irhod(idust))*q(i,ivd(idust))
        rhon = rhon+(alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*q(i,irhod(idust))
     end do
     do idust = 1,ndust
        uold(i,ivd(idust)) = q(i,irhod(idust))*(q(i,ivd(idust))/(1.0d0+alphak(idust)*dt)+ (alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*pn/rhon)
     end do 
     uold(i,iv) = pn/rhon*uold(i,irho)
      ! Direction - y
     pn   = uold(i,ivy)
     rhon = uold(i,irho)
     do idust=1,ndust
        ! alpha = 1/ts
        alphak(idust) = 1.0d0/(sqrt(pi*gamma/8.0d0)*(rhograin/unit_d)*sdust(i,idust)/(q(i,irho)*cs(i)))
        pn   = pn+(alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*q(i,irhod(idust))*q(i,ivdy(idust))
        rhon = rhon+(alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*q(i,irhod(idust))
     end do
     do idust = 1,ndust
        uold(i,ivdy(idust)) = q(i,irhod(idust))*(q(i,ivdy(idust))/(1.0d0+alphak(idust)*dt)+ (alphak(idust)*dt/(1.0d0+alphak(idust)*dt))*pn/rhon)
     end do 
     uold(i,ivy) = pn/rhon*uold(i,irho)    
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Regularisation to avoid negative dust densities
  do idust=1,ndust
     uold(:,irhod(idust)) = max(uold(:,irho)*dust_ratio_min, uold(:,irhod(idust)))
  enddo
  
  
end subroutine dust_drag


#if NDUST>0
! Routine to compute the charge according to Marchand et al., 2021. Probably the most horrible routine of the code so sit tight :).
subroutine charge
  use parameters
  use commons
  use units
  use OMP_LIB 
  implicit none
  integer :: i,idust,lowT
  real(dp), dimension(1:ndust) :: tau_k,alpha_k,n_k,l_grain_loc
  real(dp), dimension(1:ndust) :: t_sdust,sigmav_dust,sigmas_dust,omegas_dust
  real(dp) :: thetai,sigmav_ie,vi
  real(dp) :: t_sions,t_sel,vrms_i,vrms_el,sigmav_ions,sigmav_el,sigmas_ions,omegas_ions,sigmas_el,omegas_el
  real(dp) :: as_He_dust,as_He_ions,as_He_el,mu_i,mu_e
  real(dp) :: psi_loc,psi0,B_gauss,cs_eos
  real(dp) :: fpsi,dfdpsi,epsone
  real(dp) :: convergence_ionis,eps_psi,n_i_loc
  real(dp) :: dni_dpsi
  real(dp) :: eps_theta
  integer  :: niter_ionis,niter_ionis_max
  real(dp) :: T,barotrop,nH_loc
  real(dp) :: B_field,cross_sec

  as_He_dust = 1.28
  as_He_ions = 1.14
  as_He_el   = 1.16
  !Limit for epsilon cannot be 1.
  epsone     = 0.99999d0

  convergence_ionis = 1d4
  niter_ionis       = 0

  !Initial guess
  psi0      = -2.5d0
  thetai=stickeff_el*dsqrt(mu_ions*mH/m_el)
  niter_ionis_max=0
  do while(convergence_ionis>1d-6.and.niter_ionis<1000)
     psi_loc=psi0
     !f(psi)
     fpsi=(1.d0-psi_loc)/(epsone*thetai)*exp(-psi_loc)-1.0d0
     ! df/dpsi
     dfdpsi=-exp(-psi_loc)/(epsone*thetai)-(1.d0-psi_loc)/(epsone*thetai)*exp(-psi_loc)
     !Newton-Raphson step
     psi_loc=psi_loc-fpsi/dfdpsi
     !Convergence check
     convergence_ionis=abs(fpsi)
     niter_ionis=niter_ionis+1
     psi0=psi_loc
  end do

  !Now we compute the actual ionisation
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,T,l_grain_loc,sigmav_ie,tau_k,vi,alpha_k,n_k,nH_loc,convergence_ionis,psi_loc,eps_psi,eps_theta,n_i_loc,fpsi) &
  !$OMP PRIVATE(dni_dpsi,dfdpsi,niter_ionis,B_gauss,sigmav_dust,t_sdust,mu_i,mu_e, vrms_i, vrms_el,sigmav_el,sigmav_ions,t_sel,t_sions,sigmas_el)&
  !$OMP PRIVATE(sigmas_ions,sigmas_dust,omegas_el,omegas_ions,omegas_dust,lowT,cross_sec)
  !$OMP DO
  do i=1,ncells
     if(active_cell(i)==1) then
     ! Temperature
     T            = barotrop(u_prim(i,irho))
     ! Grain size in cm
     l_grain_loc=sdust(i,:)*unit_l
     ! Ion electron collisional cross-section
     sigmav_ie=(2d-7)/dsqrt(T/300.0d0)
     ! Reduced dust temperature
     tau_k = l_grain_loc*kB*T/e_el_stat**2.
     ! Ion thermal velocity
     vi=dsqrt(8.0d0*kB*T/(pi*mu_ions*mH))
     ! Variable alpha_k for convenience (used in Marchand et al., 21)
     alpha_k=dsqrt(8.0d0/(pi*tau_k))
     ! Dust number density
     do idust=1,ndust
        n_k(idust)=u_prim(i,irhod(idust))*unit_d/(mdust(i,idust)*unit_m)
     end do
     ! Gas number density
     nH_loc=u_prim(i,irho)*unit_nh

     !Psi0 determined, serves as a guess for psi unless it's not the first dt
     convergence_ionis=1d4
     niter_ionis=0
     psi_loc=psi0

     if(psi_old(i).ne.0.0d0) psi_loc=psi_old(i)
     lowT=0
     do while(convergence_ionis>epsilon_ionis .and.niter_ionis<nitermax_ionis)

        ! Electron fraction
        eps_psi  = (1.0d0-psi_loc)/thetai*exp(-psi_loc)

        ! Variable added for convenience. Should in principle be replaced everywhere but the expression are already horrible
        eps_theta= eps_psi*thetai
        ! Ion number density
        n_i_loc  = sum(-(psi_loc*tau_k(:)+(1.0d0-eps_theta**2.0)/(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0))*n_k(:)/(1.0d0-eps_psi))

        !f(psi). This is the function we try to zero.
        fpsi=sigmav_ie*eps_psi*n_i_loc**2.0/(x*nH_loc)+n_i_loc*vi/(x*nH_loc)*(sum(n_k(:)*pi*l_grain_loc(:)**2.0*((1.0d0-psi_loc)+(2.0d0/tau_k(:))*(eps_theta**2.0+eps_theta)/(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0))))-1.0d0


        ! dni_dpsi /!\ trust me it's correct :). Yes. I suffered
        dni_dpsi   = -1.0d0/(1.0d0-eps_psi)**2.*(-eps_psi*(2.0d0-psi_loc)/(1.0d0-psi_loc))&
             &*(sum(n_k(:)*(-eps_theta**4.0 + (4.0d0*thetai**2.-thetai**3.*alpha_k(:))*eps_psi**2.&
             &+(2.0d0*thetai*alpha_K(:)-4.0d0*thetai**2.0)*eps_psi +1.0d0 - thetai*alpha_k(:))/(1.0d0 + eps_theta*alpha_k(:) + eps_theta**2.0d0)**2.0))&
             &-(psi_loc/(1.0d0-eps_psi)*(-eps_psi*(2.0d0-psi_loc)/(1.0d0-psi_loc))+1.0d0)*sum(n_k(:)*tau_k(:))/(1.0d0-eps_psi)

        ! Actual derivative
        dfdpsi=sigmav_ie*n_i_loc/(x*nH_loc)*(n_i_loc*(-eps_psi*(2.0d0-psi_loc)/(1.0d0-psi_loc))+2.0d0*eps_psi*dni_dpsi)+2.0d0*n_i_loc*vi/&
             &(x*nH_loc)*(eps_theta**2.+eps_theta)*((dni_dpsi/n_i_loc+(-eps_psi*(2.0d0-psi_loc)/(1.0d0-psi_loc))*(2.0d0*eps_theta+1.0d0)/(eps_psi*(1.0+eps_theta)))&
             &*(sum(n_k(:)*pi*l_grain_loc(:)**2.0/(tau_k(:)*(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0d0))))&
             &+((-eps_psi*(2.0d0-psi_loc)/(1.0d0-psi_loc))*sum(n_k(:)*pi*l_grain_loc(:)**2.0*(2.0d0*eps_psi*thetai**2.+thetai*alpha_k(:))&
             &/(tau_k(:)*(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0d0)**2.0d0))))+n_i_loc*vi/(x*nH_loc)*(dni_dpsi/n_i_loc*(1.0d0-psi_loc)-1.0d0)*(sum(n_k(:)*pi*l_grain_loc(:)**2.0))

        !Psi^n+1=Psi^n-f(psi)/dfdpsi(psi)
        psi_loc=psi_loc-fpsi/dfdpsi
        if(psi_loc<psi0) then
           psi_loc = psi0*0.9999999
           lowT=1              
           if(convergence_ionis==1d5) then
              convergence_ionis=epsilon_ionis*0.1 !To not loop indefinitely
           else
              convergence_ionis=1d5  !We cheated so we shouldn't converge at this step
           endif 
        else if (psi_loc>0.0d0) then
           psi_loc=-1d-5
           convergence_ionis=1d4 !We cheated so we shouldn't converge at this step
        end if
        convergence_ionis=abs(fpsi)
        niter_ionis=niter_ionis+1
     end do

     niter_ionis_max=max(niter_ionis_max,niter_ionis)
     eps_psi=(1.0d0-psi_loc)/thetai*exp(-psi_loc)
     if(lowT>0)eps_psi=epsone
     ni(i)=0.0d0
     do idust=1,ndust
        zd(i,idust)=psi_loc*tau_k(idust)+(1.0d0-eps_psi**2.0*thetai**2.0)/(1.0d0+eps_psi*thetai*alpha_k(idust)+eps_psi**2*thetai**2.0)
        ni(i)= ni(i)-zd(i,idust)*n_k(idust)/(1.0d0-eps_psi)
        if(lowT>0) ni(i)=sqrt(eps_psi*x*nH_loc/sigmav_ie)
     end do
     ne(i)=eps_psi*ni(i)
     psi_old(i)=psi_loc
     if(lowT>0)then
        psi_old(i)=psi0
     end if

     ! Resistivity computation

     !B_gauss = min(B_0_lee*sqrt(u_prim(i,irho)*unit_nh/1d4),B_threshold)! Magnetic field
     B_gauss = B_0_lee*sqrt(u_prim(i,irho)*unit_nh/1d4)! Magnetic field

     sigmav_dust(:)=pi*(l_grain_loc(:))**2.*dsqrt(8.0*kB*T/(pi*2.0d0*mH))*(1.0d0+dsqrt(pi/(2.*tau_k(:))))

     t_sdust(:)=dsqrt(pi*gamma/8.0d0)*(rhograin)*(sdust(i,:)*unit_l)/(u_prim(i,irho)*unit_d*cs_eos(T)*unit_v)

     mu_i=2.0d0*mH*mu_ions*mH/(2.0d0*mH+mu_ions*mH)
     mu_e=2.0d0*mH*m_el/(m_el+2.0d0*mH)

     vrms_i =dsqrt(8.0d0*kB*T/(pi*mu_i))*1d-5  ! /!\ /!\ These velocities need to be in km/s for the Pinto & Galli 2008 fit
     vrms_el=dsqrt(8.0d0*kB*T/(pi*mu_e))*1d-5

     sigmav_el  = 3.16d-11*vrms_el**1.3
     sigmav_ions= 2.4d-9  *vrms_i**0.6

     t_sel  = 1.0d0/as_He_el*((m_el+2.0d0*mH)/(2.0d0*mH))/sigmav_el/nH_loc
     t_sions= 1.0d0/as_He_ions*((mu_ions*mH+2.0d0*mH)/(2.0d0*mH))/sigmav_ions/nH_loc

     sigmas_el     = (ne(i))*e_el_stat**2.*t_sel/m_el
     sigmas_ions   = (ni(i))*e_el_stat**2*t_sions/(mu_ions*mH)
     sigmas_dust(:)= n_k(:)*(zd(i,:)*e_el_stat)**2.*t_sdust(:)/(mdust(i,:)*unit_m)

     omegas_el     = -e_el_stat*B_gauss/clight/m_el
     omegas_ions   = e_el_stat *B_gauss/clight/(mu_ions*mH)
     omegas_dust   = zd(i,:)*e_el_stat*B_gauss/clight/(mdust(i,:)*unit_m)

     ! Conductivities
     sigma_o(i)=sum(sigmas_dust(:))+sigmas_ions+sigmas_el
     sigma_p(i)=sigmas_ions/(1.0d0+(omegas_ions*t_sions)**2.0)+sigmas_el/(1.0d0+(omegas_el*t_sel)**2.0)+sum(sigmas_dust(:)/(1.0d0+(omegas_dust(:)*t_sdust(:))**2.0))
     sigma_H(i)=-sigmas_ions*(omegas_ions*t_sions)/(1.0d0+(omegas_ions*t_sions)**2.0)-sigmas_el*(omegas_el*t_sel)/(1.0d0+(omegas_el*t_sel)**2.0)-sum(sigmas_dust(:)*(omegas_dust(:)*t_sdust(:))/(1.0d0+(omegas_dust(:)*t_sdust(:))**2.0))

     ! Resistivities
     eta_o(i)=1.0d0/sigma_o(i)
     eta_H(i)=sigma_H(i)/(sigma_p(i)*2.+sigma_h(i)**2.)
     eta_a(i)=sigma_p(i)/(sigma_p(i)**2.+sigma_H(i)**2.)-1.0d0/sigma_o(i)
     !Hall factors
     do idust=1,ndust
        gamma_d(i,idust)=t_sdust(idust)*omegas_dust(idust)
     end do
     endif
end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine charge
#endif
#if NDUST==0
subroutine charge
  use parameters
  use commons
  use units
  use OMP_LIB 
  implicit none
  integer :: i,idust,lowT
  real(dp) :: thetai,sigmav_ie,vi
  real(dp) :: t_sions,t_sel,vrms_i,vrms_el,sigmav_ions,sigmav_el,sigmas_ions,omegas_ions,sigmas_el,omegas_el
  real(dp) :: as_He_dust,as_He_ions,as_He_el,mu_i,mu_e
  real(dp) :: psi_loc,psi0,B_gauss,cs_eos
  real(dp) :: fpsi,dfdpsi,epsone
  real(dp) :: convergence_ionis,eps_psi,n_i_loc
  real(dp) :: dni_dpsi
  real(dp) :: eps_theta
  integer  :: niter_ionis,niter_ionis_max
  real(dp) :: T,barotrop,nH_loc
  real(dp) :: B_field,cross_sec

  as_He_ions = 1.14
  as_He_el   = 1.16
  !Limit for epsilon cannot be 1.
  epsone     = 0.99999d0

  convergence_ionis = 1d4
  niter_ionis=0

  !Now we compute the actual ionisation
  do i=1,ncells
  !if (active_cell(i)==1) then
     ! Temperature
     T            = barotrop(q(i,irho))
     ! Ion electron collisional cross-section
     sigmav_ie=(2d-7)/dsqrt(T/300.0d0)
     ! Reduced dust temperature
     ! Ion thermal velocity
     vi=dsqrt(8.0d0*kB*T/(pi*mu_ions*mH))

     ! Gas number density
     nH_loc=q(i,irho)*unit_nh
     ni(i)=sqrt(x*nH_loc/sigmav_ie)
     ne(i)=ni(i)

     ! Resistivity computation

     B_gauss = min(B_0_lee*sqrt(q(i,irho)*unit_nh/1d4),B_threshold)! Magnetic field

     mu_i=2.0d0*mH*mu_ions*mH/(2.0d0*mH+mu_ions*mH)
     mu_e=2.0d0*mH*m_el/(m_el+2.0d0*mH)

     vrms_i =dsqrt(8.0d0*kB*T/(pi*mu_i))*1d-5  ! /!\ /!\ These velocities need to be in km/s for the Pinto & Galli 2008 fit
     vrms_el=dsqrt(8.0d0*kB*T/(pi*mu_e))*1d-5

     sigmav_el  = 3.16d-11*vrms_el**1.3
     sigmav_ions= 2.4d-9*vrms_i**0.6

     t_sel  = 1.0d0/as_He_el*((m_el+2.0d0*mH)/(2.0d0*mH))/sigmav_el/nH_loc
     t_sions= 1.0d0/as_He_ions*((mu_ions*mH+2.0d0*mH)/(2.0d0*mH))/sigmav_ions/nH_loc

     sigmas_el     = (ne(i))*e_el_stat**2.*t_sel/m_el
     sigmas_ions   = (ni(i))*e_el_stat**2*t_sions/(mu_ions*mH)

     omegas_el     = -e_el_stat*B_gauss/clight/m_el
     omegas_ions   = e_el_stat*B_gauss/clight/(mu_ions*mH)

     ! Conductivities
     sigma_o(i)=sigmas_ions+sigmas_el
     sigma_p(i)=sigmas_ions/(1.0d0+(omegas_ions*t_sions)**2.0)+sigmas_el/(1.0d0+(omegas_el*t_sel)**2.0)
     sigma_H(i)=-sigmas_ions*(omegas_ions*t_sions)/(1.0d0+(omegas_ions*t_sions)**2.0)-sigmas_el*(omegas_el*t_sel)/(1.0d0+(omegas_el*t_sel)**2.0)
     ! Resistivities
     eta_o(i)=1.0d0/sigma_o(i)
     eta_H(i)=sigma_H(i)/(sigma_p(i)*2.+sigma_h(i)**2.)
     eta_a(i)=sigma_p(i)/(sigma_p(i)**2.+sigma_H(i)**2.)-1.0d0/sigma_o(i)

!endif
end do

end subroutine charge
#endif



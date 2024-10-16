#if NDUST>0
! Routine to compute the charge according to Marchand et al., 2021. Probably the most horrible routine of the code so sit tight :).
!Needs dust distribution with cgs properties
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
     !if(active_cell(i)==1) then
     ! Temperature --> taken from setup_commons.
         !or: T(i) if energy equation
      T            = barotrop(u_prim(i,irho))
#if TURB>0
      T = T_cloud
#endif

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

        ! print *, 'psi_loc',psi_loc
        ! print *, 'tau_k',tau_k
        ! print *, 'eps_psi',eps_psi
        ! print *, 'thetai',thetai
        ! print *, 'alpha_k',alpha_k

       

        if(lowT>0) ni(i)=sqrt(eps_psi*x*nH_loc/sigmav_ie)
     end do
     ne(i)=eps_psi*ni(i)
     psi_old(i)=psi_loc
     if(lowT>0)then
        psi_old(i)=psi0
     end if

     ! Resistivity computation

     if (res_Marchand) then

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
      end if


     !endif
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
     ! Temperature
     T            = barotrop(q(i,irho))
#if TURB==1
      T = T_cloud !from dustyTurb setup_commons
#endif
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

end do

end subroutine charge
#endif


subroutine resistivities_with_dust_inertia
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
  real(dp) :: T,barotrop,nH
  real(dp) :: B,cross_sec

  as_He_ions = 1.14
  as_He_el   = 1.16
  !Limit for epsilon cannot be 1.
  epsone     = 0.99999d0

  convergence_ionis = 1d4
  niter_ionis=0

  !Now we compute the actual ionisation
  do i=1,ncells
  !if (active_cell(i)==1) then

      !Magnetic field intensity
      !B = dsqrt(Bx(i)**2+By(i)**2+Bz(i)**2) !If computed directly by induction equation
      !Else: use analytical prescription
      B = B_0_lee*sqrt(u_prim(i,irho)*unit_nh/1d4)! Magnetic field for collapse
      T            = barotrop(q(i,irho))   !Temperature for collapse

#if TURB==1
#if MHD==1
      T = T_cloud
      B = dsqrt(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)
   
#endif
#endif
     ! Ion electron collisional cross-section
     sigmav_ie=(2d-7)/dsqrt(T/300.0d0)
     ! Reduced dust temperature
     ! Ion thermal velocity
     vi=dsqrt(8.0d0*kB*T/(pi*mu_ions*mH))

     ! Gas number density
     nH=q(i,irho)*unit_nh

     ! print*,'nh=',nH


     ! Resistivity computation


     mu_i=2.0d0*mH*mu_ions*mH/(2.0d0*mH+mu_ions*mH)
     mu_e=2.0d0*mH*m_el/(m_el+2.0d0*mH)

     vrms_i =dsqrt(8.0d0*kB*T/(pi*mu_i))*1d-5  ! /!\ /!\ These velocities need to be in km/s for the Pinto & Galli 2008 fit
     vrms_el=dsqrt(8.0d0*kB*T/(pi*mu_e))*1d-5

     sigmav_el  = 3.16d-11*vrms_el**1.3
     sigmav_ions= 2.4d-9*vrms_i**0.6

     t_sel  = 1.0d0/as_He_el*((m_el+2.0d0*mH)/(2.0d0*mH))/sigmav_el/nH
     t_sions= 1.0d0/as_He_ions*((mu_ions*mH+2.0d0*mH)/(2.0d0*mH))/sigmav_ions/nH

     sigmas_el     = (ne(i))*e_el_stat**2.*t_sel/m_el
     sigmas_ions   = (ni(i))*e_el_stat**2.*t_sions/(mu_ions*mH)

     ! print*,'ne(i)==',ne(i)
     ! print*,'ni(i)==',ni(i)


     ! print*,'sigmas_ions=',sigmas_ions
     ! print*,'sigmas_el=',sigmas_el



     if (electrons .eqv. .false.) then
           sigmas_el     = 0.0
           omegas_el = 0.0
      end if
     if (ions .eqv. .false.) then
           sigmas_ions     = 0.0
           omegas_ions = 0.0

      end if      

     omegas_el     = -e_el_stat*B/clight/m_el
     omegas_ions   = e_el_stat*B/clight/(mu_ions*mH)


     !Hall factors
     Hall_e(i) = omegas_el*t_sel
     Hall_i(i) = omegas_ions*t_sions


     ! print*,'omega_ions=',omegas_ions
     ! print*,'omega_el=',omegas_el

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

end subroutine resistivities_with_dust_inertia


#if MHD==1
subroutine effective_diffusion_coef_induction
  use parameters
  use commons
  use units
  use OMP_LIB 
  implicit none
  real(dp) :: B_norm,T,barotrop,nH,bx,by,bz,mu_i,mu_e,as_He_ions,as_He_el,vrms_i,vrms_el,sigmav_el,sigmav_ions,t_sel,t_sions,omegas_el,omegas_ions
  integer :: i



    as_He_ions = 1.14
    as_He_el   = 1.16

 do i=1,ncells
    !if (active_cell(i)==1) then


    B_norm=dsqrt(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)

    if (dusty_nonideal_MHD) then
        bx=q(i,iBx)/B_norm
        by=q(i,iBy)/B_norm
        bz=q(i,iBz)/B_norm

        ! print*,'bx=',bx
        ! print*,'by=',by
        ! print*,'bz=',bz


        eta_eff_yy(i) = (eta_o(i)+(bx**2+by**2)*eta_a(i)) !Resistivities must be in cm^2/s--> that is why we call res_units in solve.f90 first.
        eta_eff_yz(i)=(bx*eta_h(i)+by*bz*eta_a(i))

        eta_eff_zy(i) = (by*bz*eta_a(i)-eta_h(i)*bx) !Resistivities must be in cm^2/s
        eta_eff_zz(i)=(eta_o(i)+(bx**2+bz**2)*eta_a(i))

    endif

    if (dusty_nonideal_MHD_no_electron) then


    if(analytical_charging .eqv. .false.) then !!Compute Hall factors here
        
             nH=q(i,irho)*unit_nh
             B_norm = dsqrt(q(i,iBx)**2 + q(i,iBy)**2 + q(i,iBz)**2)

            T = barotrop(q(i,irho)/(mu_gas*mH))

#if TURB==1
            T = T_cloud

#endif




             mu_i=2.0d0*mH*mu_ions*mH/(2.0d0*mH+mu_ions*mH)
             mu_e=2.0d0*mH*m_el/(m_el+2.0d0*mH)

             vrms_i =dsqrt(8.0d0*kB*T/(pi*mu_i))*1d-5  ! /!\ /!\ These velocities need to be in km/s for the Pinto & Galli 2008 fit
             vrms_el=dsqrt(8.0d0*kB*T/(pi*mu_e))*1d-5

             sigmav_el  = 3.16d-11*vrms_el**1.3
             sigmav_ions= 2.4d-9*vrms_i**0.6

             t_sel  = 1.0d0/as_He_el*((m_el+2.0d0*mH)/(2.0d0*mH))/sigmav_el/nH
             t_sions= 1.0d0/as_He_ions*((mu_ions*mH+2.0d0*mH)/(2.0d0*mH))/sigmav_ions/nH

! 
             omegas_el     = -e_el_stat*B_norm/clight/m_el
             omegas_ions   = e_el_stat*B_norm/clight/(mu_ions*mH)


             Hall_e(i) = omegas_el*t_sel
             Hall_i(i) = omegas_ions*t_sions




    endif


        eta_eff_ohm(i) = -B_norm*clight/(Hall_i(i)*e_el_stat*(ni(i)+ne(i))*4*pi)
        eta_eff_Hall_y(i) = -q(i,iBx)*clight/(e_el_stat*(ni(i)+ne(i))*4*pi)
        eta_eff_Hall_z(i) = q(i,iBx)*clight/(e_el_stat*(ni(i)+ne(i))*4*pi) 


    end if
 end do

! print*,'eta_o=',eta_o
! print*,'eta_a=',eta_a
! print*,'eta_h=',eta_h

! print*,'sigma_o=',sigma_o
! print*,'sigma_a=',sigma_p
! print*,'sigma_h=',sigma_H

end subroutine effective_diffusion_coef_induction

#endif





#if MHD==1
subroutine res_units
  use parameters
  use commons
  use units
  use OMP_LIB 
  implicit none
  integer :: i

  do i=1,ncells

    eta_o(i) = eta_o(i)*clight**2/(4*pi) !from s to cm^2/s
    eta_a(i) = eta_a(i)*clight**2/(4*pi)
    eta_H(i) = eta_H(i)*clight**2/(4*pi)

end do




end subroutine res_units

#endif


#if MHD==1
#if NDUST>0
subroutine electric_field
  use parameters
  use commons
  use units
  use OMP_LIB 
  use slope_limiter

  implicit none

  real(dp) :: B_norm,bx,by,bz,E_i_x,E_ohm_x,E_H_x,E_ad_x,E_i_y,E_ohm_y,E_H_y,E_ad_y,E_i_z,E_ohm_z,E_H_z,E_ad_z,E_x_g,E_y_g,E_z_g
  real(dp) :: total_dust_current_x,total_dust_current_y,total_dust_current_z
  real(dp) :: v_e1x,v_e1y,v_e1z,v_e2x,v_e2y,v_e2z,v_e3x,v_e3y,v_e3z,v_i1x,v_i1y,v_i1z,v_i2x,v_i2y,v_i2z,v_i3x,v_i3y,v_i3z
  real(dp) :: dBy,dBz
  real(dp) :: dxBy,dxBz
  real(dp), dimension(1:ncells) ::  Bym,Byp,Bzm,Bzp 
  integer :: i,idust,ix,iy,il,ir,icell,iyy,ixx

    total_dust_current_x = 0.0d0
    total_dust_current_y = 0.0d0
    total_dust_current_z = 0.0d0


  do i=1,ncells
    if(active_cell(i)==1) then

        ix=ixx(i)
        iy=iyy(i)

        if(slope_type>0) then
            il = icell(ix-1,iy)
            ir = icell(ix+1,iy)
            dBy = slope_limit(2.0d0*(q(i,iBy) - q(il,iBy))/(dx(i,1)+dx(il,1)),2.0d0*(q(ir,iBy) - q(i,iBy))/(dx(ir,1)+dx(i,1)))
            dBz = slope_limit(2.0d0*(q(i,iBz) - q(il,iBz))/(dx(i,1)+dx(il,1)),2.0d0*(q(ir,iBz) - q(i,iBz))/(dx(ir,1)+dx(i,1)))

            !dBy = 2.0d0*(q(i,iBy) - q(il,iBy))/(dx(i,1)+dx(il,1))
            !dBz = 2.0d0*(q(i,iBz) - q(il,iBz))/(dx(i,1)+dx(il,1))

            !Infer magnetic field at cell surfaces


            Bym(i) = q(i,iBy) + half*dBy*dx(i,1)
            Byp(i) = q(i,iBy) - half*dBy*dx(i,1)
            Bzm(i) = q(i,iBz) + half*dBz*dx(i,1)
            Bzp(i) = q(i,iBz) - half*dBz*dx(i,1)

        endif
    endif

   end do


 if (dusty_nonideal_MHD) then
     do i=1,ncells
        !if (active_cell(i)==1) then

            do idust=1,ndust

                total_dust_current_x = total_dust_current_x + (1/clight)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)*(q(i,ivdx(idust))-q(i,ivx)) !Total (relative to neutral velocity) dust current
                total_dust_current_y = total_dust_current_y + (1/clight)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)*(q(i,ivdy(idust))-q(i,ivy))
                total_dust_current_z = total_dust_current_z + (1/clight)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)*(q(i,ivdz(idust))-q(i,ivz))!Charge!!

            end do



            B_norm=dsqrt(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)
            bx=q(i,iBx)/B_norm
            by=q(i,iBy)/B_norm
            bz=q(i,iBz)/B_norm


            ! dxBy=(qp(i,iBy,1)-qm(i,iBy,1))/dx(i,1)
            ! dxBz=(qp(i,iBz,1)-qm(i,iBz,1))/dx(i,1)
            dxBy=(Bym(i)-Byp(i))/dx(i,1)
            dxBz=(Bzm(i)-Bzp(i))/dx(i,1)






            E_i_x = (-q(i,ivy)*q(i,iBz)+q(i,ivz)*q(i,iBy))/clight
            E_i_y = (-q(i,ivz)*q(i,iBx)+q(i,ivx)*q(i,iBz))/clight
            E_i_z = (-q(i,ivx)*q(i,iBy)+q(i,ivy)*q(i,iBx))/clight

            E_ohm_x = eta_o(i)*(-total_dust_current_x)/clight
            E_ohm_y = eta_o(i)*(-dxBz-total_dust_current_y)/clight
            E_ohm_z = eta_o(i)*(dxBy-total_dust_current_z)/clight

            E_H_x = eta_H(i)*(-bz*dxBz-bz*total_dust_current_y)/clight
            E_H_y = eta_H(i)*(-bx*dxBy-bx*total_dust_current_z+bz*total_dust_current_x)/clight
            E_H_z = eta_H(i)*(bx*dxBz+bx*total_dust_current_y+by*total_dust_current_x)/clight

            E_ad_x = eta_a(i)*(-(by**2+bz**2)*total_dust_current_x+bx*by*(dxBz+total_dust_current_y)+bx*bz*(total_dust_current_z-dxBy))/clight
            E_ad_y = eta_a(i)*(-(bx**2+bz**2)*(dxBz+total_dust_current_y)-by*bz*dxBy+bx*by*total_dust_current_x+by*bz*total_dust_current_z)/clight
            E_ad_z = eta_a(i)*((bx**2+by**2)*dxBy+by*bz*dxBz-(bx**2+by**2)*total_dust_current_z+bx*bz*total_dust_current_x+by*bz*total_dust_current_y)/clight


            E_x(i) = E_i_x + E_ohm_x + E_H_x + E_ad_x
            E_y(i) = E_i_y + E_ohm_y + E_H_y + E_ad_y
            E_z(i) = E_i_z + E_ohm_z + E_H_z + E_ad_z

            E_x_g = E_ohm_x + E_H_x + E_ad_x !E in neutral (gas) reference frame = E + v_g/c x B
            E_y_g = E_ohm_y + E_H_y + E_ad_y
            E_z_g = E_ohm_z + E_H_z + E_ad_z

            !Ion and electron velocities

            v_e1x=Hall_e(i)/(1+Hall_e(i)**2)*(E_y_g*bz-E_z_g*by)
            v_e1y=Hall_e(i)/(1+Hall_e(i)**2)*(E_z_g*bx-E_x_g*bz)
            v_e1z=Hall_e(i)/(1+Hall_e(i)**2)*(E_x_g*by-E_y_g*bx)

            v_i1x=Hall_i(i)/(1+Hall_e(i)**2)*(E_y_g*bz-E_z_g*by)
            v_i1y=Hall_i(i)/(1+Hall_e(i)**2)*(E_z_g*bx-E_x_g*bz)
            v_i1z=Hall_i(i)/(1+Hall_e(i)**2)*(E_x_g*by-E_y_g*bx)


            v_e2x=1/(1+Hall_e(i)**2)*(by*(E_x_g*by-E_y_g*bx)-bz*(E_z_g*bx-E_x_g*bz))
            v_e2y=1/(1+Hall_e(i)**2)*(bz*(E_y_g*bz-E_z_g*by)-bx*(E_x_g*by-E_y_g*bx))
            v_e2z=1/(1+Hall_e(i)**2)*(bx*(E_z_g*bx-E_x_g*bz)-by*(E_y_g*bz-E_z_g*by))

            v_i2x=1/(1+Hall_i(i)**2)*(by*(E_x_g*by-E_y_g*bx)-bz*(E_z_g*bx-E_x_g*bz))
            v_i2y=1/(1+Hall_i(i)**2)*(bz*(E_y_g*bz-E_z_g*by)-bx*(E_x_g*by-E_y_g*bx))
            v_i2z=1/(1+Hall_i(i)**2)*(bx*(E_z_g*bx-E_x_g*bz)-by*(E_y_g*bz-E_z_g*by))


            v_e3x=(E_x_g*bx+E_y_g*by+E_z_g*bz)*bx
            v_e3y=(E_x_g*bx+E_y_g*by+E_z_g*bz)*by
            v_e3z=(E_x_g*bx+E_y_g*by+E_z_g*bz)*bz

            v_i3x=(E_x_g*bx+E_y_g*by+E_z_g*bz)*bx
            v_i3y=(E_x_g*bx+E_y_g*by+E_z_g*bz)*by
            v_i3z=(E_x_g*bx+E_y_g*by+E_z_g*bz)*bz

            v_e_x(i) = clight*Hall_e(i)/B_norm*(v_e1x + v_e2x + v_e3x) + q(i,ivx)
            v_e_y(i) = clight*Hall_e(i)/B_norm*(v_e1y + v_e2y + v_e3y) + q(i,ivy)
            v_e_z(i) = clight*Hall_e(i)/B_norm*(v_e1z + v_e2z + v_e3z) + q(i,ivz)

            v_i_x(i) = clight*Hall_i(i)/B_norm*(v_i1x + v_i2x + v_i3x) + q(i,ivx)
            v_i_y(i) = clight*Hall_i(i)/B_norm*(v_i1y + v_i2y + v_i3y) + q(i,ivy)
            v_i_z(i) = clight*Hall_i(i)/B_norm*(v_i1z + v_i2z + v_i3z) + q(i,ivz)



      end do  !endif
  endif


    if (dusty_nonideal_MHD_no_electron) then

    do i=1,ncells


        idust = i_coupled_species

        B_norm=dsqrt(q(i,iBx)**2+q(i,iBy)**2+q(i,iBz)**2)

        dxBy=(Bym(i)-Byp(i))/dx(i,1)
        dxBz=(Bzm(i)-Bzp(i))/dx(i,1)

        E_x(i) = 1/(4*pi*e_el_stat*(ni(i)+ne(i)))*half*(2*q(i,iBz)*dxBz+2*q(i,iBy)*dxBy) + (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/clight/(ni(i)+ne(i))*(q(i,ivdy(idust))*q(i,iBz)-q(i,ivdz(idust))*q(i,iBy)) + (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)*B_norm/(clight*Hall_i(i)*(ni(i)+ne(i)))*q(i,ivdx(idust)) + B_norm/(clight*Hall_i(i))*q(i,ivx)
        E_y(i) = 1/(4*pi*e_el_stat*(ni(i)+ne(i)))*(-q(i,iBx)*dxBy) + (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/clight/(ni(i)+ne(i))*(q(i,ivdz(idust))*q(i,iBx)-q(i,ivdx(idust))*q(i,iBz)) + B_norm/(4*pi*Hall_i(i)*e_el_stat*(ni(i)+ne(i)))*dxBz + (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)*B_norm/(clight*Hall_i(i)*(ni(i)+ne(i)))*q(i,ivdy(idust)) + B_norm/(clight*Hall_i(i))*q(i,ivy)
        E_z(i) = 1/(4*pi*e_el_stat*(ni(i)+ne(i)))*(-q(i,iBx)*dxBz) + (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/clight/(ni(i)+ne(i))*(q(i,ivdx(idust))*q(i,iBy)-q(i,ivdy(idust))*q(i,iBx)) - B_norm/(4*pi*Hall_i(i)*e_el_stat*(ni(i)+ne(i)))*dxBy + (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)*B_norm/(clight*Hall_i(i)*(ni(i)+ne(i)))*q(i,ivdz(idust)) + B_norm/(clight*Hall_i(i))*q(i,ivz)
        ! E_x(i) =  (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/clight/(ni(i)+ne(i))*(q(i,ivdy(idust))*q(i,iBz)-q(i,ivdz(idust))*q(i,iBy)) + (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)*B_norm/(clight*Hall_i(i)*(ni(i)+ne(i)))*q(i,ivdx(idust)) + B_norm/(clight*Hall_i(i))*q(i,ivx)
        ! E_y(i) =  (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/clight/(ni(i)+ne(i))*(q(i,ivdz(idust))*q(i,iBx)-q(i,ivdx(idust))*q(i,iBz)) + B_norm/(4*pi*Hall_i(i)*e_el_stat*(ni(i)+ne(i)))*dxBz + (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)*B_norm/(clight*Hall_i(i)*(ni(i)+ne(i)))*q(i,ivdy(idust)) + B_norm/(clight*Hall_i(i))*q(i,ivy)
        ! E_z(i) =  (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/clight/(ni(i)+ne(i))*(q(i,ivdx(idust))*q(i,iBy)-q(i,ivdy(idust))*q(i,iBx)) - B_norm/(4*pi*Hall_i(i)*e_el_stat*(ni(i)+ne(i)))*dxBy + (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)*B_norm/(clight*Hall_i(i)*(ni(i)+ne(i)))*q(i,ivdz(idust)) + B_norm/(clight*Hall_i(i))*q(i,ivz)


        v_i_x(i) =  -(q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/(ni(i)+ne(i))*q(i,ivdx(idust))
        v_i_y(i) = -clight/(4*pi*e_el_stat*(ni(i)+ne(i)))*dxBz -(q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/(ni(i)+ne(i))*q(i,ivdy(idust))
        v_i_z(i) = clight/(4*pi*e_el_stat*(ni(i)+ne(i)))*dxBy -(q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/(ni(i)+ne(i))*q(i,ivdz(idust))

    end do

    endif

end subroutine electric_field
#endif 
#endif



#if MHD==1
#if NDUST>0
subroutine Lorentz_force
  use parameters
  use commons
  use units
  use OMP_LIB 

  implicit none
  integer :: i,idust


!if (active_cell)

  if (dusty_nonideal_MHD) then
      do i=1,ncells
          do idust=1,ndust
             ! FLor_x_d(i,idust)=zd(i,idust)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*(E_x(i) + q(i,ivdy(idust))/clight*q(i,iBz) - q(i,ivdz(idust))/clight*q(i,iBy)) !Check expression
             ! FLor_y_d(i,idust)=zd(i,idust)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*(E_y(i) + q(i,ivdz(idust))/clight*q(i,iBx) - q(i,ivdx(idust))/clight*q(i,iBz))
             ! FLor_z_d(i,idust)=zd(i,idust)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*(E_z(i) + q(i,ivdx(idust))/clight*q(i,iBy) - q(i,ivdy(idust))/clight*q(i,iBx))

             ! FLor_x_d(i,idust)=zd(i,idust)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*(q(i,ivdy(idust))/clight*q(i,iBz) - q(i,ivdz(idust))/clight*q(i,iBy)) !Check expression
             ! FLor_y_d(i,idust)=zd(i,idust)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*(q(i,ivdz(idust))/clight*q(i,iBx) - q(i,ivdx(idust))/clight*q(i,iBz))
             ! FLor_z_d(i,idust)=zd(i,idust)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*(q(i,ivdx(idust))/clight*q(i,iBy) - q(i,ivdy(idust))/clight*q(i,iBx))

          end do

         ! FLor_x(i) = (e_el_stat*ni(i)*(E_x(i) + v_i_y(i)/clight*q(i,iBz) - v_i_z(i)/clight*q(i,iBy)) - e_el_stat*ne(i)*(E_x(i) + v_e_y(i)/clight*q(i,iBz) - v_e_z(i)/clight*q(i,iBy))) !ions + electrons. Check units
         ! FLor_y(i) = (e_el_stat*ni(i)*(E_y(i) + v_i_z(i)/clight*q(i,iBx) - v_i_x(i)/clight*q(i,iBz)) - e_el_stat*ne(i)*(E_y(i) + v_e_z(i)/clight*q(i,iBx) - v_e_x(i)/clight*q(i,iBz))) !ions + electrons. Check units
         ! FLor_z(i) = (e_el_stat*ni(i)*(E_z(i) + v_i_x(i)/clight*q(i,iBy) - v_i_y(i)/clight*q(i,iBx)) - e_el_stat*ne(i)*(E_z(i) + v_e_x(i)/clight*q(i,iBy) - v_e_y(i)/clight*q(i,iBx)))!ions + electrons. Check units
        ! FLor_x(i) = q(i,ivy)*q(i,iBz) - q(i,ivz)*q(i,iBy)
        ! FLor_y(i) = q(i,ivz)*q(i,iBx) - q(i,ivx)*q(i,iBz)
        ! FLor_z(i) = q(i,ivx)*q(i,iBy) - q(i,ivy)*q(i,iBx)



      end do
   endif


  if (dusty_nonideal_MHD_no_electron) then

        do i=1,ncells
          do idust=1,ndust
             FLor_x_d(i,idust)=zd(i,idust)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*(E_x(i) + q(i,ivdy(idust))/clight*q(i,iBz) - q(i,ivdz(idust))/clight*q(i,iBy)) !Check expression
             FLor_y_d(i,idust)=zd(i,idust)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*(E_y(i) + q(i,ivdz(idust))/clight*q(i,iBx) - q(i,ivdx(idust))/clight*q(i,iBz))
             FLor_z_d(i,idust)=zd(i,idust)*e_el_stat*(q(i,irhod(idust))/mdust(i,idust))*(E_z(i) + q(i,ivdx(idust))/clight*q(i,iBy) - q(i,ivdy(idust))/clight*q(i,iBx))

          end do

         FLor_x(i) = e_el_stat*ni(i)*(E_x(i) + v_i_y(i)/clight*q(i,iBz) - v_i_z(i)/clight*q(i,iBy))  !ions 
         FLor_y(i) = e_el_stat*ni(i)*(E_y(i) + v_i_z(i)/clight*q(i,iBx) - v_i_x(i)/clight*q(i,iBz))  !ions
         FLor_z(i) = e_el_stat*ni(i)*(E_z(i) + v_i_x(i)/clight*q(i,iBy) - v_i_y(i)/clight*q(i,iBx))  !ions 



      end do



  endif


end subroutine Lorentz_force
#endif 
#endif





#if NDUST>0

subroutine analytical_charge  !(Fujii et. al 2011) and see Lebreuilly 2020. 


    use parameters
    use commons
    use units
    use precision
    use functions_NR
    use newton_raphson

    implicit none

    integer :: i,idust,iters,maxiter
    logical :: debug
    real(dp) :: k_i_d,k_e_d,m_i,m_e,T,zeta,barotrop,y,fy,mu_i,mu_e,vrms_el,vrms_i,B_norm,nH,sigmav_el,sigmav_ions,t_sel,t_sions,omegas_ions,omegas_el,as_He_el,as_He_ions

    maxiter = 100000
    debug=.false.

    zeta = x !CR ionisation rate

    !print*, 'zeta', zeta
    m_i=mu_ions*mH
    m_e=m_el

    as_He_ions = 1.14
    as_He_el   = 1.16


    do idust=1,ndust
        do i=1,ncells

            T = barotrop(q(i,irho)/(mu_gas*mH))

#if TURB==1
            T = T_cloud

#endif


            !call function_dust_charge(zd(i,idust),zeta,q(i,irho),mdust(i,idust),q(i,irhod(idust)),sdust(i,idust),T)
            ! call function_prime_dust_charge(zd(i,idust),zeta,q(i,irho),mdust(i,idust),q(i,irhod(idust)),sdust(i,idust),T)

            call setValues(zeta,q(i,irho),mdust(i,idust),q(i,irhod(idust)),sdust(i,idust),T)
            call solve_newton_raphson(wrapper_dust_charge, wrapper_prime_dust_charge,-1.0d0,y,iters , maxiter, debug) !Update zd(i,idust)

            zd(i,idust) = y
            ! print *,'s_d',sdust(i,idust)
            ! print *,'m_d',mdust(i,idust)

            ! print*, 'zeta', zeta


            !print *,'iter=',iters
            !y = -1.0d0
            !fy = wrapper_dust_charge(y)
            !print*, 'f(y) = ',fy
            ! fy = wrapper_dust_charge(y=-1.0d0)
            ! print*, 'f(-1) = ',fy
            ! fy = wrapper_dust_charge(y=-1.0d0)/wrapper_prime_dust_charge(y=-1.0d0)
            ! print*, 'deltay(-1) = ',fy


            k_e_d = pi*sdust(i,idust)**2*dsqrt(8*kB*T/(pi*m_e))*exp((e_el_stat)**2*y/(sdust(i,idust)*kB*T))
            k_i_d = pi*sdust(i,idust)**2*dsqrt(8*kB*T/(pi*m_i))*(1 - (e_el_stat)**2*y/(sdust(i,idust)*kB*T))

            ni(i) = zeta*q(i,irho)*mdust(i,idust)/(k_i_d*mu_gas*mH*q(i,irhod(idust)))
            ne(i) = zeta*q(i,irho)*mdust(i,idust)/(k_e_d*mu_gas*mH*q(i,irhod(idust)))

            if(electrons .eqv. .false.) then

                ne(i) = 0.0d0

            endif


#if MHD==1

             !Hall factors

             nH=q(i,irho)*unit_nh
             B_norm = dsqrt(q(i,iBx)**2 + q(i,iBy)**2 + q(i,iBz)**2)



             mu_i=2.0d0*mH*mu_ions*mH/(2.0d0*mH+mu_ions*mH)
             mu_e=2.0d0*mH*m_el/(m_el+2.0d0*mH)

             vrms_i =dsqrt(8.0d0*kB*T/(pi*mu_i))*1d-5  ! /!\ /!\ These velocities need to be in km/s for the Pinto & Galli 2008 fit
             vrms_el=dsqrt(8.0d0*kB*T/(pi*mu_e))*1d-5

             sigmav_el  = 3.16d-11*vrms_el**1.3
             sigmav_ions= 2.4d-9*vrms_i**0.6

             t_sel  = 1.0d0/as_He_el*((m_el+2.0d0*mH)/(2.0d0*mH))/sigmav_el/nH
             t_sions= 1.0d0/as_He_ions*((mu_ions*mH+2.0d0*mH)/(2.0d0*mH))/sigmav_ions/nH


             omegas_el     = -e_el_stat*B_norm/clight/m_el
             omegas_ions   = e_el_stat*B_norm/clight/(mu_ions*mH)


             Hall_e(i) = omegas_el*t_sel
             Hall_i(i) = omegas_ions*t_sions

#endif


        end do
    end do



end subroutine analytical_charge

#endif

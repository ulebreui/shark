module smoluchowski
	contains

	subroutine dust_growth_shark(dt,K_0,CFL_growth,dust_dens,ndust,sdust,mdust,massmin,&
	   &eta,dvij,rho_gas,sticking_efficiency,eps_threshold,rhodust_min,&
	   &frag_test,Ebr_mono,m_mono,redistribute_fragments,eps_threshold_frag,kernel,frag_model,vfrag_code)

	   ! Smoluchowski coagulation/fragmentation in one cell (explicit, sub-cycled).
	   ! frag_model 0 (energy):   a fraction f_frag(Ecol/Ebr) of the colliding mass is fragmented, the rest coagulates.
	   ! frag_model 1 (velocity): relative velocities are Maxwellian around dv; the fraction of collisions with
	   !                          v > vfrag fragment completely, the others coagulate:
	   !                          p_frag = (1 + x) exp(-x), x = 1.5 (vfrag/dv)^2

	   use precision

	   implicit none
	   integer  :: ndust,ic1,ic2,idust,jdust,ifrag2
	   real(dp) :: dt,dt_growth,time_growth,massmin,eta,rho_gas,K_0,Ebr,Ecol

	   real(dp) :: CFL_growth
	   real(dp), dimension(1:ndust)           :: dust_dens,sdust,mdust
	   real(dp), dimension(1:ndust,1:ndust)   :: dvij
	   real(dp), dimension(1:ndust, 1:ndust)  :: redistribute_fragments

	   real(dp) :: f_frag,p_frag,p_coag,sticking_efficiency,eps_threshold_frag,eps_threshold

	   integer  :: niter_growth
	   integer  :: frag_test
	   integer  :: kernel          ! 0 physical, 1 constant (K = 1), 2 additive (K = m1 + m2)
	   integer  :: frag_model      ! 0 energy, 1 velocity threshold with Maxwellian distribution
	   real(dp) :: vfrag_code      ! fragmentation velocity (code units, used if frag_model = 1)
	   real(dp) :: Ebr_mono,m_mono,s1,s2,m1,m2
	   real(dp) :: epsilon_mass,rhodust_min,xfrag
	   real(dp) :: m_coag,frac_coag,frac_frag,m_frag_max

	   real(dp) , dimension(1:ndust) :: drhodt
	   real(dp) :: dndt
	   real(dp), parameter :: binmax_frag = 0.1d0   ! largest fragment / fragmented mass

	   niter_growth = 0
	   time_growth  = 0.0d0
	   ! Time loop
	   do while (time_growth < dt)
	      drhodt = 0.0d0
	      niter_growth = niter_growth + 1
	      do idust = 1, ndust
	         m1 = mdust(idust)
	         s1 = sdust(idust)
	         do jdust = 1, idust
	            m2 = mdust(jdust)
	            s2 = sdust(jdust)

	            ! Collision rate per unit volume: dndt = K n1 n2
	            select case (kernel)
	            case (1)
	               dndt = dust_dens(idust)*dust_dens(jdust)/m1/m2
	            case (2)
	               dndt = (m1 + m2)*dust_dens(idust)*dust_dens(jdust)/m1/m2
	            case default
	               dndt = K_0*(s1+s2)**2.*dvij(idust, jdust)*dust_dens(idust)*dust_dens(jdust)/m1/m2
	            end select
	            if (idust == jdust) dndt = dndt/2.0d0

	            ! Outcome of the collision
	            f_frag = 0.0d0                ! fragmented mass fraction per collision (energy model)
	            p_frag = 0.0d0                ! fraction of the collisions that fragment (velocity model)
	            p_coag = 1.0d0
	            if (frag_test == 1) then
	               if (frag_model == 1) then
	                  if (dvij(idust, jdust) > 0.0d0) then
	                     xfrag  = 1.5d0*(vfrag_code/dvij(idust, jdust))**2
	                     p_frag = (1.0d0 + xfrag)*exp(-xfrag)
	                     if (p_frag < 1d-10) p_frag = 0.0d0
	                  end if
	               else
	                  Ecol   = 0.5d0*(m1*m2)/(m1 + m2)*dvij(idust, jdust)**2
	                  Ebr    = (m1 + m2)/m_mono*Ebr_mono
	                  f_frag = max(min((Ecol - 0.1d0*Ebr)/(4.9d0*Ebr), 1.0d0), 0.0d0)
	               end if
	            end if

	            !Coagulation/frag conditions
	            if (eps_threshold > 0.0d0 .and. (dust_dens(idust) < eps_threshold*rho_gas)) p_coag = 0.0d0
	            if (eps_threshold > 0.0d0 .and. (dust_dens(jdust) < eps_threshold*rho_gas)) p_coag = 0.0d0
	            if (eps_threshold_frag > 0.0d0 .and. (dust_dens(idust) < eps_threshold_frag*rho_gas)) then
	               f_frag = 0.0d0; p_frag = 0.0d0
	            end if
	            if (eps_threshold_frag > 0.0d0 .and. (dust_dens(jdust) < eps_threshold_frag*rho_gas)) then
	               f_frag = 0.0d0; p_frag = 0.0d0
	            end if

	            ! Mass budget of the colliding pair
	            frac_coag  = (1.0d0 - p_frag)*(1.0d0 - f_frag)          ! fraction of the mass that coagulates
	            frac_frag  = p_frag + (1.0d0 - p_frag)*f_frag           ! fraction of the mass that is fragmented
	            m_coag     = (1.0d0 - f_frag)*(m1 + m2)                 ! mass of the coagulation product
	            m_frag_max = binmax_frag*(m1 + m2)                      ! largest fragment
	            if (frag_model /= 1) m_frag_max = m_frag_max*f_frag

	            ! Receiving bins of the coagulation product, bracketed by the representative masses:
	            ! mdust(ic1) <= m_coag < mdust(ic2), so that epsilon_mass is in [0,1]
	            ic1 = 1
	            ic2 = 1
	            epsilon_mass = 1.0d0
	            if (frac_coag > 0.0d0) then
	               ic1 = max(min(floor(dlog(m_coag/mdust(1))/dlog(eta)) + 1, ndust), 1)
	               do while (ic1 > 1 .and. mdust(ic1) > m_coag)
	                  ic1 = ic1 - 1
	               end do
	               do while (ic1 < ndust .and. mdust(min(ic1 + 1, ndust)) <= m_coag)
	                  ic1 = ic1 + 1
	               end do
	               ic2 = min(ic1 + 1, ndust)
	               if (ic1 < ic2) epsilon_mass = max(min((mdust(ic2) - m_coag)/(mdust(ic2) - mdust(ic1)), 1.0d0), 0.0d0)
	            end if

	            if (p_coag > 0.0d0) then
	               ! Losses
	               drhodt(idust) = drhodt(idust) - m1*dndt*(frac_coag*sticking_efficiency + frac_frag)
	               drhodt(jdust) = drhodt(jdust) - m2*dndt*(frac_coag*sticking_efficiency + frac_frag)
	               ! Coagulation
	               if (frac_coag > 0.0d0) then
	                  drhodt(ic1) = drhodt(ic1) + epsilon_mass*frac_coag*(m1 + m2)*dndt*sticking_efficiency
	                  drhodt(ic2) = drhodt(ic2) + (1.0d0 - epsilon_mass)*frac_coag*(m1 + m2)*dndt*sticking_efficiency
	               end if
	               ! Fragmentation: fragments up to m_frag_max (edges), at most the bin of the larger grain (idust >= jdust)
	               if (frac_frag > 0.0d0 .and. m_frag_max > 0.0d0) then
	                  ifrag2 = max(min(floor(dlog(m_frag_max/massmin)/dlog(eta)) + 1, idust), 1)
	                  drhodt(1:ifrag2) = drhodt(1:ifrag2) + frac_frag*(m1 + m2)*dndt*redistribute_fragments(1:ifrag2, ifrag2)
	               end if
	            end if
	         end do
	      end do

	      dt_growth = 1d32
	      do idust = 1, ndust
	         if (abs(drhodt(idust)) .ne. 0.0d0) then
	            dt_growth = min(dt_growth, dust_dens(idust)/abs(drhodt(idust))*CFL_growth)
	         end if
	      end do
	      dt_growth=min(dt_growth,dt-time_growth)
	      do idust = 1, ndust
	         dust_dens(idust) = max(dust_dens(idust) + drhodt(idust)*dt_growth, rhodust_min)
	      end do
	      time_growth = time_growth + dt_growth

	   end do

	end subroutine dust_growth_shark

	function dv_ormel(alpha_turb,cs,ts1,ts2,Reynolds,t_L)
	  use precision
	  implicit none
	  real(dp)  :: dv_ormel
	  real(dp)  :: alpha_turb,cs,ts1,ts2,Reynolds,t_L,x_stokes,f_Stokes,vclass2,vclass3,vclass1,St1,St2,t_eta  

       t_eta    = t_L/sqrt(Reynolds)                 ! Kolmogorov time
       x_stokes = min(ts2/ts1,ts1/ts2)
       f_Stokes = 3.2d0 - 1.0d0 - x_stokes + 2.0d0/(1.0d0+x_stokes)*(1.0d0/2.6d0 + x_stokes**3/(1.6d0 + x_stokes))
       St1 = max(ts1/t_l,ts2/t_l)
       St2 = min(ts1/t_l,ts2/t_l)

       vclass1 = sqrt(alpha_turb)*cs*sqrt((St1-St2)/(St1+St2))*dsqrt(St1**2/(St1+Reynolds**(-0.5))-St2**2/(St2+Reynolds**(-0.5)))
       vclass2 = sqrt(alpha_turb)*cs*sqrt(f_Stokes*St1)
       vclass3 = sqrt(alpha_turb)*cs*sqrt(1.0d0/(1.0d0 + St1) + 1.0d0/(1.0d0 + St2))

       dv_ormel = vclass2
       if (max(ts1,ts2) < t_eta) dv_ormel = vclass1   ! both grains in the small-eddy regime
       if (max(ts1,ts2) > t_L)   dv_ormel = vclass3   ! heavy grains

	end function dv_ormel

	function dv_brownian(w,m1,m2)
	  use precision
	  implicit none
	  real(dp)  :: dv_brownian
	  real(dp)  :: w,m1,m2 

      dv_brownian=  w*sqrt((m1+m2)/(m1*m2))

	end function dv_brownian


end module smoluchowski
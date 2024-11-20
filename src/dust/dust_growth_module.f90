module smoluchowski
	contains

	subroutine dust_growth_shark(dt,K_0,CFL_growth,dust_dens,ndust,sdust,mdust,massmin,&
	   &eta,dvij,rho_gas,sticking_efficiency,eps_threshold,rhodust_min,&
	   &frag_test,Ebr_mono,m_mono,redistribute_fragments,eps_threshold_frag)

	   use precision

	   implicit none
	   integer  :: ndust,ic1,ic2,idust,jdust,kdust,ifrag2
	   real(dp) :: dt,dt_growth,time_growth,massmin,eta,rho_gas,K_0,Ebr,Ecol

	   real(dp) :: CFL_growth
	   real(dp), dimension(1:ndust)           :: dust_dens,sdust,mdust
	   real(dp), dimension(1:ndust,1:ndust)   :: dvij
	   real(dp), dimension(1:ndust, 1:ndust)  :: redistribute_fragments

	   real(dp) :: f_frag,p_coag,sticking_efficiency,eps_threshold_frag,eps_threshold

	   integer  :: niter_growth
	   integer  :: frag_test
	   real(dp) :: Ebr_mono,m_mono,s1,s2,m1,m2
	   real(dp) :: epsilon_mass,rhodust_min

	   real(dp) , dimension(1:ndust) :: drhodt
	   real(dp) :: dndt

	   ! K_0 = pi

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
	            !#######################################
	            !#######################################
	            !#######################################
	            !#######################################
	            !Physical kernel
	            !#######################################
	            !#######################################
	            !#######################################
	            !#######################################

	            f_frag = 0.0d0
	            p_coag = 1.0d0

	            If (frag_test == 1) then

	               Ecol = 0.5d0*(m1*m2)/(m1 + m2)*dvij(idust, jdust)**2
	               Ebr = (m1 + m2)/m_mono*Ebr_mono
	               f_frag = max(min((Ecol - 0.1d0*Ebr)/(4.9d0*Ebr), 1.0d0), 0.0d0)

	            end if

	            dndt = K_0*(s1+s2)**2.*dvij(idust, jdust)*dust_dens(idust)*dust_dens(jdust)/m1/m2 ! K n1 n2
	            if (idust == jdust) dndt = dndt/2.0d0


	            !Coagulation/frag conditions
	            if (eps_threshold > 0.0d0 .and. (dust_dens(idust) < eps_threshold*rho_gas)) p_coag = 0.0d0
	            if (eps_threshold > 0.0d0 .and. (dust_dens(jdust) < eps_threshold*rho_gas)) p_coag = 0.0d0

	            if (eps_threshold_frag > 0.0d0 .and. (dust_dens(idust) < eps_threshold_frag*rho_gas)) f_frag = 0.0d0
	            if (eps_threshold_frag > 0.0d0 .and. (dust_dens(jdust) < eps_threshold_frag*rho_gas)) f_frag = 0.0d0

	            ic1 = max(min(floor(dlog((1.0d0 - f_frag)*(m1 + m2)/massmin)/dlog(eta) + 1), ndust), 1)
	            ic2 = max(min(floor(dlog((1.0d0 - f_frag)*(m1 + m2)/massmin)/dlog(eta) + 1) + 1, ndust), 1)

	            epsilon_mass = 1.0d0
	            if (ic1 < ic2) epsilon_mass = min((mdust(ic2) - (1.0d0 - f_frag)*(m1 + m2))/(mdust(ic2) - mdust(ic1)), 1.0d0)

	            ! Coagulation/sticking
	            if (p_coag .ge. 0.0d0) then
	               drhodt(idust) = drhodt(idust) - p_coag*(1.0d0 - f_frag)*m1*dndt*sticking_efficiency - p_coag*f_frag*m1*dndt
	               drhodt(jdust) = drhodt(jdust) - p_coag*(1.0d0 - f_frag)*m2*dndt*sticking_efficiency - p_coag*f_frag*m2*dndt
	               if (f_frag < 1.0d0) then
	                  drhodt(ic1) = drhodt(ic1) + epsilon_mass*p_coag*(1.0d0 - f_frag)*(m1 + m2)*dndt*sticking_efficiency
	                  drhodt(ic2) = drhodt(ic2) + (1.0d0-epsilon_mass)*p_coag*(1.0d0 - f_frag)*(m1 + m2)*dndt*sticking_efficiency
	               end if
	            end if

	            ! Fragmentation
	            if ((frag_test == 1 .and. f_frag > 0.0d0) .and. p_coag > 0.) then
	               ifrag2 = max(min(floor(dlog(0.1d0*f_frag*(m1 + m2)/massmin)/dlog(eta) + 1), ic1), 1)
	               do kdust = 1, ifrag2
	                  drhodt(kdust) = drhodt(kdust) + p_coag*f_frag*(m1 + m2)*dndt*redistribute_fragments(kdust, ifrag2)
	               end do
	            end if
	         end do
	      end do

	      dt_growth = 1d32!(dt-time_growth)
	      do idust = 1, ndust
	         if (abs(drhodt(idust)) .ne. 0.0d0) then
	            dt_growth = min(dt_growth, dust_dens(idust)/abs(drhodt(idust))*CFL_growth)
	         end if
	      end do
	      !stop
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

       x_stokes = min(ts2/ts1,ts1/ts2)
       f_Stokes = 3.2 - 1.0d0 - x_stokes + 2.0d0/(1.+x_stokes)*(1./2.6 + x_stokes**3./(1.6 + x_stokes))
       St1 = max(ts1/t_l,ts2/t_l)
       St2 = min(ts1/t_l,ts2/t_l)

       vclass1 = alpha_turb*cs*sqrt((St1-St2)/(St1+St2))*dsqrt(St1**2/(St1+Reynolds**(-0.5))-St2**2/(St2+Reynolds**(-0.5)))
       vclass2 = alpha_turb*cs*sqrt(f_Stokes*St1)
       vclass3 = alpha_turb*cs*sqrt(1.0d0/(1.0d0 + St1) + 1.0d0/(1.0d0 + St2))

       dv_ormel = vclass2
       if (ts1 < t_eta) dv_ormel  = vclass1
       if (ts1 > t_L)   dv_ormel  = vclass3

	end function dv_ormel

	function dv_brownian(w,m1,m2)
	  use precision
	  implicit none
	  real(dp)  :: dv_brownian
	  real(dp)  :: w,m1,m2 

      dv_brownian=  w*sqrt((m1+m2)/(m1*m2))

	end function dv_brownian


end module smoluchowski
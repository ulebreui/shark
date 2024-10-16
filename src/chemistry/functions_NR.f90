module functions_NR !!!!Function (to find root of with Newton-Raphson method) needed for analytical ionization prescription (Fujii et. al 2011) and see Lebreuilly 2020


	use commons
	use units
	use precision


	public :: function_dust_charge, wrapper_dust_charge, function_prime_dust_charge, wrapper_prime_dust_charge
	real(dp) :: zeta_fixed,rho_fixed,m_d_fixed,rho_d_fixed,s_d_fixed,T_fixed

	contains


	 real(dp) function function_dust_charge(y,zeta,rho,m_d,rho_d,s_d,T) 

	    use commons
	    use units
	    use precision

	    implicit none
	    real(dp), intent(in) :: y !!Variable == Zd
	    real(dp) :: zeta,rho,m_d,rho_d,s_d,T,m_i,m_e,k_e_d,k_i_d


	    m_i=mu_ions*mH
	    m_e=m_el

	    k_e_d = pi*s_d**2*dsqrt(8*kB*T/(pi*m_e))*exp((e_el_stat)**2*y/(s_d*kB*T))
	    k_i_d = pi*s_d**2*dsqrt(8*kB*T/(pi*m_i))*(1 - (e_el_stat)**2*y/(s_d*kB*T))


	    !f_dust_charge = mu_gas*mH*zeta/rho*(m_d*rho/rho_d/mu_gas/mH)**2*(1/k_i_d - 1/k_e_d) - y
    	function_dust_charge = mu_gas*mH*zeta/rho*(m_d*rho/(rho_d*mu_gas*mH))**2*((1/k_i_d) - (1/k_e_d)) - y

    	if (electrons .eqv. .false.) then

    		function_dust_charge = mu_gas*mH*zeta/rho*(m_d*rho/(rho_d*mu_gas*mH))**2*((1/k_i_d)) - y

    	endif

	    return 

	end function function_dust_charge

	real(dp) function wrapper_dust_charge(y) !Within which variables not to be passed in Newton_Raphen solver are fixed.
	    use commons
	    use units
	    use precision

	    implicit none
	    real(dp), intent(in) :: y 

	    wrapper_dust_charge = function_dust_charge(y,zeta_fixed,rho_fixed,m_d_fixed,rho_d_fixed,s_d_fixed,T_fixed)

	    return

	end function wrapper_dust_charge

	real(dp) function function_prime_dust_charge(y,zeta,rho,m_d,rho_d,s_d,T)


	    use commons
	    use units
	    use precision

	    implicit none
	    real(dp), intent(in) :: y !!Variable == Zd
	    real(dp) :: zeta,rho,m_d,rho_d,s_d,T,m_i,m_e,k_e_d,k_i_d,k_prime_e_d,k_prime_i_d

	    m_i=mu_ions*mH
	    m_e=m_el



	    k_e_d = pi*s_d**2*dsqrt(8*kB*T/(pi*m_e))*exp((e_el_stat)**2*y/(s_d*kB*T))
	    k_i_d = pi*s_d**2*dsqrt(8*kB*T/(pi*m_i))*(1 - (e_el_stat)**2*y/(s_d*kB*T))

	    k_prime_e_d = pi*s_d*dsqrt(8*kB*T/pi/m_e)*(e_el_stat)**2/kB/T*exp((e_el_stat)**2*(y)/s_d/kB/T)
	    k_prime_i_d = -pi*s_d*dsqrt(8*kB*T/pi/m_i)*(e_el_stat)**2/kB/T

	    !fprime_dust_charge = mu_gas*mH*zeta/rho*(m_d*rho/rho_d/mu_gas/mH)**2*(k_prime_e_d/k_e_d**2 - k_prime_i_d/k_i_d**2) - 1
	    function_prime_dust_charge = mu_gas*mH*zeta/rho*(m_d*rho/rho_d/mu_gas/mH)**2*(k_prime_e_d/k_e_d**2 - k_prime_i_d/k_i_d**2) - 1

    	if (electrons .eqv. .false.) then

	    	function_prime_dust_charge = mu_gas*mH*zeta/rho*(m_d*rho/rho_d/mu_gas/mH)**2*(- k_prime_i_d/k_i_d**2) - 1

    	endif

	    return 


	end function function_prime_dust_charge


	real(dp) function wrapper_prime_dust_charge(y)
	    use commons
	    use units
	    use precision

	    implicit none
	    real(dp), intent(in) :: y 


	    wrapper_prime_dust_charge = function_prime_dust_charge(y,zeta_fixed,rho_fixed,m_d_fixed,rho_d_fixed,s_d_fixed,T_fixed)

	    return

	end function wrapper_prime_dust_charge


	subroutine setValues(zeta,rho,m_d,rho_d,s_d,T)

	    use commons
	    use units
	    use precision
		implicit none
		real(dp) :: zeta,rho,m_d,rho_d,s_d,T

		zeta_fixed = zeta
		rho_fixed = rho
		m_d_fixed = m_d
		rho_d_fixed = rho_d
		s_d_fixed = s_d
		T_fixed = T

	end subroutine




end module functions_NR
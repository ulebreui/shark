double precision function barotrop(nH)
   use parameters
   use commons
   use units, ONLY : unit_dcol
   implicit none

   real(dp) :: nH,nstar,n1,n2,n3

   !barotrop = T0 *(1.0d0+(nH*unit_dcol/sigma_crit)**(gamma-1.))
   barotrop = T0 *(1.0d0+(nH*unit_dcol/sigma_crit)**(gamma-1.))
   !barotrop = T0*sqrt(1.0d0+(nH*unit_dcol/sigma_crit)**(0.8))!*(1.0d0+(nH*unit_nH/1d16))**(-0.3d0)*(1.0d0+(nH*unit_nH/1d21))**(0.56667d0)

end function barotrop

double precision function cs_eos(T)
  use parameters
  use commons
  use units
  implicit none
   real(dp) :: T
   cs_eos= sqrt(gamma*kB*T/(mu_gas*mH))/unit_v
end function cs_eos
double precision function barotrop(nH)
   use parameters
   use units, ONLY : unit_nH
   implicit none

   real(dp) :: nH,nstar
   nstar=1.0d-13/(mu_gas*mH)
   barotrop= T0_cloud*(1.0+(nH*unit_nH/nstar)**(gamma-1.0d0))!/((1.0+(nH*unit_nH/1e16))**(0.3))*((1.0+(nH*unit_nH/1e21))**0.56667)
end function barotrop

double precision function cs_eos(T)
  use parameters
  use units
  implicit none
   real(dp) :: T
   cs_eos= sqrt(gamma*kB*T/(mu_gas*mH))/unit_v

end function cs_eos

double precision function barotrop(nH)
   use parameters
   use units, ONLY : unit_nH
   implicit none

   real(dp) :: nH,nstar,n1,n2,n3
   if(ramses_baro) then
      nstar=1.0d-13/(mu_gas*mH)
      barotrop= T0_cloud*(1.0+(nH*unit_nH/nstar)**(gamma-1.0d0))
   else
      barotrop= T0_cloud!*dsqrt(1.0d0+(nH*unit_nH/1d11)**(0.8))*(1.0d0+(nH*unit_nH/1d16))**(-0.3d0)*(1.0d0+(nH*unit_nH/1d21))**(0.56667d0)
   endif
end function barotrop

double precision function cs_eos(T)
  use parameters
  use units
  implicit none
   real(dp) :: T
   cs_eos= dsqrt(gamma*kB*T/(mu_gas*mH))/unit_v

end function cs_eos

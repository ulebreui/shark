double precision function barotrop(nH)
   use parameters
   use units, ONLY : unit_nH
   implicit none

   real(dp) :: nH,nstar,n1,n2,n3

end function barotrop

double precision function cs_eos(T)
  use parameters
  use commons
  use units
  implicit none
   real(dp) :: T
   cs_eos= sqrt(gamma*kB*T/(mu_gas*mH))
end function cs_eos

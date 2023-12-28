module phys_const
  use precision
  implicit none

  real(dp),parameter:: e_el_stat     = 4.8d-10      ! Electron charge /!\ in statcoulomb here.
  real(dp),parameter:: m_el          = 9.109d-28    ! Electron mass in grams
  real(dp),parameter:: clight        = 3d10         ! Speed of light
  real(dp),parameter:: mH            = 1.67d-24     ! Mass of an hydrogen atom
  real(dp),parameter:: kB            = 1.38d-16     ! Boltzmann constant
  real(dp),parameter:: Grav          = 6.67d-8      ! Gravity constant
  real(dp),parameter:: au            = 1.5d13       ! Astronomical unit
  real(dp),parameter:: pc            = 3d18         ! Parsec
  real(dp),parameter:: Msun          = 2d33         ! Solar mass
  real(dp),parameter:: pi            = 3.141592653589793238462643383279d0  ! Pi

end module phys_const

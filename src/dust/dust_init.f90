!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This routine computes the dust distribution (either as IC or during the run depending on
! initi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine distribution_dust(initi)
   use parameters
   use commons
   use units
   use OMP_LIB

   implicit none
   logical  ::initi
   real(dp) :: zeta, eta, m_min, sum_dust, m_sum
   integer  :: idust, i, jdust, kdust, icell, ix, iy

   zeta = (smax/smin)**(1.0d0/ndust)
   eta = zeta**3.
   m_min = 4./3.*pi*rhograin*smin**3./unit_m
   if (kernel_type == 1) m_min = smin**3.
   if (kernel_type == 2) m_min = smin**3.

   if (initi) then
      do idust = 1, ndust
         aminus(idust) = smin*zeta**(idust - 1)/unit_l
         aplus(idust) = smin*zeta**(idust)/unit_l
         mminus(idust) = m_min*eta**(idust - 1)
         mplus(idust) = m_min*eta**(idust)
      end do
      do idust = 1, ndust
         do i = 1, ncells
            sdust(idust) = sqrt(aminus(idust)*aplus(idust))
         end do
      end do
   end if

   If (initi .and. (nrestart .eq. 0)) then !Else we use the updated dust distribution
      if (trim(dust_distribution) .eq. 'mrn' .and. kernel_type == 0) then
         call mrn_distri
      else if (trim(dust_distribution) .eq. 'paruta' .and. kernel_type == 0) then
         call paruta_distri
      else if (trim(dust_distribution) .eq. 'themis' .and. kernel_type == 0) then
         call themis_distri
      else if (trim(dust_distribution) .eq. 'logn' .and. kernel_type == 0) then
         call lognorm_distri
      else if (trim(dust_distribution) .eq. 'equal' .and. kernel_type == 0) then
         call equal_distri
      else if (trim(dust_distribution) .eq. 'table' .and. kernel_type == 0) then
         call table_distri
      else if (kernel_type == 1) then
         print *, 'constant kernel'
         mdust = sdust**3.
         do idust = 1, ndust
            do i = 1, ncells
               epsilondust(i, idust) = mdust(idust)*exp(-mdust(idust))*(mplus(idust) - mminus(idust))
            end do
         end do
      else if (kernel_type == 2) then
         print *, 'additive kernel'
         mdust = sdust**3.
         do idust = 1, ndust
            do i = 1, ncells
               epsilondust(i, idust) = mdust(idust)*exp(-mdust(idust))*(mplus(idust) - mminus(idust))
            end do
         end do
      else
         print *, 'unknown dust distribution'
         stop
      end if

   else !Else we use the updated dust distribution
      do iy = 1, ny_max
         do ix = 1, ny_max
            do idust = 1, ndust
               epsilondust(icell(ix, iy), idust) = q(ix, iy, irhod(idust))/q(ix, iy, irho)
            end do
         end do
      end do
      !end do
   end if

   !Add ice mantle
   if ((initi) .and. kernel_type == 0) then
      sdust = sdust + ice_mantle/unit_l
      aplus = aplus!+ice_mantle/unit_l
      aminus = aminus!+ice_mantle/unit_l
   end if
   if ((initi)) then
      print *, 'Initial dust distribution, size :'
      print *, sdust(:)*unit_l
      print *, 'Initial dust distribution, epsilon :'
      print *, epsilondust(1, :)
      print *, 'Initial dust distribution, mass :'
      print *, mdust(:)*unit_m
   end if

end subroutine distribution_dust

subroutine allocate_dust
   use parameters
   use commons
   use units
   implicit none

   allocate (aplus(1:ndust))
   allocate (aminus(1:ndust))
   allocate (mplus(1:ndust))
   allocate (mminus(1:ndust))
   allocate (epsilondust(1:ncells, 1:ndust))
   allocate (sdust(1:ndust))
   allocate (mdust(1:ndust))
   allocate (tstop(1:nx_max, 1:ny_max, 1:ndust))
   allocate (tcoag(1:nx_max, 1:ny_max, 1:ndust))

   allocate (force_dust_x(1:nx_max, 1:ny_max, 1:ndust))
   allocate (force_dust_y(1:nx_max, 1:ny_max, 1:ndust))
   allocate (force_dust_z(1:nx_max, 1:ny_max, 1:ndust))

   allocate (irhod(1:ndust))
   allocate (ivdx(1:ndust))
   allocate (ivdy(1:ndust))
   allocate (ivdz(1:ndust))

#if NDUSTPSCAL>0
   allocate (idust_pscal(1:ndust, 1:ndustpscal))
#endif
   !Size : bin edges
   aplus = 0.0d0
   aminus = 0.0d0
   !Mass : bin edges

   mplus = 0.0d0
   mminus = 0.0d0

   epsilondust = 0.0d0
   sdust = 0.0d0
   mdust = 0.0d0

   tstop = 0.0d0
   tcoag = 0.0d0
   force_dust_x = 0.0d0
   force_dust_y = 0.0d0
   force_dust_z = 0.0d0

end subroutine allocate_dust

subroutine mrn_distri
   use parameters
   use commons
   use units
   use OMP_LIB
   implicit none
   logical  :: initi
   real(dp) :: zeta, eta, m_min, sum_dust, m_sum, massmin, mono1, mono2, normalise
   integer  :: idust, i, jdust, kdust, icutmin, icutmax

   zeta = (smax/smin)**(1.0d0/ndust)
   eta = zeta**3.
   massmin = mminus(1)
   mono1 = 4./3.*pi*rhograin*scutmin**3./unit_m
   mono2 = 4./3.*pi*rhograin*scut**3./unit_m
   icutmin = max(1, floor(dlog(mono1/massmin)/log(eta) + 1))
   icutmax = min(floor(dlog(mono2/massmin)/log(eta) + 1), ndust)

   do i = 1, ncells
      normalise = 0.0d0
      do idust = icutmin, icutmax
         normalise = normalise + (aplus(idust)**(4.0 - mrn) - aminus(idust)**(4.0 - mrn))
         epsilondust(i, idust) = (aplus(idust)**(4.0 - mrn) - aminus(idust)**(4.0 - mrn))
      end do
      epsilondust(i, :) = dust2gas*epsilondust(i, :)/normalise

   end do

   epsilondust = max(epsilondust, dust_ratio_min)
   mdust = 4./3.*pi*rhograin/unit_d*sdust**3.

end subroutine mrn_distri

subroutine paruta_distri
   use parameters
   use commons
   use units
   use OMP_LIB
   implicit none
   logical  ::initi
   real(dp) :: zeta, eta, m_min, sum_dust, m_sum, massmin, mono1, mono2, normalise
   integer  :: idust, i, jdust, kdust, icutmin, icutmax

   zeta = (smax/smin)**(1.0d0/ndust)
   eta = zeta**3.
   massmin = mminus(1)
   mono1 = 4./3.*pi*rhograin*scutmin**3./unit_m
   mono2 = 4./3.*pi*rhograin*scut**3./unit_m
   icutmin = max(1, floor(dlog(mono1/massmin)/log(eta) + 1))
   icutmax = min(floor(dlog(mono2/massmin)/log(eta) + 1), ndust)
   mdust = 4./3.*pi*rhograin/unit_d*sdust**3.
   do i = 1, ncells
      normalise = 0.0d0
      do idust = icutmin, icutmax
         normalise = normalise + mdust(idust)**(-5./6.)*(mplus(idust) - mminus(idust))
         epsilondust(i, idust) = mdust(idust)**(-5./6.)*(mplus(idust) - mminus(idust))
      end do
      epsilondust(i, :) = dust2gas*epsilondust(i, :)/normalise
   end do

   epsilondust = max(epsilondust, dust_ratio_min)

end subroutine paruta_distri

subroutine lognorm_distri
   use parameters
   use commons
   use units
   use OMP_LIB
   implicit none
   logical  ::initi
   real(dp) :: zeta, eta, m_min, sum_dust, m_sum, renorm
   integer  :: idust, i, jdust, kdust
   real(dp), dimension(1:ndust):: distrib_no_norm, distrib_no_norm1

   do idust = 1, ndust
      distrib_no_norm(idust)=(derf(dlog(aplus(idust)*unit_l/aO_themis)/dsqrt(2.0d0*sigma_themis))-derf(dlog(aminus(idust)*unit_l/aO_themis)/dsqrt(2.0d0*sigma_themis)))
   end do
   distrib_no_norm = dust2gas*distrib_no_norm/sum(distrib_no_norm)

   do i = 1, ncells
      do idust = 1, ndust
         epsilondust(i, idust) = distrib_no_norm(idust)
      end do
   end do
   epsilondust = max(epsilondust, dust_ratio_min)
   mdust = 4./3.*pi*rhograin/unit_d*sdust**3.

end subroutine lognorm_distri

subroutine themis_distri
   use parameters
   use commons
   use units
   use OMP_LIB
   implicit none
   logical  ::initi
   real(dp) :: zeta, eta, m_min, sum_dust, m_sum, renorm
   integer  :: idust, i, jdust, kdust
   real(dp), dimension(1:ndust):: distrib_no_norm, distrib_no_norm1
   renorm = 0.0d0
   do idust = 1, ndust
     distrib_no_norm(idust)=half*(aplus(idust)-aminus(idust))*(aplus(idust)**(2)*exp(-(log(aplus(idust)*unit_l/aO_themis)/sigma_themis)**2.0d0)+aminus(idust)**(2)*exp(-(log(aminus(idust)*unit_l/aO_themis)/sigma_themis)**2.0d0))
      renorm = renorm + distrib_no_norm(idust)
   end do
   distrib_no_norm = dust2gas*distrib_no_norm/renorm*0.67d0
   renorm = 0.0d0
   do idust = 1, ndust
      distrib_no_norm1(idust)=half*(aplus(idust)-aminus(idust))*(aplus(idust)**(3.0-themis_slope)+aminus(idust)**(3.0-themis_slope))
     if(aplus(idust)>acut_themis/unit_l) distrib_no_norm1(idust)=  half*(aplus(idust)-aminus(idust))*(aplus(idust)**(3.0-themis_slope)*exp(-(aplus(idust)-acut_themis/unit_l/(awidthcut_themis/unit_l))**3.0)+aminus(idust)**(3.0-themis_slope)*exp(-(aminus(idust)-acut_themis/unit_l/(awidthcut_themis/unit_l))**3.0))
      renorm = renorm + distrib_no_norm1(idust)
   end do
   distrib_no_norm1 = dust2gas*distrib_no_norm1/renorm*0.33d0

   do i = 1, ncells
      do idust = 1, ndust
         epsilondust(i, idust) = distrib_no_norm(idust) + distrib_no_norm1(idust)
      end do
   end do
   epsilondust = max(epsilondust, dust_ratio_min)
   mdust = 4./3.*pi*rhograin/unit_d*sdust**3.

end subroutine themis_distri

subroutine equal_distri

   use parameters
   use commons
   use units
   use OMP_LIB
   implicit none
   logical  ::initi
   real(dp) :: zeta, eta, m_min, sum_dust, m_sum, renorm
   integer  :: idust, i, jdust, kdust

   do i = 1, ncells
      do idust = 1, ndust
         epsilondust(i, idust) = dust2gas/real(ndust)
      end do
   end do
   epsilondust = max(epsilondust, dust_ratio_min)
   mdust = 4./3.*pi*rhograin/unit_d*sdust**3.

end subroutine equal_distri

subroutine table_distri

   use parameters
   use commons
   use units
   use OMP_LIB
   implicit none
   integer :: idust, i
   character(len=200):: tabfile
   real(dp), dimension(1:ndust) :: epsi
   CALL getarg(2, tabfile)
   open (11, FILE=trim(tabfile))
   do idust = 1, ndust
      read (11, *) epsi(idust)
      print *, epsi(idust)
   end do
   do i = 1, ncells
      do idust = 1, ndust
         epsilondust(i, idust) = epsi(idust)
      end do
   end do
   close (11)
   epsilondust = max(epsilondust, dust_ratio_min)

   mdust = 4./3.*pi*rhograin/unit_d*sdust**3.

end subroutine table_distri

subroutine read_dust_params(ilun, nmlfile)

   use parameters
   use commons
   implicit none
   character(len=70):: nmlfile
   integer :: io, ilun
   logical::nml_ok
   namelist /dust_params/ frag_thre, vfrag, drag, dust_back_reaction, smin, smax, scut, scutmin, mrn, rhograin&
   &, dust2gas, growth, fragmentation, eps_threshold, eps_threshold_frag, growth_step &
   &, CFL_growth, rhodust_threshold, dust_ratio_min, dust_distribution, aO_themis, acut_themis, awidthcut_themis,&
   & themis_slope, sigma_themis, kernel_type, turb_in_growth, drift_in_growth, brownian_in_growth,&
   & slope_mono, ice_mantle, gamma_grains, estar_grains, sticking_efficiency, &
   & dtcontrol_growth, clustered_fraction, alpha_turb
  print *, "########################################################################################################################################"
  print *, "########################################################################################################################################"
   print *, "Dust namelist reading  !"
   read (13, dust_params, IOSTAT=io)
   rewind (13)
   if (io /= 0) then
      write (*, *) 'Invalid line in the dust namelist'
      stop
   end if
   print *, " Minimal grain size possible : ", smin
   print *, " Maximal grain size possible : ", smax
   print *, " Grain density               :", rhograin
   print *, " Dust-to-gas ratio           :", dust2gas
   print *, " Ice mantle                  :", ice_mantle

   if (trim(dust_distribution) .eq. 'mrn') then
      print *, ' MRN dust distribution'
      print *, " MRN index             : ", mrn
      print *, " Minimal grain size    : ", scutmin
      print *, " Maximal grain size    : ", scut
   else if (trim(dust_distribution) .eq. 'paruta') then
      print *, ' Paruta dust distribution'
      print *, " Minimal grain size    : ", scutmin
      print *, " Maximal grain size    : ", scut

   else if (trim(dust_distribution) .eq. 'themis') then
      print *, 'THEMIS dust distribution'
      print *, 'Average grain size       :', aO_themis
      print *, 'Standard dev             :', sigma_themis
      print *, 'Max monomer size         :', acut_themis
      print *, 'Width Max monomer size   :', awidthcut_themis
      print *, 'Slope of monomer distrib :', themis_slope

   else if (trim(dust_distribution) .eq. 'logn') then
      print *, ' THEMIS dust distribution'
      print *, ' Average grain size       :', aO_themis
      print *, ' Standard dev             :', sigma_themis
   else if (trim(dust_distribution) .eq. 'equal') then
      print *, ' Equal dust bins'
   else if (trim(dust_distribution) .eq. 'table') then
   else
      print *, 'unknown dust distribution'
      stop
   end if
   if (growth) then
      print *, "Growth is activated"
      print *, "Dust ratio threshold for growth :", eps_threshold
      print *, "Density threshold for growth    :", rhodust_threshold
      print *, "CFL for growth                  :", CFL_growth
      print *, 'gamma_grains                    :', gamma_grains
      print *, 'estar_grains                    :', estar_grains
      print *, 'sticking_efficiency             :', sticking_efficiency
      print *, 'clustered_fraction              :', clustered_fraction

      if (sticking_efficiency > 1 .or. sticking_efficiency <= 0) then
         print *, "sticking_efficiency must be larger than zero and smaller than 1"
         stop
      end if
      if (clustered_fraction < 1.) then
         print *, "clustered_fraction must at least be 1 "
         stop
      end if
      if (drift_in_growth) print *, "Growth by drift activated"
      if (turb_in_growth) print *, "Growth by turbulence (Ormel) activated"
      if (brownian_in_growth) print *, "Growth by brownian motion activated"
   else
      print *, "Growth is deactivated"
   end if
   if (fragmentation) then
      print *, "Fragmentation is activated"
      if (.not. growth) then
         print *, "You can't have fragmentation without growth"
         stop
      end if
      print *, "Power law index of monomers ", slope_mono
   else
      print *, "Fragmentation is deactivated"
   end if
   if (kernel_type == 0) print *, 'physical kernel'
   if (kernel_type == 1) then
      print *, 'Constant kernel'
   end if
   if (kernel_type == 2) then
      print *, 'Additive kernel'
   end if
   if (kernel_type == 3) then
      print *, 'wrong kernel'
      stop
   end if
  print *, "########################################################################################################################################"
  print *, "########################################################################################################################################"

end subroutine read_dust_params


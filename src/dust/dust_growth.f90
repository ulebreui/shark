
! This files contain the dust growth routines.
subroutine dust_growth(verbose)
   use parameters
   use commons
   use smoluchowski
   use units
   use OMP_LIB
   implicit none
   logical :: verbose

   integer :: ix, iy

   integer :: idust, jdust, kdust

   integer  :: frag_test, turbgrow, driftgrow, browgrow
   real(dp) :: eta, zeta, massmin, m1, m2, s1, s2, m_mono, a_mu
   real(dp) :: f_frag, p_frag, p_coag

   real(dp) :: T
   real(dp) :: Ecol, Ebr, Ebr_mono
   real(dp) :: t_L, t_eta, Reynolds
   real(dp) :: vdrift_turb, vdrift_brow, vdrift_hydro
   
   real(dp), dimension(1:ndust, 1:ndust)  :: dvij        
   real(dp), dimension(1:ndust, 1:ndust)  :: redistribute_fragments
   real(dp), dimension(1:ndust)           :: t_sdust
   real(dp), dimension(1:ndust)           :: dust_dens
   real(dp):: dt_growth, time_growth

   ! Flags
   frag_test = 0
   turbgrow  = 0
   driftgrow = 0
   browgrow  = 0

   if (fragmentation) frag_test     = 1
   if (turb_in_growth) turbgrow     = 1
   if (drift_in_growth) driftgrow   = 1
   if (brownian_in_growth) browgrow = 1

   ! Usefull dust grid quantities
   zeta    = (smax/smin)**(1.0d0/ndust)
   eta     = zeta**3.
   massmin = mminus(1)

   ! Monomer size and binding energy
   a_mu     = size_mono/2.0d0
   Ebr_mono = Abr*gamma_grains**(5.0d0/3.0d0)*(a_mu)**(4.0d0/3.0d0)/estar_grains**(2.0d0/3.0d0)/(unit_v**2*unit_m)
   m_mono   = 4.0d0/3.0d0*pi*rhograin*(size_mono)**3./unit_m

   ! Mass redistribution of fragments
   redistribute_fragments = 0.0d0
   do idust = 1, ndust
      do jdust = 1, ndust
         if(jdust<=idust) redistribute_fragments(jdust,idust) = (aplus(jdust)**(4.0d0-slope_mono)-aminus(jdust)**(4.0d0-slope_mono))/(aplus(idust)**(4.0d0-slope_mono)-aminus(1)**(4.0d0-slope_mono))
      end do
   end do
   do idust = 1, ndust
      redistribute_fragments(:, idust) = redistribute_fragments(:, idust)/sum(redistribute_fragments(:, idust))
   end do

   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         T        = (cs(ix,iy)*unit_v)**2*sqrt(mu_gas*mH/gamma/kB)
         Reynolds = 6.2d7*dsqrt(q(ix, iy, irho)*unit_d/(mu_gas*mH)/1d5)*dsqrt(T/10.0d0)
         t_L      = sqrt(3.*pi/(32.*grav*q(ix, iy, irho)*unit_d))/unit_t
         t_eta    = t_L/dsqrt(Reynolds)

         ! Differential velocity loop
         dvij = 0.0d0
         do idust = 1, ndust
            do jdust = 1, ndust
               if (turbgrow  == 1)  dvij(idust, jdust)  = dv_ormel(alpha_turb,cs(ix,iy),tstop(ix, iy, idust),tstop(ix, iy, jdust),Reynolds,t_L)
               if (browgrow  == 1)  dvij(idust, jdust)  = dsqrt(dvij(idust, jdust)**2.&
                  &+(dv_brownian(cs(ix,iy)*sqrt(mu_gas*mh/unit_m)/sqrt(pi*gamma/8.0d0),mdust(idust),mdust(jdust)))**2.)
               if (driftgrow == 1)  dvij(idust, jdust)  = dsqrt(dvij(idust, jdust)**2.&
                  &+(q(ix,iy,ivdx(idust))-q(ix,iy,ivdx(jdust)))**2&
                  &+(q(ix,iy,ivdy(idust))-q(ix,iy,ivdy(jdust)))**2&
                  &+(q(ix,iy,ivdz(idust))-q(ix,iy,ivdz(jdust)))**2)
            end do
         end do
      end do
   end do

   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         do idust = 1, ndust
              dust_dens(idust) = u_prim(ix, iy, irhod(idust))
         end do

         call dust_growth_shark(dt,pi,CFL_growth,dust_dens,ndust,sdust,mdust,massmin,&
         &eta,dvij,u_prim(ix, iy, irho),sticking_efficiency,eps_threshold,u_prim(ix, iy, irho)*dust_ratio_min,&
         &frag_test,Ebr_mono,m_mono,redistribute_fragments,eps_threshold_frag)

         do idust = 1, ndust
            u_prim(ix, iy, irhod(idust)) = dust_dens(idust)
         end do         
      end do
   end do

end subroutine dust_growth


! This files contain the dust growth routines.
subroutine dust_growth_stepinski
   use parameters
   use commons
   use units
   use OMP_LIB
   implicit none
   logical :: verbose
   integer :: idust, jdust, ipscal, ix, iy

   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         ! Differential velocity loop
         do idust = 1, ndust
           u_prim(ix,iy,idust_pscal(idust,1))=max(u_prim(ix,iy,idust_pscal(idust,1))*(1.0d0+dt/tcoag(ix,iy,idust)),q(ix,iy,irhod(idust))*sminstep/unit_l)
         end do
      end do
   end do

end subroutine dust_growth_stepinski

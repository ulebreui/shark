
! This files contain the dust growth routines.
subroutine dust_growth(verbose)
   use parameters
   use commons
   use smoluchowski
   use units
   use OMP_LIB
   implicit none
   logical :: verbose

   integer :: ix,iy,icell

   integer :: idust, jdust, kdust

   integer  :: frag_test, bouncing_test, turbgrow, driftgrow, browgrow
   real(dp) :: eta, zeta, massmin, m1, m2, s1, s2, m_mono, a_mu
   real(dp) :: f_frag, p_frag, p_coag

   real(dp) :: T
   real(dp) :: Ecol, Ebr, Ebr_mono
   real(dp) :: t_L, t_eta, Reynolds
   real(dp) :: vdrift_turb, vdrift_brow, vdrift_hydro, cs_modified
   
   real(dp), dimension(1:ndust, 1:ndust)  :: dvij        
   real(dp), dimension(1:ndust, 1:ndust)  :: redistribute_fragments
   real(dp), dimension(1:ndust)           :: t_sdust
   real(dp), dimension(1:ndust)           :: dust_dens
   real(dp), dimension(1:ndust)           :: v_dust_x_tampon
   real(dp), dimension(1:ndust)           :: v_dust_y_tampon
   real(dp), dimension(1:ndust)           :: v_dust_z_tampon

   real(dp):: dt_growth, time_growth

   ! Flags
   frag_test = 0
   bouncing_test = 0
   turbgrow  = 0
   driftgrow = 0
   browgrow  = 0

   if (fragmentation) frag_test     = 1
   if (bouncing) bouncing_test     = 1
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

!$omp parallel do default(shared) schedule(RUNTIME) private(ix,iy,idust, jdust, T, t_L, t_eta, Reynolds,dvij,dust_dens)
   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         T        = (cs(ix,iy)*unit_v)**2*sqrt(mu_gas*mH/gamma/kB)
         Reynolds = 6.2d7*dsqrt(q(irho,ix,iy)*unit_d/(mu_gas*mH)/1d5)*dsqrt(T/10.0d0)
         t_L      = sqrt(3.*pi/(32.*grav*q(irho,ix,iy)*unit_d))/unit_t
         t_eta    = t_L/dsqrt(Reynolds)


         if (SI .eqv. .true. .and. turbgrow  == 1) then


            t_L = 1/Omega_shear !Use the right dynamical time for SI to compute St correctly for Ormel. Alternatively, provide directly the St array defined in SI setup
            if (modified_Ormel .eqv. .true.) then !Modify soundspeed due to dust backreaction

               cs_modified = cs(ix,iy)/dsqrt(1+SUM(epsilondust(icell(ix,iy),:)/(1+St(:,ix,iy))))

            endif

         endif

         ! Differential velocity loop
         dvij = 0.0d0
         do idust = 1, ndust
            do jdust = 1, ndust
               if (turbgrow  == 1)  dvij(idust,jdust)  = dv_ormel(alpha_turb,cs(ix,iy),epsilondust(icell(ix,iy),idust),epsilondust(icell(ix,iy),jdust),tstop(idust,ix,iy),tstop(jdust,ix,iy),Reynolds,t_L,SI,modified_Ormel,cs_modified)

               if (browgrow  == 1)  dvij(idust,jdust)  = dsqrt(dvij(idust,jdust)**2.&
                  &+(dv_brownian(cs(ix,iy)*sqrt(mu_gas*mh/unit_m)/sqrt(pi*gamma/8.0d0),mdust(idust),mdust(jdust)))**2.)
               if (driftgrow == 1)  dvij(idust, jdust)  = dsqrt(dvij(idust, jdust)**2.&
                  &+(q(ivdx(idust),ix,iy)-q(ivdx(jdust),ix,iy))**2&
                  &+(q(ivdy(idust),ix,iy)-q(ivdy(jdust),ix,iy))**2&
                  &+(q(ivdz(idust),ix,iy)-q(ivdz(jdust),ix,iy))**2)
            end do
         end do
         do idust = 1, ndust
              dust_dens(idust) = u_prim(irhod(idust),ix,iy)
              ! if (SI) then 
              !  v_dust_x_tampon(idust) = u_prim(ivdx(idust),ix,iy)/u_prim(irhod(idust),ix,iy)
              !  v_dust_y_tampon(idust) = u_prim(ivdy(idust),ix,iy)/u_prim(irhod(idust),ix,iy)
              !  v_dust_z_tampon(idust) = u_prim(ivdz(idust),ix,iy)/u_prim(irhod(idust),ix,iy)
              ! endif

         end do

         call dust_growth_shark(dt,pi,CFL_growth,dust_dens,ndust,sdust,mdust,massmin,&
         &eta,dvij,u_prim(irho,ix,iy),sticking_efficiency,eps_threshold,u_prim(irho,ix,iy)*dust_ratio_min,&
         &frag_test,bouncing_test,Ebr_mono,m_mono,redistribute_fragments,eps_threshold_frag,SI,vfrag,v_bouncing)

         do idust = 1, ndust
            u_prim(irhod(idust),ix,iy) = dust_dens(idust)
            ! if (SI) then
            !    u_prim(ivdx(idust),ix,iy)=u_prim(irhod(idust),ix,iy)*v_dust_x_tampon(idust)
            !    u_prim(ivdy(idust),ix,iy)=u_prim(irhod(idust),ix,iy)*v_dust_y_tampon(idust)
            !    u_prim(ivdz(idust),ix,iy)=u_prim(irhod(idust),ix,iy)*v_dust_z_tampon(idust)

            ! endif
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
   integer :: idust, jdust, ipscal, ix,iy

   do iy = first_active_y, last_active_y
      do ix = first_active, last_active
         ! Differential velocity loop
         do idust = 1, ndust
           u_prim(idust_pscal(idust,1),ix,iy)=max(u_prim(idust_pscal(idust,1),ix,iy)*(1.0d0+dt/tcoag(ix,iy,idust)),q(irhod(idust),ix,iy)*sminstep/unit_l)
         end do
      end do
   end do

end subroutine dust_growth_stepinski

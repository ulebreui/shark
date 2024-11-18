
! This files contain the dust growth routines.
subroutine dust_growth(verbose)
  use parameters
  use commons
  use units
  use OMP_LIB
  implicit none
  logical :: verbose
  integer :: i,idust,jdust,kdust,niter_growth,ifrag1,ifrag2,ic1,ic2,ic3,ic4,ix,iy,icell
  integer :: frag_test,turbgrow,ambigrow,driftgrow,browgrow

  real(dp) :: eta,zeta,massmin,mono1,mono2,m_add,m1,m2,s1,s2,m_mono,a_mu,barotrop
  real(dp) :: mcoag1,mcoag2,f_frag
  real(dp) :: Kernel,dv,T,Ecol,Eroll,Ebr
  real(dp) :: cs_eos
  real(dp) :: Ebr_mono,Erol_mono
  real(dp) :: t_L,t_eta,Reynolds,St1,St2,vclass1,vclass2,vclass3,f_Stokes,x_stokes,vdrift_turb,vdrift_brow,vdrift_hydro,vdrift_ad
  real(dp) :: dndt,dndt2,p_frag,p_coag,epsilon_mass
  real(dp), dimension(1:ndust)          :: drhodt        ! Coagulation rate
  real(dp), dimension(1:ndust,1:ndust)  :: K_coag , dvij             ! Kernel and differential velocity
  real(dp), dimension(1:ndust,1:ndust)  :: redistribute_fragments
  real(dp), dimension(1:ndust)          :: t_sdust,sigmav_dust,tau_k,omega_dust,gamma_d2
  real(dp):: dt_growth,time_growth

  ! Flags 
  frag_test = 0
  turbgrow  = 0
  driftgrow = 0
  browgrow  = 0
  if(fragmentation)       frag_test = 1
  if(turb_in_growth)      turbgrow  = 1
  if(drift_in_growth)     driftgrow = 1
  if(brownian_in_growth)  browgrow  = 1

  ! Usefull dust grid quantities
  zeta   =(smax/smin)**(1.0d0/ndust)
  eta    = zeta**3.
  massmin= mminus(1)
  
  ! Monomer size and binding energy
  a_mu     = size_mono/2.0d0
  Ebr_mono = Abr*gamma_grains**(5.0d0/3.0d0)*(a_mu)**(4.0d0/3.0d0)/estar_grains**(2.0d0/3.0d0)/(unit_v**2*unit_m)
  m_mono   = 4.0d0/3.0d0*pi*rhograin*(size_mono)**3./unit_m

  ! Mass redistribution of fragments
  redistribute_fragments=0.0d0
  do idust=1,ndust
     do jdust=1,ndust
        if(jdust<=idust)redistribute_fragments(jdust,idust)= (aplus(jdust)**(4.0d0-slope_mono)-aminus(jdust)**(4.0d0-slope_mono))/(aplus(idust)**(4.0d0-slope_mono)-aminus(1)**(4.0d0-slope_mono))
     end do
  end do
  do idust=1,ndust
    redistribute_fragments(:,idust)= redistribute_fragments(:,idust)/sum(redistribute_fragments(:,idust))
  end do

   do iy = first_active_y,last_active_y
   do ix = first_active,last_active
      i =icell(ix,iy)
     ! Cell by cell-cycling
     time_growth  = 0.0d0
     niter_growth = 0
     T            = barotrop(q(ix,iy,irho))
     Reynolds     = 6.2d7*dsqrt(q(ix,iy,irho)*unit_d/(mu_gas*mH)/1d5)*dsqrt(T/10.0d0)
     t_L          = sqrt(3.*pi/(32.*grav*q(ix,iy,irho)*unit_d))/unit_t
     t_eta        = t_L/dsqrt(Reynolds)

     ! Differential velocity loop
     dvij=0.0d0           
     do idust=1,ndust
        do jdust=1,idust

           x_stokes = tstop(ix,iy,jdust)/tstop(ix,iy,idust)
           f_Stokes = 3.2-1.0d0-x_stokes+2.0d0/(1.+x_stokes)*(1./2.6+x_stokes**3./(1.6+x_stokes))
           St1 = tstop(ix,iy,idust)/t_l
           St2 = tstop(ix,iy,jdust)/t_l
           
           vclass1 = alpha_turb*cs_eos(T)*dsqrt((St1-St2)/(St1+St2))*dsqrt(St1**2/(St1+Reynolds**(-0.5))-St2**2/(St2+Reynolds**(-0.5)))
           vclass2 = alpha_turb*cs_eos(T)*dsqrt(f_Stokes*St1)           
           vclass3 = alpha_turb*cs_eos(T)*dsqrt(1.0d0/(1.0d0+St1)+1.0d0/(1.0d0+St2))
           
           vdrift_turb                            = vclass2
           if(tstop(ix,iy,idust)<t_eta)vdrift_turb    = vclass1
           if(tstop(ix,iy,idust)>t_L)vdrift_turb      = vclass3
           
           vdrift_brow  = dsqrt(dSQRT((8.0d0*kb*T/pi)*(mdust(idust)*unit_m+mdust(jdust)*unit_m)/(mdust(idust)*mdust(jdust)*unit_m**2)/unit_v**2)**2)
           vdrift_hydro = dsqrt((q(ix,iy,ivdx(idust))-q(ix,iy,ivdx(jdust)))**2+(q(ix,iy,ivdy(idust))-q(ix,iy,ivdy(jdust)))**2+(q(ix,iy,ivdz(idust))-q(ix,iy,ivdz(jdust)))**2)
           if(turbgrow ==1)  dvij(idust,jdust) = vdrift_turb
           if(browgrow ==1)  dvij(idust,jdust) = dsqrt(dvij(idust,jdust)**2.+vdrift_brow**2.)
           if(driftgrow==1)  dvij(idust,jdust) = dsqrt(dvij(idust,jdust)**2.+vdrift_hydro**2.)
        end do
     end do

     ! Time loop
     do while(time_growth<dt)
        drhodt       = 0.0d0
        niter_growth = niter_growth+1
        do idust=1,ndust
           m1 = mdust(idust)           
           s1 = sdust(idust)
           do jdust=1,idust
              m2  = mdust(jdust)
              s2  = sdust(jdust)
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

                 If(frag_test==1) then

                    Ecol   = 0.5d0*(m1*m2)/(m1+m2)*dvij(idust,jdust)**2
                    Ebr    = (m1+m2)/m_mono*Ebr_mono
                    f_frag = max(min((Ecol-0.1d0*Ebr)/(4.9d0*Ebr),1.0d0),0.0d0)

                 endif
                 dndt = clustered_fraction**2.0*pi*(s1+s2)**2.*dvij(idust,jdust)*u_prim(ix,iy,irhod(idust))*u_prim(ix,iy,irhod(jdust))/m1/m2 ! K n1 n2

              if(kernel_type==1) then ! Constant Kernel

                 p_coag = 1.0d0
                 f_frag = 0.0d0
                 dndt   = u_prim(ix,iy,irhod(idust))*u_prim(ix,iy,irhod(jdust))/m1/m2

              else if(kernel_type==2) then ! Additive Kernel

                 p_coag = 1.0d0
                 f_frag = 0.0d0
                 dndt   = (m1+m2)*u_prim(ix,iy,irhod(idust))*u_prim(ix,iy,irhod(jdust))/m1/m2

              endif

              if(idust==jdust) dndt=dndt/2.0d0


              !#######################################
              !#######################################
              !#######################################
              !#######################################
              !Coagulation/frag conditions
              !#######################################
              !#######################################
              !#######################################
              !#######################################
              if(kernel_type==0) then
                 if(eps_threshold>0.0d0.and.(u_prim(ix,iy,irhod(idust))<eps_threshold*u_prim(ix,iy,irho)))           p_coag=0.0d0
                 if(eps_threshold>0.0d0.and.(u_prim(ix,iy,irhod(jdust))<eps_threshold*u_prim(ix,iy,irho)))           p_coag=0.0d0
                 if(eps_threshold_frag>0.0d0.and.(u_prim(ix,iy,irhod(idust))<eps_threshold_frag*u_prim(ix,iy,irho))) f_frag=0.0d0
                 if(eps_threshold_frag>0.0d0.and.(u_prim(ix,iy,irhod(jdust))<eps_threshold_frag*u_prim(ix,iy,irho))) f_frag=0.0d0
              endif

              ic1 = max(min(floor(dlog((1.0d0-f_frag)*(m1+m2)/massmin)/dlog(eta)+1),ndust),1)
              ic2 = max(min(floor(dlog((1.0d0-f_frag)*(m1+m2)/massmin)/dlog(eta)+1)+1,ndust),1)

              epsilon_mass            = 1.0d0
              if(ic1<ic2)epsilon_mass = min((mdust(ic2)-(1.0d0-f_frag)*(m1+m2))/(mdust(ic2)-mdust(ic1)),1.0d0)

              !#######################################
              !#######################################
              !#######################################
              !#######################################
              ! Coagulation/sticking
              !#######################################
              !#######################################
              !#######################################
              !#######################################

              if(p_coag.ge.0.0d0) then
                 drhodt(idust) = drhodt(idust) - p_coag*(1.0d0-f_frag)*m1*dndt*sticking_efficiency - p_coag*f_frag*m1*dndt
                 drhodt(jdust) = drhodt(jdust) - p_coag*(1.0d0-f_frag)*m2*dndt*sticking_efficiency - p_coag*f_frag*m2*dndt
                 if(f_frag<1.0d0) then
                    drhodt(ic1)  = drhodt(ic1)   +     epsilon_mass * p_coag*(1.0d0-f_frag)*(m1+m2)*dndt*sticking_efficiency
                    drhodt(ic2)  = drhodt(ic2)   + (1.-epsilon_mass)* p_coag*(1.0d0-f_frag)*(m1+m2)*dndt*sticking_efficiency
                 endif
              endif
           
   
              !#######################################
              !#######################################
              !#######################################
              !#######################################
              ! Fragmentation
              !#######################################
              !#######################################
              !#######################################
              !#######################################
              
              if((frag_test==1 .and.f_frag>0.0d0).and.p_coag>0.) then
                 ifrag2 = max(min(floor(dlog(0.1d0*f_frag*(m1+m2)/massmin)/dlog(eta)+1),ic1),1)
                 do kdust = 1,ifrag2
                    drhodt(kdust) = drhodt(kdust)+p_coag*f_frag*(m1+m2)*dndt*redistribute_fragments(kdust,ifrag2)
                 end do
              endif
           end do
        end do

        dt_growth=1d32!(dt-time_growth)
         do idust=1,ndust
            if(abs(drhodt(idust)).ne.0.0d0) then
               dt_growth=min(dt_growth,u_prim(ix,iy,irhod(idust))/abs(drhodt(idust))*CFL_growth)
            endif
         end do
         dt_growth=min(dt_growth,dt-time_growth)

         if(dtcontrol_growth>0.0d0)then
            dt_growth=min(dtcontrol_growth*365.25*3600*24/unit_t,dt-time_growth)
            print *,'time',time_growth,'dt', dt_growth, 'full dt', dt
         endif
         !stop
         do idust=1,ndust
            u_prim(ix,iy,irhod(idust))=max(u_prim(ix,iy,irhod(idust))+drhodt(idust)*dt_growth,u_prim(ix,iy,irho)*dust_ratio_min)
         end do
         time_growth=time_growth+dt_growth
         !print *,time_growth, dt_growth,i
         !print *, 'number of iterations', niter_growth,dt_growth, dt
         
      end do
      if(kernel_type>0) print *, 'number of iterations', niter_growth
      end do
   end do

   
 end subroutine dust_growth



! This files contain the dust growth routines.
subroutine dust_growth_stepinski(verbose)
  use parameters
  use commons
  use units
  use OMP_LIB
  implicit none
  logical :: verbose
  integer :: i,idust,jdust,ipscal,icell,ix,iy

   do iy = first_active_y,last_active_y
      do ix = first_active,last_active
         ! Differential velocity loop
         do idust=1,ndust
           u_prim(ix,iy,idust_pscal(idust,1))=max(u_prim(ix,iy,idust_pscal(idust,1))*(1.0d0+dt/tcoag(ix,iy,idust)),q(ix,iy,irhod(idust))*sminstep/unit_l)
         end do
      end do
   end do

   
 end subroutine dust_growth_stepinski

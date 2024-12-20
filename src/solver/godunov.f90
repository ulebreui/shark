
! This routine performs the predictor operation ! 

! Reference for cylindrical geom :
! The Athena Astrophysical MHD Code in Cylindrical Geometry
! A Proper Discretization of Hydrodynamic Equations in Cylindrical Coordinates for Astrophysical Simulations
subroutine predictor
  use parameters
  use commons
  use units
  use OMP_LIB
  use slope_limiter
  implicit none
  integer  :: i,ix,iy,icell,iymin,iymax,idust,il,ir
  integer  :: ixx,iyy
  integer  :: ivar,idim,ix0,iy0
  real(dp) :: slope_lft,slope_rgt,slope_lim,barotrop,cs_eos,ddxp,ddxm,dx_loc
  real(dp) :: drx,dry,dpx,dpy,dux,duy,dvx,dvy,dwx,dwy,r_rho,u,v,w,p,sr0,sp0,su0,sv0,sw0,dcen,dsgn,dlim,slop,radius_polar


  real(dp) :: dBx_x,dBy_x,dBz_x,Bx,By,Bz,sBx,sBy,sBz,B_norm,total_dust_current_x,total_dust_current_y,total_dust_current_z,derivative_total_dust_current_x,derivative_total_dust_current_y,derivative_total_dust_current_z
  real(dp) :: dB_norm, dBx_over_Bnorm, dBy_over_Bnorm, dBz_over_Bnorm,dbybz,dbybx,dbxbz,deta_a,deta_h,deta_o,dzd
  real(dp) :: total_nd_zd_vx,total_nd_zd_vy,total_nd_zd_vz,vx_derivative_total_nd_zd,vy_derivative_total_nd_zd,vz_derivative_total_nd_zd
  real(dp) :: vz_derivative_total_nd_zd_eta_o,vy_derivative_total_nd_zd_eta_o,vx_derivative_total_nd_zd_eta_h,vy_derivative_total_nd_zd_eta_h,vz_derivative_total_nd_zd_eta_h,vx_derivative_total_nd_zd_eta_a,vy_derivative_total_nd_zd_eta_a,vz_derivative_total_nd_zd_eta_a
  real(dp) :: dHall_i,dni,dne,dB_over_hall,dndzd_over_ni,nd_zd_over_ni
  integer  :: irho_spe,ivx_spe,ivy_spe,ivz_spe,ipscal


  real(dp), dimension(:,:,:),allocatable :: dq
  if(static) return


  ! Initialise to zero
  allocate(dq(1:ncells,1:nvar,1:ndim))

  dq    = 0.0d0
  qp    = 0.0d0
  qm    = 0.0d0

  ! Computes primitive variables
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,il,ir,ivar,ipscal,ix,iy,idim,ix0,iy0,irho_spe,ivx_spe,ivy_spe,ivz_spe,idust,slope_lft,slope_rgt,ddxp,ddxm,drx,dry,dpx,dpy,dux,duy,dvx,dvy,dwx,dwy,r_rho,u,v,w,p,sr0,sp0,su0,sv0,sw0,radius_polar,dBx_x,dBy_x,dBz_x,Bx,By,Bz,sBx,sBy,sBz,B_norm,total_dust_current_x,total_dust_current_y,total_dust_current_z,derivative_total_dust_current_x,derivative_total_dust_current_y,derivative_total_dust_current_z,dB_norm, dBx_over_Bnorm, dBy_over_Bnorm, dBz_over_Bnorm,dbybz,dbybx,dbxbz,deta_a,deta_h,deta_o,dzd,total_nd_zd_vx,total_nd_zd_vy,total_nd_zd_vz,vx_derivative_total_nd_zd,vy_derivative_total_nd_zd,vz_derivative_total_nd_zd,vz_derivative_total_nd_zd_eta_o,vy_derivative_total_nd_zd_eta_o,vx_derivative_total_nd_zd_eta_h,vy_derivative_total_nd_zd_eta_h,vz_derivative_total_nd_zd_eta_h,vx_derivative_total_nd_zd_eta_a,vy_derivative_total_nd_zd_eta_a,vz_derivative_total_nd_zd_eta_a,nd_zd_over_ni,dHall_i,dni,dne,dB_over_hall,dndzd_over_ni)
  drx   = 0.0d0 
  dpx   = 0.0d0  
  dry   = 0.0d0 
  dpy   = 0.0d0
  duy   = 0.0d0
  dvx   = 0.0d0
  dvy   = 0.0d0
  dwx   = 0.0d0
  dwy   = 0.0d0
  u     = 0.0d0
  v     = 0.0d0
  w     = 0.0d0
  dBx_x = 0.0d0
  dBy_x = 0.0d0
  dBz_x = 0.0d0
  Bx    = 0.0d0
  By    = 0.0d0
  Bz    = 0.0d0
  sBx   = 0.0d0
  sBy   = 0.0d0
  sBz   = 0.0d0
  B_norm = 0.0d0
  dB_norm = 0.0d0
  dBx_over_Bnorm = 0.0d0
  dBy_over_Bnorm = 0.0d0
  dBz_over_Bnorm = 0.0d0
  dbybz = 0.0d0
  dbybx = 0.0d0
  dbxbz = 0.0d0
  deta_a = 0.0d0
  deta_h = 0.0d0
  deta_o = 0.0d0
  dzd = 0.0d0
  vz_derivative_total_nd_zd_eta_o = 0.0d0
  vy_derivative_total_nd_zd_eta_o = 0.0d0
  vx_derivative_total_nd_zd_eta_h = 0.0d0
  vy_derivative_total_nd_zd_eta_h = 0.0d0
  vz_derivative_total_nd_zd_eta_h = 0.0d0
  vx_derivative_total_nd_zd_eta_a = 0.0d0
  vy_derivative_total_nd_zd_eta_a = 0.0d0
  vz_derivative_total_nd_zd_eta_a = 0.0d0




  !$OMP DO
  do i=1,ncells
    if(active_cell_predictor(i)==1) then
      ix=ixx(i)
      iy=iyy(i)

      radius_polar=1.0d0     
#if GEOM==2
      radius_polar=radii_c(i)
#endif
#if GEOM==4
      radius_polar=1.0d0
#endif  
    if(slope_type>0) then
    do ivar = 1, nvar
            il = icell(ix-1,iy)
            ir = icell(ix+1,iy)
            dq(i,ivar,1) = slope_limit(2.0d0*(q(i,ivar) - q(il,ivar))/(dx(i,1)+dx(il,1)),2.0d0*(q(ir,ivar) - q(i,ivar))/(dx(ir,1)+dx(i,1)))
#if NY>1
            il = icell(ix,iy-1)
            ir = icell(ix,iy+1)
            dq(i,ivar,2) = slope_limit(2.0d0*(q(i,ivar) - q(il,ivar))/(dx(i,2)+dx(il,2)),2.0d0*(q(ir,ivar) - q(i,ivar))/(dx(ir,2)+dx(i,2)))/radius_polar
#endif   
    end do
    endif

      r_rho = q(i,irho)
      u     = q(i,ivx)
      dux   = dq(i,ivx,1)
      drx   = dq(i,irho,1)
      w     = q(i,ivz)
      dwx   = dq(i,ivz,1)
      v     = q(i,ivy)
      dvx   = dq(i,ivy,1)
      p     = q(i,iP)
      dPx   = dq(i,iP,1)
      if(iso_cs==1) then
       P     = q(i,irho)*cs(i)**2
       dPx   = dq(i,irho,1)*cs(i)**2
      endif 
#if NY>1
      duy   = dq(i,ivx,2)
      dry   = dq(i,irho,2)
      dvy   = dq(i,ivy,2)
      dpy   = dq(i,iP,2)
      dwy   = dq(i,ivz,2)
      if(iso_cs==1) then
        dPy   = dq(i,irho,2)*cs(i)**2
      endif
#endif
      
    
    if(force_kick) then
        u    = u + force(i,1)*half*dt
        v    = v + force(i,2)*half*dt
        w    = w + force(i,3)*half*dt
    endif


#if MHD==1
      Bx       = q(i,iBx)
      By       = q(i,iBy)
      Bz       = q(i,iBz)
      dBx_x    = dq(i,iBx,1)
      dBy_x    = dq(i,iBy,1)
      dBz_x    = dq(i,iBz,1)
      B_norm = dsqrt(Bx**2+By**2+Bz**2)
      dB_norm = (2*Bx*dBx_x+2*By*dBy_x+2*Bz*dBz_x)/(2*SQRT(Bx**2+By**2+Bz**2))
      dBx_over_Bnorm = (dBx_x*B_norm-Bx*dB_norm)/B_norm**2
      dBy_over_Bnorm = (dBy_x*B_norm-By*dB_norm)/B_norm**2
      dBz_over_Bnorm = (dBz_x*B_norm-Bz*dB_norm)/B_norm**2
      dbybz = dBy_over_Bnorm*Bz/B_norm + By/B_norm*dBz_over_Bnorm
      dbybx = dBy_over_Bnorm*Bx/B_norm + By/B_norm*dBx_over_Bnorm
      dbxbz = dBx_over_Bnorm*Bz/B_norm + Bx/B_norm*dBz_over_Bnorm

      if (dusty_nonideal_MHD) then
          deta_o = slope_limit(2.0d0*(eta_o(i) - eta_o(il))/(dx(i,1)+dx(il,1)),2.0d0*(eta_o(ir) - eta_o(i))/(dx(ir,1)+dx(i,1)))
          deta_h = slope_limit(2.0d0*(eta_h(i) - eta_h(il))/(dx(i,1)+dx(il,1)),2.0d0*(eta_h(ir) - eta_h(i))/(dx(ir,1)+dx(i,1)))
          deta_a = slope_limit(2.0d0*(eta_a(i) - eta_a(il))/(dx(i,1)+dx(il,1)),2.0d0*(eta_a(ir) - eta_a(i))/(dx(ir,1)+dx(i,1)))
      endif

       if (dusty_nonideal_MHD_no_electron) then
          dHall_i = slope_limit(2.0d0*(Hall_i(i) - Hall_i(il))/(dx(i,1)+dx(il,1)),2.0d0*(Hall_i(ir) - Hall_i(i))/(dx(ir,1)+dx(i,1)))
          dni = slope_limit(2.0d0*(ni(i) - ni(il))/(dx(i,1)+dx(il,1)),2.0d0*(ni(ir) - ni(i))/(dx(ir,1)+dx(i,1)))
          dne = slope_limit(2.0d0*(ne(i) - ne(il))/(dx(i,1)+dx(il,1)),2.0d0*(ne(ir) - ne(i))/(dx(ir,1)+dx(i,1)))
          if(electrons .eqv. .false.) dne = 0.0d0  

          idust = i_coupled_species

          dzd = slope_limit(2.0d0*(zd(i,idust) - zd(il,idust))/(dx(i,1)+dx(il,1)),2.0d0*(zd(ir,idust) - zd(i,idust))/(dx(ir,1)+dx(i,1)))
          dB_over_hall = (dB_norm*Hall_i(i)-B_norm*dHall_i)/Hall_i(i)**2 !Initialize Hall_i =/ 0 to avoid Nan at first timestep
          dndzd_over_ni = ((dq(i,irhod(idust),1)/mdust(i,idust)*zd(i,idust)+q(i,irhod(idust))/mdust(i,idust)*dzd)*(ni(i)+ne(i)) - (q(i,irhod(idust))/mdust(i,idust)*zd(i,idust))*(dni+dne))/(ni(i)+ne(i))**2
          ! print *,'dzd',dzd
          ! print *,'dni',dni
          ! print *,'dne',dne
          ! print *,'dHall_i',dHall_i
          ! print *,'dB_over_hall',dB_over_hall



      endif
           


#endif

!Cartesian geometry
#if GE0M==0
      sr0    = -u*drx-v*dry - (dux+dvy)*r_rho
      sp0    = -u*dpx-v*dpy - (dux+dvy)*gamma*p
      su0    = -u*dux-v*duy - (dpx        )/r_rho
      sv0    = -u*dvx-v*dvy - (dpy        )/r_rho
      sw0    = -u*dwx-v*dwy
#endif

!Spherical geometry
#if GEOM==1
      sr0    = -u*drx-v*dry - (dux+dvy)*r_rho     - r_rho*u/radii_c(i)
      sp0    = -u*dpx-v*dpy - (dux+dvy)*gamma*p
      su0    = -u*dux-v*duy - (dpx        )/r_rho
      sv0    = -u*dvx-v*dvy - (dpy        )/r_rho
      sw0    = -u*dwx-v*dwy
#endif       

!Disk (face-on) geometry
#if GEOM==2
      sr0    = -u*drx-v*dry - (dux+dvy)*r_rho      - r_rho*u    / radius_polar
      sp0    = -u*dpx-v*dpy - (dux+dvy)*gamma*p    - gamma*p*u  / radius_polar 
      su0    = -u*dux-v*duy - (dpx        )/r_rho  + (v**2.)    / radius_polar
      sv0    = -u*dvx-v*dvy - (dpy        )/r_rho  - u*v        / radius_polar
      sw0    = -u*dwx-v*dwy
#endif


!Disk (edge-on) geometry 
#if GEOM==4
      sr0    = -u*drx-v*dry - (dux+dvy)*r_rho      - r_rho*u    / radii_c(i)
      sp0    = -u*dpx-v*dpy - (dux+dvy)*gamma*p    - gamma*p*u  / radii_c(i)
      su0    = -u*dux-v*duy - (dpx        )/r_rho  + (w**2.)    / radii_c(i)
      sv0    = -u*dvx-v*dvy - (dpy        )/r_rho  
      sw0    = -u*dwx-v*dwy - u*w        / radii_c(i)
#endif

#if MHD==1
!If no dust, couple B to the gas
#if NDUST==0
    sBy   = -u*dBy_x + Bx*dvx - By*dux
    sBz   = -u*dBz_x + Bx*dwx - Bz*dux
    su0   = su0 - By/r_rho*dBy_x - Bz/r_rho*dBz_x - Bx*dBx_x/r_rho + 2.0d0*Bx*dBx_x/r_rho 
    sv0   = sv0 + Bx*dBy_x/r_rho
    sw0   = sw0 + Bx*dBz_x/r_rho
#endif
#endif

#if GRAVITY==1
#if NY==1
     !Gravity source term in 1D
     su0 = su0 + (-Mc(i)/(radii_c(i)**2.+(l_soft/unit_l)**2.))
#endif
#endif      
    
    ! if(force_kick) then
    !     su0    = su0 + force(i,1)
    !     sv0    = sv0 + force(i,2)
    !     sw0    = sw0 + force(i,3)
    ! endif


    !direction x
    dx_loc=dx(i,1)
    qm(i,irho,1)   = r_rho     + half*dt*sr0  + half*drx   *dx_loc
    qm(i,ivx,1)    = u         + half*dt*su0  + half*dux   *dx_loc
    qm(i,iP,1)     = max(p     + half*dt*sP0  + half*dpx   *dx_loc,smallP)
    qm(i,ivy,1)    = v         + half*dt*sv0  + half*dvx   *dx_loc
    qm(i,ivz,1)    = w         + half*dt*sw0  + half*dwx   *dx_loc 
    qp(i,irho,1)   =r_rho      + half*dt*sr0  - half*drx   *dx_loc
    qp(i,ivx,1)    = u         + half*dt*su0  - half*dux   *dx_loc
    qp(i,iP,1)     = max(p     + half*dt*sP0  - half*dpx   *dx_loc,smallP)
    qp(i,ivy,1)    = v         + half*dt*sv0  - half*dvx   *dx_loc
    qp(i,ivz,1)    = w         + half*dt*sw0  - half*dwx   *dx_loc

    ! direction y
#if NY>1
    dx_loc=radius_polar*dx(i,2)
    qm(i,irho,2)   = r_rho     + half*dt*sr0  + half*dry   *dx_loc
    qm(i,ivx,2)    = u         + half*dt*su0  + half*duy   *dx_loc
    qm(i,iP,2)     = max(p     + half*dt*sP0  + half*dpy   *dx_loc,smallP)
    qm(i,ivy,2)    = v         + half*dt*sv0  + half*dvy   *dx_loc
    qm(i,ivz,2)    = w         + half*dt*sw0  + half*dwy   *dx_loc 
    qp(i,irho,2)   = r_rho     + half*dt*sr0  - half*dry   *dx_loc
    qp(i,ivx,2)    = u         + half*dt*su0  - half*duy   *dx_loc
    qp(i,iP,2)     = max(p     + half*dt*sP0  - half*dpy   *dx_loc,smallP)
    qp(i,ivy,2)    = v         + half*dt*sv0  - half*dvy   *dx_loc
    qp(i,ivz,2)    = w         + half*dt*sw0  - half*dwy   *dx_loc
 
#endif

    ! Dust terms: same remark as for the gas
#if NDUST>0

      total_dust_current_x = 0.0d0
      total_dust_current_y = 0.0d0
      total_dust_current_z = 0.0d0
      derivative_total_dust_current_x = 0.0d0
      derivative_total_dust_current_y = 0.0d0
      derivative_total_dust_current_z = 0.0d0

      total_nd_zd_vx = 0.0d0
      total_nd_zd_vy = 0.0d0
      total_nd_zd_vz = 0.0d0

      vx_derivative_total_nd_zd = 0.0d0
      vy_derivative_total_nd_zd = 0.0d0
      vz_derivative_total_nd_zd = 0.0d0


      if (dusty_nonideal_MHD) then    
          do idust=1,ndust  !!Quantities needed for non-ideal dusty MHD (1D only) 
            irho_spe = irhod(idust)
            ivx_spe  = ivdx(idust)
            ivy_spe  = ivdy(idust) 
            ivz_spe  = ivdz(idust)
            r_rho = q(i,irho_spe)
            u     = q(i,ivx_spe)
            v     = q(i,ivy_spe)
            w     = q(i,ivz_spe)
            drx   = dq(i,irho_spe,1)
            dux   = dq(i,ivx_spe,1)
            dvx   = dq(i,ivy_spe,1)
            dwx   = dq(i,ivz_spe,1)

            dzd = slope_limit(2.0d0*(zd(i,idust) - zd(il,idust))/(dx(i,1)+dx(il,1)),2.0d0*(zd(ir,idust) - zd(i,idust))/(dx(ir,1)+dx(i,1)))


            !!RQ: zd et eta are space dependent and as such corresponding derivatives should be included. The current state of this routine is approximate. 


            total_dust_current_x = total_dust_current_x + (1/clight)*(e_el_stat)*(r_rho/mdust(i,idust))*zd(i,idust)*(u-q(i,ivx)) !Total (relative to neutral velocity) dust current
            total_dust_current_y = total_dust_current_y + (1/clight)*(e_el_stat)*(r_rho/mdust(i,idust))*zd(i,idust)*(v-q(i,ivy))
            total_dust_current_z = total_dust_current_z + (1/clight)*(e_el_stat)*(r_rho/mdust(i,idust))*zd(i,idust)*(w-q(i,ivz))

            derivative_total_dust_current_x = derivative_total_dust_current_x + (1/clight)*(e_el_stat)*(r_rho/mdust(i,idust))*zd(i,idust)*(dux-dq(i,ivx,1)) 
            derivative_total_dust_current_y = derivative_total_dust_current_y + (1/clight)*(e_el_stat)*(r_rho/mdust(i,idust))*zd(i,idust)*(dvx-dq(i,ivy,1))
            derivative_total_dust_current_z = derivative_total_dust_current_z + (1/clight)*(e_el_stat)*(r_rho/mdust(i,idust))*zd(i,idust)*(dwx-dq(i,ivz,1))

            total_nd_zd_vx = total_nd_zd_vx + (r_rho/mdust(i,idust))*zd(i,idust)*(e_el_stat)*(u-q(i,ivx))
            total_nd_zd_vy = total_nd_zd_vy + (r_rho/mdust(i,idust))*zd(i,idust)*(e_el_stat)*(v-q(i,ivy))
            total_nd_zd_vz = total_nd_zd_vz + (r_rho/mdust(i,idust))*zd(i,idust)*(e_el_stat)*(w-q(i,ivz))

            vx_derivative_total_nd_zd = vx_derivative_total_nd_zd + ( (e_el_stat)*dzd*(r_rho/mdust(i,idust)) + (e_el_stat)*zd(i,idust)*drx/mdust(i,idust) )*(u-q(i,ivx)) 
            vy_derivative_total_nd_zd = vy_derivative_total_nd_zd + ( (e_el_stat)*dzd*(r_rho/mdust(i,idust)) + (e_el_stat)*zd(i,idust)*drx/mdust(i,idust) )*(v-q(i,ivy)) 
            vz_derivative_total_nd_zd = vz_derivative_total_nd_zd + ( (e_el_stat)*dzd*(r_rho/mdust(i,idust)) + (e_el_stat)*zd(i,idust)*drx/mdust(i,idust) )*(w-q(i,ivz))

          end do

          vy_derivative_total_nd_zd_eta_o = (1/clight)*(deta_o*total_nd_zd_vy + eta_o(i)*vy_derivative_total_nd_zd)
          vz_derivative_total_nd_zd_eta_o = (1/clight)*(deta_o*total_nd_zd_vz + eta_o(i)*vz_derivative_total_nd_zd)

          vx_derivative_total_nd_zd_eta_h = (1/clight)*(deta_h*total_nd_zd_vx + eta_h(i)*vx_derivative_total_nd_zd)
          vy_derivative_total_nd_zd_eta_h = (1/clight)*(deta_h*total_nd_zd_vy + eta_h(i)*vy_derivative_total_nd_zd)
          vz_derivative_total_nd_zd_eta_h = (1/clight)*(deta_h*total_nd_zd_vz + eta_h(i)*vz_derivative_total_nd_zd)


          vx_derivative_total_nd_zd_eta_a = (1/clight)*(deta_a*total_nd_zd_vx + eta_a(i)*vx_derivative_total_nd_zd)
          vy_derivative_total_nd_zd_eta_a = (1/clight)*(deta_a*total_nd_zd_vy + eta_a(i)*vy_derivative_total_nd_zd)
          vz_derivative_total_nd_zd_eta_a = (1/clight)*(deta_a*total_nd_zd_vz + eta_a(i)*vz_derivative_total_nd_zd)

      endif


      do idust=1,ndust
        irho_spe = irhod(idust)
        ivx_spe  = ivdx(idust)
        ivy_spe  = ivdy(idust) 
        ivz_spe  = ivdz(idust)

        r_rho = q(i,irho_spe)
        u     = q(i,ivx_spe)
        dux   = dq(i,ivx_spe,1)
        drx   = dq(i,irho_spe,1)
        v     = q(i,ivy_spe)
        dvx   = dq(i,ivy_spe,1)
        w     = q(i,ivz_spe)
        dwx   = dq(i,ivz_spe,1)
#if NY>1   
        dry   = dq(i,irho_spe,2)   
        duy   = dq(i,ivx_spe,2)  
        dvy   = dq(i,ivy_spe,2)
        dwy   = dq(i,ivz_spe,2)
#endif
    if(force_kick) then
        u    = u + force_dust(i,1,idust)*half*dt
        v    = v + force_dust(i,2,idust)*half*dt
        w    = w + force_dust(i,3,idust)*half*dt
    endif

        sr0    = -u*drx-v*dry - (dux+dvy)*r_rho
#if GEOM==1
        !Spherical geometry source term
        sr0    = sr0-r_rho*u/radii_c(i)
#endif

        su0    = -u*dux-v*duy       
        sv0    = -u*dvx-v*dvy 
        sw0    = -u*dwx-v*dwy

#if GEOM==2
        !Polar geometry source terms
        sr0    = sr0 - r_rho*u/radius_polar
        su0    = su0 + v**2.  /radius_polar
        sv0    = sv0 - u*v    /radius_polar
#endif

#if GEOM==4
        !Polar geometry source terms -- TODO add the missing source terms
        sr0    = sr0 - r_rho*u/radii_c(i)
        su0    = su0 + w**2.  /radii_c(i)
        sw0    = sw0 - u*w    /radii_c(i)
#endif
#if MHD==1
!Valid only for GEOM==0 since the s are overwritten
    if (idust==i_coupled_species) then
        if(dusty_nonideal_MHD .eqv. .false.) then
        if(dusty_nonideal_MHD_no_electron .eqv. .false.) then
        sBy   = -u*dBy_x + Bx*dvx - By*dux
        sBz   = -u*dBz_x + Bx*dwx - Bz*dux

        su0   = -u*dux-v*duy - By/r_rho*dBy_x - Bz/r_rho*dBz_x - Bx*dBx_x/r_rho + 2.0d0*Bx*dBx_x
        sv0   = -u*dvx-v*dvy + Bx*dBy_x/r_rho
        sw0   = -u*dwx-v*dwy + Bx*dBz_x/r_rho



        endif
        endif

    endif

    if (dusty_nonideal_MHD .or. dusty_nonideal_MHD_no_electron) then 

        su0   = -u*dux-v*duy !Lorentz forces are source terms --> do not appear here
        sv0   = -u*dvx-v*dvy
        sw0   = -u*dwx-v*dwy 


        ! su0   = -u*dux-v*duy - By/r_rho*dBy_x - Bz/r_rho*dBz_x - Bx*dBx_x/r_rho + 2.0d0*Bx*dBx_x
        ! sv0   = -u*dvx-v*dvy + Bx*dBy_x/r_rho
        ! sw0   = -u*dwx-v*dwy + Bx*dBz_x/r_rho
    endif



#endif
#if GRAVITY==1
#if NY==1
        !Gravity source term in 1D
        su0=su0+(-Mc(i)/(radii_c(i)**2.+(l_soft/unit_l)**2.))
#endif
#endif   

    ! if(force_kick) then
    !     su0    = su0 + force_dust(i,1,idust)
    !     sv0    = sv0 + force_dust(i,2,idust)
    !     sw0    = sw0 + force_dust(i,3,idust)
    ! endif

    !Direction x
    dx_loc=dx(i,1)
    qm(i,irho_spe,1)   = r_rho     + half*dt*sr0  + half*drx   *dx_loc
    qm(i,ivx_spe,1)    = u         + half*dt*su0  + half*dux   *dx_loc
    qm(i,ivy_spe,1)    = v         + half*dt*sv0  + half*dvx   *dx_loc
    qm(i,ivz_spe,1)    = w         + half*dt*sw0  + half*dwx   *dx_loc 
    qp(i,irho_spe,1)   = r_rho     + half*dt*sr0  - half*drx   *dx_loc
    qp(i,ivx_spe,1)    = u         + half*dt*su0  - half*dux   *dx_loc
    qp(i,ivy_spe,1)    = v         + half*dt*sv0  - half*dvx   *dx_loc
    qp(i,ivz_spe,1)    = w         + half*dt*sw0  - half*dwx   *dx_loc
    
    !direction y
#if NY>1
    dx_loc=radius_polar*dx(i,2)
    qm(i,irho_spe,2)   = r_rho     + half*dt*sr0  + half*dry   *dx_loc
    qm(i,ivx_spe,2)    = u         + half*dt*su0  + half*duy   *dx_loc
    qm(i,ivy_spe,2)    = v         + half*dt*sv0  + half*dvy   *dx_loc
    qm(i,ivz_spe,2)    = w         + half*dt*sw0  + half*dwy   *dx_loc 
    qp(i,irho_spe,2)   = r_rho     + half*dt*sr0  - half*dry   *dx_loc
    qp(i,ivx_spe,2)    = u         + half*dt*su0  - half*duy   *dx_loc
    qp(i,ivy_spe,2)    = v         + half*dt*sv0  - half*dvy   *dx_loc
    qp(i,ivz_spe,2)    = w         + half*dt*sw0  - half*dwy   *dx_loc

#endif

#if NDUSTPSCAL>0
    do ipscal= 1, ndustpscal
        r_rho = q(i,idust_pscal(idust,ipscal))

        drx   = dq(i,idust_pscal(idust,ipscal),1)
        dry   = dq(i,idust_pscal(idust,ipscal),2)

        sr0   = -u*drx-v*dry - (dux+dvy)*r_rho
#if GEOM==2
        !Polar geometry source terms
        sr0    = sr0 - r_rho*u/radius_polar
#endif
#if GEOM==4
        !Polar geometry source terms
        sr0    = sr0 - r_rho*u/radii_c(i)
#endif
        ! Direction x

        dx_loc=dx(i,1)
        qm(i,idust_pscal(idust,ipscal),1)   = r_rho + half*dt*sr0  + half*drx   * dx_loc
        qp(i,idust_pscal(idust,ipscal),1)   = r_rho + half*dt*sr0  - half*drx   * dx_loc
        
        ! Direction y
#if NY>1
       dx_loc=radius_polar*dx(i,2)
       qm(i,idust_pscal(idust,ipscal),2)   = r_rho + half*dt*sr0  + half*dry   * dx_loc
       qp(i,idust_pscal(idust,ipscal),2)   = r_rho + half*dt*sr0  - half*dry   * dx_loc
#endif
    end do
#endif
      end do !dust loop
#endif

#if MHD==1


    if (dusty_nonideal_MHD) then
        !There are more conservative terms in the induction equation (a lot more!)
        !This part accounts for multiple dust species


        sBy = q(i,ivx)*dBy_x + Bx*dq(i,ivy,1) - By*dq(i,ivx,1) !"ideal term" --> features gas velocity and not dust velocity
        ! sBy = sBy - eta_o(i)*derivative_total_dust_current_z - vz_derivative_total_nd_zd_eta_o!conservative term associated to Ohm resistivity
        ! sBy = sBy + eta_H(i)*Bx/B_norm*derivative_total_dust_current_y + eta_H(i)*total_dust_current_y*dBx_over_Bnorm + bx*vy_derivative_total_nd_zd_eta_h + eta_H(i)*By/B_norm*derivative_total_dust_current_x + eta_H(i)*total_dust_current_x*dBy_over_Bnorm + by*vx_derivative_total_nd_zd_eta_h !Hall related
        ! sBy = sBy + eta_a(i)*dbybz*total_dust_current_y + eta_a(i)*By/B_norm*Bz/B_norm*derivative_total_dust_current_y + by*bz*vy_derivative_total_nd_zd_eta_a - 2*eta_a(i)*(Bx/B_norm)*dBx_over_Bnorm*total_dust_current_z - eta_a(i)*(Bx/B_norm)**2*derivative_total_dust_current_z  - bx**2*vz_derivative_total_nd_zd_eta_a - 2*eta_a(i)*By/B_norm*dBy_over_Bnorm*total_dust_current_z - eta_a(i)*(By/B_norm)**2*derivative_total_dust_current_z - by**2*vz_derivative_total_nd_zd_eta_a + eta_a(i)*dbxbz*total_dust_current_x + eta_a(i)*Bx/B_norm*Bz/B_norm*derivative_total_dust_current_x + bx*bz*vx_derivative_total_nd_zd_eta_a

        !sBz --> substitute y and z indices + change signs for AD and Ohm terms and second Hall term
        sBz = -q(i,ivx)*dBz_x + Bx*dq(i,ivz,1) - Bz*dq(i,ivx,1) !"ideal term" --> features gas velocity and not dust velocity
        ! sBz = sBz + eta_o(i)*derivative_total_dust_current_y + vy_derivative_total_nd_zd_eta_o!conservative term associated to Ohm resistivity
        ! sBz = sBz + eta_H(i)*Bx/B_norm*derivative_total_dust_current_z + eta_H(i)*total_dust_current_z*dBx_over_Bnorm + bx*vz_derivative_total_nd_zd_eta_h - eta_H(i)*Bz/B_norm*derivative_total_dust_current_x - eta_H(i)*total_dust_current_x*dBz_over_Bnorm - bz*vx_derivative_total_nd_zd_eta_h!Hall related
        ! sBz = sBz - eta_a(i)*dbybz*total_dust_current_z - eta_a(i)*By/B_norm*Bz/B_norm*derivative_total_dust_current_z + bz*by*vz_derivative_total_nd_zd_eta_a + 2*eta_a(i)*(Bx/B_norm)*dBx_over_Bnorm*total_dust_current_y + eta_a(i)*(Bx/B_norm)**2*derivative_total_dust_current_y + bx**2*vy_derivative_total_nd_zd_eta_a + 2*eta_a(i)*Bz/B_norm*dBz_over_Bnorm*total_dust_current_y + eta_a(i)*(Bz/B_norm)**2*derivative_total_dust_current_y + bz**2*vy_derivative_total_nd_zd_eta_a - eta_a(i)*dbybx*total_dust_current_x - eta_a(i)*Bx/B_norm*By/B_norm*derivative_total_dust_current_x - bx*by*vx_derivative_total_nd_zd_eta_a
    
    endif


    if (dusty_nonideal_MHD_no_electron) then !A single grain only

        sBy   = -u*dBy_x + Bx*dvx - By*dux
        sBz   = -u*dBz_x + Bx*dwx - Bz*dux

        idust = i_coupled_species
        nd_zd_over_ni = (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/(ni(i)+ne(i))

        ! print *, 'nd_zd', nd_zd_over_ni
        ! print *, 'dnd_zd', dndzd_over_ni
        !   print *,'rhod',q(i,irhod(:))
        !   print *,'zd',zd(i,:)
        !   print *,'ni(:)',ni(i)
        !   print *,'ne(:)',ne(i)
           ! print *,'Hall',Hall_i(i)
           ! print *,'dB_Hall',dB_over_hall
           ! print *,'dHall',dHall_i
           ! print *,'B_norm',B_norm
           ! print *,'dB_norm',dB_norm

           ! print *,'u',u
           ! print *,'v',v
           ! print *,'w',w
           ! print *,'dux',dux
           ! print *,'dvx',dvx
           ! print *,'dq(i,ivz,1)',dq(i,ivz,1)
           ! print *,'dq(i,ivz,1)',dq(i,ivz,1)



        sBy = dndzd_over_ni*u*By + nd_zd_over_ni*(u*dBy_x + By*dux) - Bx*v*dndzd_over_ni - nd_zd_over_ni*Bx*dvx !"ideal term" 
        sBy = sBy + dB_over_hall*(nd_zd_over_ni*w + q(i,ivz)) + B_norm/Hall_i(i)*(dndzd_over_ni*w + nd_zd_over_ni*dwx + dq(i,ivz,1)) 

        ! !!change signs for extra non-ideal term
        sBz = dndzd_over_ni*u*Bz + nd_zd_over_ni*(u*dBz_x + Bz*dux) - Bx*w*dndzd_over_ni - nd_zd_over_ni*Bx*dwx !"ideal term" 
        sBz = sBz - dB_over_hall*(nd_zd_over_ni*v + q(i,ivy)) - B_norm/Hall_i(i)*(dndzd_over_ni*v + nd_zd_over_ni*dvx + dq(i,ivy,1))

        ! nd_zd_over_ni = (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/(ni(i)+ne(i))
        !sBy = nd_zd_over_ni*(u*dBy_x + By*dux) - nd_zd_over_ni*Bx*dvx !"ideal term" 
        ! !sBy = sBy + dB_over_hall*(nd_zd_over_ni*w + q(i,ivz)) + B_norm/Hall_i(i)*(dndzd_over_ni*w + nd_zd_over_ni*dwx + dq(i,ivz,1)) 

        ! !!change signs for extra non-ideal term
        !sBz =  nd_zd_over_ni*(u*dBz_x + Bz*dux) - nd_zd_over_ni*Bx*dwx !"ideal term" 
        ! !sBz = sBz - dB_over_hall*(nd_zd_over_ni*v + q(i,ivy)) - B_norm/Hall_i(i)*(dndzd_over_ni*v + nd_zd_over_ni*dvx + dq(i,ivy,1))  

    endif




    do idim=1,ndim
        ddxm = dx(i,idim)
        ddxp = dx(i,idim)
#if GEOM==1 
        !THis is a correction for log grids dx_l =!= dx_r
#if NY==1
        ddxm=dx_l_cell(i)
        ddxp=dx_r_cell(i)
#endif
#endif

        qm(i,iBx,idim)    = q(i,iBx)     + half*dt*sBx + half*dq(i,iBx,idim)  * ddxm 
        qm(i,iBy,idim)    = q(i,iBy)     + half*dt*sBy + half*dq(i,iBy,idim)  * ddxm 
        qm(i,iBz,idim)    = q(i,iBz)     + half*dt*sBz + half*dq(i,iBz,idim)  * ddxm 
        qp(i,iBx,idim)    = q(i,iBx)     + half*dt*sBx - half*dq(i,iBx,idim)  * ddxp
        qp(i,iBy,idim)    = q(i,iBy)     + half*dt*sBy - half*dq(i,iBy,idim)  * ddxp
        qp(i,iBz,idim)    = q(i,iBz)     + half*dt*sBz - half*dq(i,iBz,idim)  * ddxp

    end do
#endif
! We recompute the thermal pressure if we don't solve for the NRJ equation
if(iso_cs==1) then
    do idim =1,ndim
        qp(i,iP,idim)   = cs(i)**2*qp(i,irho,idim)
        qm(i,iP,idim)   = cs(i)**2*qm(i,irho,idim)
    end do
end if
if(non_standard_eos==1) then
    do idim=1,ndim
        qp(i,iP,idim)   = cs_eos(barotrop(qp(i,irho,idim)))**2*qp(i,irho,idim)
        qm(i,iP,idim)   = cs_eos(barotrop(qm(i,irho,idim)))**2*qm(i,irho,idim)
   end do
end if
end if !End of very first if loop
end do !End of very first do loop (i)
  !$OMP END DO

  !$OMP END PARALLEL
  deallocate(dq)

#if GRAVITY==1  
  call mtot
#endif 

end subroutine predictor


! This is the Riemmann solver : llf + Godunov scheme
subroutine add_delta_u
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust,ivar,ix,iy,il,ily,icell
  integer :: ixx,iyy

  real(dp), dimension(1:nvar) :: qleft,qright,flx
  real(dp) :: csr,csl,barotrop,cs_eos


  if(static) return

  flux       = 0.0d0
! #if GEOM==2
!     allocate(source_loc(1:ncells,1:nvar))
!     source_loc=0.0d0
! #endif

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,ivar,ix,iy,il,ily,qleft,qright,flx,csr,csl)
  !$OMP DO
  do i = 1, ncells
    if(active_cell_predictor(i)==1) then
      ix=ixx(i)
      iy=iyy(i)

      !First direction x
      il = icell(ix-1,iy)
      do ivar=1,nvar
            qleft(ivar)  = qm(il,ivar,1)
            qright(ivar) = qp(i,ivar,1)
            flx(ivar)    = 0.0d0
      end do

      if(iso_cs==1) then
            csl=cs(il)
            csr=cs(i)
      else if(non_standard_eos==1) then
            csl=cs_eos(barotrop(qleft(irho)))
            csr=cs_eos(barotrop(qright(irho)))
      else      
          csl=sqrt(gamma*qleft(iP)/qleft(irho))
          csr=sqrt(gamma*qright(iP)/qright(irho))
      endif

      call solve_wrapper(qleft,qright,flx,csl,csr,1,i)

      do ivar=1,nvar
            flux(i,ivar,1)=flx(ivar) 
      end do

    ! Then direction y
#if NY>1
    il = icell(ix,iy-1)

        do ivar=1,nvar
            qleft(ivar)  = qm(il,ivar,2)
            qright(ivar) = qp(i,ivar, 2)
            flx(ivar)    = 0.0d0
        end do

        if(iso_cs==1) then
            csl = cs(il)
            csr = cs(i)
        else if(non_standard_eos==1) then

            csl = cs_eos(barotrop(qleft(irho)))
            csr = cs_eos(barotrop(qright(irho)))

        else      

          csl = sqrt(gamma*qleft(iP)/qleft(irho))
          csr = sqrt(gamma*qright(iP)/qright(irho))
          
        endif

        call solve_wrapper(qleft,qright,flx,csl,csr,2,i)

        do ivar=1,nvar
            flux(i,ivar,2)=flx(ivar) 
        end do
#endif

    end if
  end do
  !$OMP END DO
  !$OMP BARRIER

  !Update state vector 
  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
        ix = ixx(i)
        iy = iyy(i)
#if NY==1
        il  = icell(ix+1,iy)
        ! 1D case
        do ivar = 1, nvar
            u_prim(i,ivar)=u_prim(i,ivar)+(flux(i,ivar,1)*surf(i,1)-flux(il,ivar,1)*surf(il,1))/vol(i)*dt
        end do
#endif
#if NY>1
        ! 2D case : we do both fluxes at the same time to be more efficient
        il  = icell(ix+1,iy)
        ily = icell(ix,iy+1)
        do ivar = 1, nvar
              u_prim(i,ivar)=u_prim(i,ivar) + (flux(i,ivar,1)*surf(i,1)-flux(il,ivar,1)*surf(il,1))  /vol(i)*dt&
              &                             + (flux(i,ivar,2)*surf(i,2)-flux(ily,ivar,2)*surf(ily,2))/vol(i)*dt
        end do
#endif

      endif
    end do

  !$OMP END DO
  !$OMP END PARALLEL


! #if GEOM==2
!   deallocate(source_loc)
! #endif

end subroutine add_delta_u







subroutine solve_wrapper(qleft,qright,flx,csl,csr,idim,i)
 use hydro_solvers
 use parameters
 use commons

 implicit none

 real(dp),dimension(1:nvar),intent(in) :: qright,qleft
 real(dp),dimension(1:nvar),intent(inout) :: flx
 real(dp) :: csl,csr

 integer  :: idim,i


 ! First the gas
#if SOLVER==0        
    call solver_llf(qleft,qright,flx,csl,csr,idim,i)
#endif
#if SOLVER==1       
    call solver_hll(qleft,qright,flx,csl,csr,idim,i)
#endif
#if SOLVER==2        
    call solver_hllc(qleft,qright,flx,csl,csr,idim,i)
#endif

 ! Then the dust

#if NDUST>0

#if SOLVERDUST==0

    call solver_dust_Huang_Bai(qleft,qright,flx,idim,i)

#endif

#if SOLVERDUST==1
    call solver_dust_llf(qleft,qright,flx,idim,i)
#endif

#if SOLVERDUST==2
    call solver_dust_hll(qleft,qright,flx,idim,i)
#endif

#if SOLVERDUST==3
    call solver_hllc_dust(qleft,qright,flx,csl,csr,idim)
#endif
#endif

! Then te magnetic field
#if MHD==1
     
#if SOLVERB==0
    call solver_induction_llf(qleft,qright,flx,csl,csr,idim,i)
#endif

#if SOLVERB==1

    call solver_induction_Huang_Bai(qleft,qright,flx,idim,i)
#endif

#if SOLVERB==2
    call solver_induction_hll(qleft,qright,flx,csl,csr,idim,i)
#endif

#endif
end subroutine solve_wrapper













! ! Old solver

! subroutine predictor
!   use parameters
!   use commons
!   use units
!   use OMP_LIB
!   implicit none
!   integer  :: i,ix,iy,icell,iymin,iymax,idust
!   integer  :: ixx,iyy
!   integer  :: ivar,idim,ix0,iy0
!   real(dp) :: slope_lft,slope_rgt,slope_lim
!   real(dp) :: drx,dry,dpx,dpy,dux,duy,dvx,dvy,dwx,dwy,r_rho,u,v,w,p,sr0,sp0,su0,sv0,sw0,dcen,dsgn,dlim,slop
!   real(dp) :: slope_theta=1.5d0
!   real(dp), dimension(:,:,:),allocatable :: dq
!   real(dp), dimension(:,:)  ,allocatable :: qpred
!   if(static) return


!   ! Initialise to zero
!   allocate(dq(1:ncells,1:nvar,1:ndim))
!   allocate(qpred(1:ncells,1:nvar))

!   dq    = 0.0d0
!   qp    = 0.0d0
!   qm    = 0.0d0
!   ! Computes primitive variables
!   qpred = q
!   !call dust_drag_pred(qpred) ! Predictor step for the drag
!   !$OMP PARALLEL &
!   !$OMP DEFAULT(SHARED)&
!   !$OMP PRIVATE(i,ivar,ix,iy,idim,ix0,iy0,idust,slope_lft,slope_rgt,slop,dcen,dsgn,dlim,drx,dry,dpx,dpy,dux,duy,dvx,dvy,dwx,dwy,r_rho,u,v,w,p,sr0,sp0,su0,sv0,sw0)
  
!   dry = 0.0d0 
!   dpy = 0.0d0
!   duy = 0.0d0
!   dvx = 0.0d0
!   dvy = 0.0d0
!   dwx = 0.0d0
!   dwy = 0.0d0
!   v   = 0.0d0
!   w   = 0.0d0
!   !$OMP DO
!   do i=1,ncells
!     if(active_cell_predictor(i)==1) then
!       ix=ixx(i)
!       iy=iyy(i)

!       ! We compute the slope with a slope limiter
!       ! Note that in cylindrical coordinates all idim=2 derivative are divided by r (r=1 if GEOM=0)
!       if(slope_type==1) then
!       do ivar = 1,nvar
!         slope_lft    = 2.0d0*(q(i,ivar)   - q(icell(ix-1,iy),ivar))/(dx(i,1)+dx(icell(ix-1,iy),1))
!         slope_rgt    = 2.0d0*(q(icell(ix+1,iy),ivar) - q(i,ivar))/(dx(icell(ix+1,iy),1)+dx(i,1))
!         dq(i,ivar,1) = slope_lft
!         if(abs(slope_rgt)<abs(slope_lft)) dq(i,ivar,1) = slope_rgt
!         if(slope_rgt*slope_lft<0.0d0) dq(i,ivar,1) = 0.0d0
! #if NY>1
!         slope_lft    = 2.0d0*(q(i,ivar) - q(icell(ix,iy-1),ivar))/(dx(i,2)+dx(icell(ix,iy-1),2))
!         slope_rgt    = 2.0d0*(q(icell(ix,iy+1),ivar) - q(i,ivar))/(dx(icell(ix,iy+1),2)+dx(i,2))
!         dq(i,ivar,2) = slope_lft
!         if(abs(slope_rgt)<abs(slope_lft)) dq(i,ivar,2) = slope_rgt
!         if(slope_rgt*slope_lft<0.0d0) dq(i,ivar,2) = 0.0d0
! #endif        
!     end do
!       else if (slope_type==2) then ! VL
!       do ivar=1,nvar
!         slope_lft    = 2.0d0*(q(i,ivar)   - q(icell(ix-1,iy),ivar))/(dx(i,1)+dx(icell(ix-1,iy),1))
!         slope_rgt    = 2.0d0*(q(icell(ix+1,iy),ivar) - q(i,ivar))/(dx(icell(ix+1,iy),1)+dx(i,1))
!         dq(i,ivar,1) = 2.0d0*slope_lft*slope_rgt/(slope_lft+slope_rgt)
!         if(slope_rgt*slope_lft<=0.0d0) dq(i,ivar,1) = 0.0d0
! #if NY>1
!         slope_lft    = 2.0d0*(q(i,ivar) - q(icell(ix,iy-1),ivar))/(dx(i,2)+dx(icell(ix,iy-1),2))
!         slope_rgt    = 2.0d0*(q(icell(ix,iy+1),ivar) - q(i,ivar))/(dx(icell(ix,iy+1),2)+dx(i,2))
!         dq(i,ivar,2) = 2.0d0*slope_lft*slope_rgt/(slope_lft+slope_rgt)
!         if(slope_rgt*slope_lft<=0.0d0) dq(i,ivar,2) = 0.0d0
! #endif    
!     end do
!     else if (slope_type==3) then ! Moncen
!      do ivar=1,nvar
!         slope_lft  = 2.0d0*(q(i,ivar)   - q(icell(ix-1,iy),ivar))/(dx(i,1)+dx(icell(ix-1,iy),1))
!         slope_rgt  = 2.0d0*(q(icell(ix+1,iy),ivar) - q(i,ivar))/(dx(icell(ix+1,iy),1)+dx(i,1))
!         dcen = half*(slope_lft+slope_rgt)
!         dsgn = sign(1.0d0,dcen)
!         slop = min(slope_theta*abs(slope_lft),slope_theta*abs(slope_rgt))
!         dlim = slop
!         if((slope_lft*slope_rgt)<=0.0d0)dlim=0.d0
!         dq(i,ivar,1)  = dsgn*min(dlim,abs(dcen))
! #if NY>1
!         slope_lft  = 2.0d0*(q(i,ivar)   - q(icell(ix,iy-1),ivar))/(dx(i,2)+dx(icell(ix,iy-1),2))
!         slope_rgt = 2.0d0*(q(icell(ix,iy+1),ivar) - q(i,ivar))/(dx(icell(ix,iy+1),2)+dx(i,2))
!         dcen = half*(slope_lft+slope_rgt)
!         dsgn = sign(1.0d0,dcen)
!         slop = min(slope_theta*abs(slope_lft),slope_theta*abs(slope_rgt))
!         dlim = slop
!         if((slope_lft*slope_rgt)<=0.0d0)dlim=0.d0
!         dq(i,ivar,2)  = dsgn*min(dlim,abs(dcen))
! #endif    
!      end do
!     endif    
!       r_rho = q(i,irho)
!       u     = q(i,ivx)
!       dux   = dq(i,ivx,1)
!       drx   = dq(i,irho,1)
! #if NY>1
!       v     = q(i,ivy)
!       dvx   = dq(i,ivy,1)
!       duy   = dq(i,ivx,2)
!       dry   = dq(i,irho,2)
!       dvy   = dq(i,ivy,2)
!       dpy   = dq(i,iP,2)

!       w     = q(i,ivz)
!       dwx   = dq(i,ivz,1)
!       dwy   = dq(i,ivz,2)
! #endif
!       p     = q(i,iP)
!       dPx   = dq(i,iP,1)

!       sr0    = -u*drx-v*dry - (dux+dvy)*r_rho
!       sp0    = -u*dpx-v*dpy - (dux+dvy)*gamma*p
!       su0    = -u*dux-v*duy - (dpx        )/r_rho
! #if NY>1      
!       sv0    = -u*dvx-v*dvy - (dpy        )/r_rho
!       sw0    = -u*dwx-v*dwy
! #endif
!       !First, we take care of the 1D terms 
!       qpred(i,irho) = qpred(i,irho)   + half*dt*sr0
!       qpred(i,ivx)   = qpred(i,ivx)   + half*dt*su0
!       qpred(i,iP)   = qpred(i,iP)     + half*dt*sP0
!       !We now add the terms that are specific to 2D problems
! #if NY>1
!       qpred(i,ivy)   = qpred(i,ivy)   + half*dt*sv0
!       qpred(i,ivz)   = qpred(i,ivz)   + half*dt*sw0
! #endif 
!     ! Dust terms: same remark as for the gas
! #if NDUST>0
!       do idust=1,ndust
!         r_rho = q(i,irhod(idust))
!         u     = q(i,ivdx(idust))
!         dux   = dq(i,ivdx(idust),1)
!         drx   = dq(i,irhod(idust),1)
! #if NY>1
!         v     = q(i,ivdy(idust))
!         dvx   = dq(i,ivdy(idust),1)
!         dvy   = dq(i,ivdy(idust),2)
!         duy   = dq(i,ivdx(idust),2)
!         dry   = dq(i,irhod(idust),2)

!         w     = q(i,ivdz(idust))
!         dwx   = dq(i,ivdz(idust),1)
!         dwy   = dq(i,ivdz(idust),2)
! #endif
!         sr0    = -u*drx-v*dry - (dux+dvy)*r_rho
!         su0    = -u*dux-v*duy 
! #if NY>1      
!         sv0    = -u*dvx-v*dvy 
!         sw0    = -u*dwx-v*dwy
! #endif
!         qpred(i,irhod(idust)) = qpred(i,irhod(idust))    + half*dt*sr0
!         qpred(i,ivdx(idust))   = qpred(i,ivdx(idust))    + half*dt*su0
! #if NY>1
!         qpred(i,ivdy(idust))  = qpred(i,ivdy(idust))   + half*dt*sv0
!         qpred(i,ivdz(idust))  = qpred(i,ivdz(idust))   + half*dt*sw0
! #endif 

!       end do
! #endif
!   if(iso_cs==1) qpred(i,iP)   = cs(i)**2*qpred(i,irho)

!   end if
!   end do
!   !$OMP END DO

!   !$OMP BARRIER

!   !$OMP DO
!       do i=1,ncells
!         do ivar = 1,nvar
!           qm(i,ivar,1)=qpred(i,ivar) - half*dq(i,ivar,1)*dx(i,1)
!           qp(i,ivar,1)=qpred(i,ivar) + half*dq(i,ivar,1)*dx(i,1)
! #if NY>1    
!           qm(i,ivar,2)=qpred(i,ivar) - half*dq(i,ivar,2)*dx(i,2) 
!           qp(i,ivar,2)=qpred(i,ivar) + half*dq(i,ivar,2)*dx(i,2)
! #endif          
!       end do

!     end do
!   !$OMP END DO
!   !$OMP END PARALLEL
!   deallocate(dq)
!   deallocate(qpred)
! end subroutine predictor



! This is the Riemmann solver : llf + Godunov scheme
! subroutine add_delta_u
!   use parameters
!   use commons
!   use units
!   use OMP_LIB
!   use hydro_solvers
  
!   implicit none
!   integer :: i,idust,ivar,ix,iy,il,icell,idim
!   integer :: ixx,iyy

!   real(dp), dimension(:,:)  , allocatable  :: delta_U
!   real(dp), dimension(1:nvar) :: qleft,qright,flx
!   real(dp) :: csr,csl


!   if(static) return

!   !Initialisation of the flux,lambda_llf and delta U: they must be allocatable for 2D simus with h-res


!   allocate(delta_U(1:ncells,1:nvar))
!   flux       = 0.0d0
!   delta_U    = 0.0d0

!   !$OMP PARALLEL &
!   !$OMP DEFAULT(SHARED)&
!   !$OMP PRIVATE(i,idust,ivar,ix,iy,il,idim,qleft,qright,flx,csr,csl)
!   !$OMP DO
!   do i = 1, ncells
!     if(active_cell_predictor(i)==1) then
!       ix=ixx(i)
!       iy=iyy(i)
!       do idim = 1,ndim 
!         il = icell(ix+1,iy)
!         if(idim==2) then
!           il = icell(ix,iy+1)
!         endif

!         do ivar=1,nvar
!             qleft(ivar)  = qp(i,ivar,idim)
!             qright(ivar) = qm(il,ivar,idim)
!             flx(ivar)=0.0d0
!         end do

!         csl=sqrt(gamma*qleft(iP)/qleft(irho))
!         csr=sqrt(gamma*qright(iP)/qright(irho))

!         if(iso_cs==1) then
!             csl=cs(i)
!             csr=cs(il)
!         endif

!         call solver_hllc(qleft,qright,flx,csl,csr,idim)
! #if NDUST>0
!         call solver_dust_Huang_Bai(qleft,qright,flx,idim)
! #endif
!         do ivar=1,nvar
!             flux(i,ivar,idim)=flx(ivar)
!         end do
!       end do
!     end if
!   end do
!   !$OMP END DO

!   !$OMP BARRIER

!   !$OMP DO
!   do i=1,ncells
!     if(active_cell(i)==1) then
!         ix = ixx(i)
!         iy = iyy(i)
!         do ivar = 1, nvar
!             il = icell(ix-1,iy)
!             delta_U(i,ivar)=(flux(il,ivar,1)*surf(il,1)-flux(i,ivar,1)*surf(i,1))/vol(i)*dt
!             if(ndim==2) then
!               il = icell(ix,iy-1)
!               delta_U(i,ivar)=delta_U(i,ivar)+(flux(il,ivar,2)*surf(il,2)-flux(i,ivar,2)*surf(i,2))/vol(i)*dt
!             endif
!         end do
!       endif
!     end do
!   !end do
!   !$OMP END DO
!   !$OMP END PARALLEL
!   !Update state vector 
!   u_prim=u_prim+delta_U


!   deallocate(delta_U)
! end subroutine add_delta_u




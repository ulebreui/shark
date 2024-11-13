
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
  integer  :: irho_spe,ivx_spe,ivy_spe,ivz_spe,ipscal


  real(dp), dimension(:,:,:),allocatable :: dq
  if(static) return


  ! Initialise to zero
  allocate(dq(1:ncells,1:nvar,1:ndim))

  dq    = 0.0d0
  qp    = 0.0d0
  qm    = 0.0d0
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

  do iy=2,ny_max-1
    do ix= 2,nx_max-1
        i=icell(ix,iy)
      radius_polar=1.0d0     
#if GEOM==2
      radius_polar=radii_c(i)
#endif
#if GEOM==4
      radius_polar=1.0d0
#endif  
    if(slope_type>0) then
    do ivar = 1, nvar
            dq(i,ivar,1) = slope_limit(2.0d0*(q(ix,iy,ivar) - q(ix-1,iy,ivar))/(dx(ix,iy,1)+dx(ix-1,iy,1)),2.0d0*(q(ix+1,iy,ivar) - q(ix,iy,ivar))/(dx(ix+1,iy,1)+dx(ix,iy,1)))
#if NY>1
            dq(i,ivar,2) = slope_limit(2.0d0*(q(ix,iy,ivar) - q(ix,iy-1,ivar))/(dx(ix,iy,2)+dx(ix,iy-1,2)),2.0d0*(q(ix,iy+1,ivar) - q(ix,iy,ivar))/(dx(ix,iy+1,2)+dx(ix,iy,2)))/radius_polar
#endif   
    end do
    endif

      r_rho = q(ix,iy,irho)
      u     = q(ix,iy,ivx)
      dux   = dq(i,ivx,1)
      drx   = dq(i,irho,1)
      w     = q(ix,iy,ivz)
      dwx   = dq(i,ivz,1)
      v     = q(ix,iy,ivy)
      dvx   = dq(i,ivy,1)
      p     = q(ix,iy,iP)
      dPx   = dq(i,iP,1)
      if(iso_cs==1) then
       P     = q(ix,iy,irho)*cs(i)**2
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


!Cartesian geometry
#if GE0M==0
      sr0    = -u*drx-v*dry - (dux+dvy)*r_rho
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


    !direction x
    dx_loc=dx(ix,iy,1)
    qm(ix,iy,irho,1)   = r_rho     + half*dt*sr0  + half*drx   *dx_loc
    qm(ix,iy,ivx,1)    = u         + half*dt*su0  + half*dux   *dx_loc
    qm(ix,iy,iP,1)     = max(p     + half*dt*sP0  + half*dpx   *dx_loc,smallP)
    qm(ix,iy,ivy,1)    = v         + half*dt*sv0  + half*dvx   *dx_loc
    qm(ix,iy,ivz,1)    = w         + half*dt*sw0  + half*dwx   *dx_loc 
    qp(ix,iy,irho,1)   =r_rho      + half*dt*sr0  - half*drx   *dx_loc
    qp(ix,iy,ivx,1)    = u         + half*dt*su0  - half*dux   *dx_loc
    qp(ix,iy,iP,1)     = max(p     + half*dt*sP0  - half*dpx   *dx_loc,smallP)
    qp(ix,iy,ivy,1)    = v         + half*dt*sv0  - half*dvx   *dx_loc
    qp(ix,iy,ivz,1)    = w         + half*dt*sw0  - half*dwx   *dx_loc

    ! direction y
#if NY>1
    dx_loc=radius_polar*dx(ix,iy,2)
    qm(ix,iy,irho,2)   = r_rho     + half*dt*sr0  + half*dry   *dx_loc
    qm(ix,iy,ivx,2)    = u         + half*dt*su0  + half*duy   *dx_loc
    qm(ix,iy,iP,2)     = max(p     + half*dt*sP0  + half*dpy   *dx_loc,smallP)
    qm(ix,iy,ivy,2)    = v         + half*dt*sv0  + half*dvy   *dx_loc
    qm(ix,iy,ivz,2)    = w         + half*dt*sw0  + half*dwy   *dx_loc 
    qp(ix,iy,irho,2)   = r_rho     + half*dt*sr0  - half*dry   *dx_loc
    qp(ix,iy,ivx,2)    = u         + half*dt*su0  - half*duy   *dx_loc
    qp(ix,iy,iP,2)     = max(p     + half*dt*sP0  - half*dpy   *dx_loc,smallP)
    qp(ix,iy,ivy,2)    = v         + half*dt*sv0  - half*dvy   *dx_loc
    qp(ix,iy,ivz,2)    = w         + half*dt*sw0  - half*dwy   *dx_loc
 
#endif

    ! Dust terms: same remark as for the gas
#if NDUST>0


      do idust=1,ndust
        irho_spe = irhod(idust)
        ivx_spe  = ivdx(idust)
        ivy_spe  = ivdy(idust) 
        ivz_spe  = ivdz(idust)

        r_rho = q(ix,iy,irho_spe)
        u     = q(ix,iy,ivx_spe)
        dux   = dq(i,ivx_spe,1)
        drx   = dq(i,irho_spe,1)
        v     = q(ix,iy,ivy_spe)
        dvx   = dq(i,ivy_spe,1)
        w     = q(ix,iy,ivz_spe)
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

    !Direction x
    dx_loc=dx(ix,iy,1)
    qm(ix,iy,irho_spe,1)   = r_rho     + half*dt*sr0  + half*drx   *dx_loc
    qm(ix,iy,ivx_spe,1)    = u         + half*dt*su0  + half*dux   *dx_loc
    qm(ix,iy,ivy_spe,1)    = v         + half*dt*sv0  + half*dvx   *dx_loc
    qm(ix,iy,ivz_spe,1)    = w         + half*dt*sw0  + half*dwx   *dx_loc 
    qp(ix,iy,irho_spe,1)   = r_rho     + half*dt*sr0  - half*drx   *dx_loc
    qp(ix,iy,ivx_spe,1)    = u         + half*dt*su0  - half*dux   *dx_loc
    qp(ix,iy,ivy_spe,1)    = v         + half*dt*sv0  - half*dvx   *dx_loc
    qp(ix,iy,ivz_spe,1)    = w         + half*dt*sw0  - half*dwx   *dx_loc
    
    !direction y
#if NY>1
    dx_loc=radius_polar*dx(ix,iy,2)
    qm(ix,iy,irho_spe,2)   = r_rho     + half*dt*sr0  + half*dry   *dx_loc
    qm(ix,iy,ivx_spe,2)    = u         + half*dt*su0  + half*duy   *dx_loc
    qm(ix,iy,ivy_spe,2)    = v         + half*dt*sv0  + half*dvy   *dx_loc
    qm(ix,iy,ivz_spe,2)    = w         + half*dt*sw0  + half*dwy   *dx_loc 
    qp(ix,iy,irho_spe,2)   = r_rho     + half*dt*sr0  - half*dry   *dx_loc
    qp(ix,iy,ivx_spe,2)    = u         + half*dt*su0  - half*duy   *dx_loc
    qp(ix,iy,ivy_spe,2)    = v         + half*dt*sv0  - half*dvy   *dx_loc
    qp(ix,iy,ivz_spe,2)    = w         + half*dt*sw0  - half*dwy   *dx_loc

#endif

#if NDUSTPSCAL>0
    do ipscal= 1, ndustpscal
        r_rho = q(ix,iy,idust_pscal(idust,ipscal))

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

        dx_loc=dx(ix,iy,1)
        qm(ix,iy,idust_pscal(idust,ipscal),1)   = r_rho + half*dt*sr0  + half*drx   * dx_loc
        qp(ix,iy,idust_pscal(idust,ipscal),1)   = r_rho + half*dt*sr0  - half*drx   * dx_loc
        
        ! Direction y
#if NY>1
       dx_loc=radius_polar*dx(ix,iy,2)
       qm(ix,iy,idust_pscal(idust,ipscal),2)   = r_rho + half*dt*sr0  + half*dry   * dx_loc
       qp(ix,iy,idust_pscal(idust,ipscal),2)   = r_rho + half*dt*sr0  - half*dry   * dx_loc
#endif
    end do
#endif
      end do !dust loop
#endif

! We recompute the thermal pressure if we don't solve for the NRJ equation
if(iso_cs==1) then
    do idim =1,ndim
        qp(ix,iy,iP,idim)   = cs(i)**2*qp(ix,iy,irho,idim)
        qm(ix,iy,iP,idim)   = cs(i)**2*qm(ix,iy,irho,idim)
    end do
end if
if(non_standard_eos==1) then
    do idim=1,ndim
        qp(ix,iy,iP,idim)   = cs_eos(barotrop(qp(ix,iy,irho,idim)))**2*qp(ix,iy,irho,idim)
        qm(ix,iy,iP,idim)   = cs_eos(barotrop(qm(ix,iy,irho,idim)))**2*qm(ix,iy,irho,idim)
   end do
end if
end do 
end do 

  deallocate(dq)

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


  do iy=2,ny_max-1
    do ix= 2,nx_max-1
      i=icell(ix,iy)

      !First direction x
      il = icell(ix-1,iy)
      do ivar=1,nvar
            qleft(ivar)  = qm(ix-1,iy,ivar,1)
            qright(ivar) = qp(ix,iy,ivar,1)
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
            qleft(ivar)  = qm(ix,iy-1,ivar,2)
            qright(ivar) = qp(ix,iy,ivar, 2)
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

    end do
  end do

  do iy=first_active_y,last_active_y
    do ix= first_active,last_active
        i=icell(ix,iy)
        ! 2D case : we do both fluxes at the same time to be more efficient
        il  = icell(ix+1,iy)
        ily = icell(ix,iy+1)
        do ivar = 1, nvar
              u_prim(ix,iy,ivar)=u_prim(ix,iy,ivar) + (flux(i,ivar,1)*surf(ix,iy,1)-flux(il,ivar,1) *surf(ix+1,iy,1))  /vol(ix,iy)*dt&
              &                             + (flux(i,ivar,2)*surf(ix,iy,2)-flux(ily,ivar,2)*surf(ix,iy+1,2))/vol(ix,iy)*dt
        end do

      end do
    end do
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

end subroutine solve_wrapper


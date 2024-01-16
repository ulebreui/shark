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
  integer  :: i,ix,iy,icell,iymin,iymax,idust
  integer  :: ixx,iyy
  integer  :: ivar,idim,ix0,iy0
  real(dp) :: slope_lft,slope_rgt,slope_lim,barotrop,cs_eos,ddxp,ddxm
  real(dp) :: drx,dry,dpx,dpy,dux,duy,dvx,dvy,dwx,dwy,r_rho,u,v,w,p,sr0,sp0,su0,sv0,sw0,dcen,dsgn,dlim,slop,radius_polar


  real(dp) :: dBx_x,dBy_x,dBz_x,Bx,By,Bz,sBx,sBy,sBz


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
  !$OMP PRIVATE(i,ivar,ix,iy,idim,ix0,iy0,idust,slope_lft,slope_rgt,ddxp,ddxm,drx,dry,dpx,dpy,dux,duy,dvx,dvy,dwx,dwy,r_rho,u,v,w,p,sr0,sp0,su0,sv0,sw0,radius_polar,dBx_x,dBy_x,dBz_x,Bx,By,Bz,sBx,sBy,sBz)
  
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

  !$OMP DO
  do i=1,ncells
    if(active_cell_predictor(i)==1) then
      ix=ixx(i)
      iy=iyy(i)
#if GEOM==2
      radius_polar=radii_c(i)
#endif
      ! We compute the slope with a slope limiter
      ! Note that in cylindrical coordinates all idim=2 derivative are divided by r (r=1 if GEOM=0)
      if(slope_type>0) then
        do ivar = 1,nvar
#if GEOM<2
            dq(i,ivar,1) = slope_limit(2.0d0*(q(i,ivar) - q(icell(ix-1,iy),ivar))/(dx(i,1)+dx(icell(ix-1,iy),1)),2.0d0*(q(icell(ix+1,iy),ivar) - q(i,ivar))/(dx(icell(ix+1,iy),1)+dx(i,1)))
#if NY>1
            dq(i,ivar,2) = slope_limit(2.0d0*(q(i,ivar) - q(icell(ix,iy-1),ivar))/(dx(i,2)+dx(icell(ix,iy-1),2)),2.0d0*(q(icell(ix,iy+1),ivar) - q(i,ivar))/(dx(icell(ix,iy+1),2)+dx(i,2)))
#endif  
#endif

! ! Special collapse modif for spherical grid
! #if GEOM==1
! #if NY==1
!             dq(i,ivar,1) = slope_limit(2.0d0*(q(i,ivar) - q(icell(ix-1,iy),ivar))/dx_l(i),2.0d0*(q(icell(ix+1,iy),ivar) - q(i,ivar))/dx_r(i))
!             !if(i<=first_active .or. i>=last_active) dq(i,ivar,1)= 0.0d0
! #endif  
! #endif  

#if GEOM==2
            dq(i,ivar,1) = slope_limit(2.0d0*(q(i,ivar) - q(icell(ix-1,iy),ivar))/(dx(i,1)+dx(icell(ix-1,iy),1)),2.0d0*(q(icell(ix+1,iy),ivar) - q(i,ivar))/(dx(icell(ix+1,iy),1)+dx(i,1)))
#if NY>1
            dq(i,ivar,2) = slope_limit(2.0d0*(q(i,ivar) - q(icell(ix,iy-1),ivar))/(dx(i,2)+dx(icell(ix,iy-1),2)),2.0d0*(q(icell(ix,iy+1),ivar) - q(i,ivar))/(dx(icell(ix,iy+1),2)+dx(i,2)))/radius_polar
#endif  
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
#if NY>1
      duy   = dq(i,ivx,2)
      dry   = dq(i,irho,2)
      dvy   = dq(i,ivy,2)
      dpy   = dq(i,iP,2)
      dwy   = dq(i,ivz,2)
      
#endif
      

#if MHD==1
      Bx    = q(i,iBx)
      By    = q(i,iBy)
      Bz    = q(i,iBz)
      dBx_x    = dq(i,iBx,1)
      dBy_x    = dq(i,iBy,1)
      dBz_x    = dq(i,iBz,1)

#endif

      sr0    = -u*drx-v*dry - (dux+dvy)*r_rho

#if GEOM==1
      !Spherical geometry source term
      sr0    = sr0 - r_rho*u/radii_c(i)
#endif       

      sp0    = -u*dpx-v*dpy - (dux+dvy)*gamma*p
      su0    = -u*dux-v*duy - (dpx        )/r_rho
      sv0    = -u*dvx-v*dvy - (dpy        )/r_rho
      sw0    = -u*dwx-v*dwy
#if GEOM==2
      !Polar geometry source terms
      sr0    = sr0 - r_rho*u  / radius_polar
      su0    = su0 + (v**2.)  / radius_polar
      sv0    = sv0 - u*v      / radius_polar
      sp0    = sp0 -gamma*p*u / radius_polar 
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
    
    if(force_kick) then
        su0    = su0 + force(i,1)
        sv0    = sv0 + force(i,2)
        sw0    = sw0 + force(i,3)
    endif


      do idim=1,ndim

#if GEOM==0
        ddxm=dx(i,idim)
        ddxp=dx(i,idim)
#endif
#if GEOM==1 
!       This is a correction for log grids dx_l =!= dx_r
#if NY==1
        ! ddxm=dx_l_cell(i)
        ! ddxp=dx_r_cell(i)
        ddxm=dx(i,idim)
        ddxp=dx(i,idim)
#endif
#endif
#if GEOM==2
        if(idim==1) then
            ddxm=dx(i,idim)
            ddxp=dx(i,idim)
        else
            ddxm=radius_polar*dx(i,idim)
            ddxp=radius_polar*dx(i,idim)
        endif 
#endif
        qm(i,irho,idim)   = max(q(i,irho)    + half*dt*sr0  + half*dq(i,irho,idim)   *ddxm,smallr)
        qm(i,ivx,idim)    = q(i,ivx)         + half*dt*su0  + half*dq(i,ivx,idim)    *ddxm
        qm(i,iP,idim)     = max(q(i,iP)      + half*dt*sP0  + half*dq(i,iP,idim)     *ddxm,smallP)
        qm(i,ivy,idim)    = q(i,ivy)         + half*dt*sv0  + half*dq(i,ivy,idim)    *ddxm
        qm(i,ivz,idim)    = q(i,ivz)         + half*dt*sw0  + half* dq(i,ivz,idim)   *ddxm 
        qp(i,irho,idim)   = max(q(i,irho)    + half*dt*sr0  - half*dq(i,irho,idim)   *ddxp,smallr)
        qp(i,ivx,idim)    = q(i,ivx)         + half*dt*su0  - half*dq(i,ivx,idim)    *ddxp
        qp(i,iP,idim)     = max(q(i,iP)      + half*dt*sP0  - half*dq(i,iP,idim)     *ddxp,smallP)
        qp(i,ivy,idim)    = q(i,ivy)         + half*dt*sv0  - half*dq(i,ivy,idim)    *ddxp
        qp(i,ivz,idim)    = q(i,ivz)         + half*dt*sw0  - half* dq(i,ivz,idim)   *ddxp

    end do

    ! Dust terms: same remark as for the gas
#if NDUST>0
      do idust=1,ndust
        r_rho = q(i,irhod(idust))
        u     = q(i,ivdx(idust))
        dux   = dq(i,ivdx(idust),1)
        drx   = dq(i,irhod(idust),1)
        v     = q(i,ivdy(idust))
        dvx   = dq(i,ivdy(idust),1)
        w     = q(i,ivdz(idust))
        dwx   = dq(i,ivdz(idust),1)
#if NY>1        
        dvy   = dq(i,ivdy(idust),2)
        duy   = dq(i,ivdx(idust),2)
        dry   = dq(i,irhod(idust),2)
        dwy   = dq(i,ivdz(idust),2)
#endif
        sr0    = -u*drx-v*dry - (dux+dvy)*r_rho
#if GEOM==1
        !Spherical geometry source term
        sr0    = sr0-r_rho*u/radii_c(i)
#endif         
#if MHD==0
        su0    = -u*dux-v*duy       
        sv0    = -u*dvx-v*dvy 
        sw0    = -u*dwx-v*dwy
#endif
#if GEOM==2
        !Polar geometry source terms
        sr0    = sr0 - r_rho*u/radius_polar
        su0    = su0 + v**2.  /radius_polar
        sv0    = sv0 - u*v    /radius_polar
#endif
#if MHD==1

    if (idust==i_coupled_species) then
        sBy   = -u*dBy_x + Bx*dvx - By*dux
        sBz   = -u*dBz_x + Bx*dwx - Bz*dux
    endif
        su0   = -u*dux-v*duy - By/r_rho*dBy_x - Bz/r_rho*dBz_x - Bx*dBx_x/r_rho + 2.0d0*Bx*dBx_x
        sv0   = -u*dvx-v*dvy + Bx*dBy_x/r_rho
        sw0   = -u*dwx-v*dwy + Bx*dBz_x/r_rho

#endif
#if GRAVITY==1
#if NY==1
        !Gravity source term in 1D
        su0=su0+(-Mc(i)/(radii_c(i)**2.+(l_soft/unit_l)**2.))
#endif
#endif   

    if(force_kick) then
        su0    = su0 + force_dust(i,idust,1)
        sv0    = sv0 + force_dust(i,idust,2)
        sw0    = sw0 + force_dust(i,idust,3)
    endif

      do idim=1,ndim
#if GEOM==0
        ddxm=dx(i,idim)
        ddxp=dx(i,idim)
#endif
#if GEOM==1 
!       This is a correction for log grids dx_l =!= dx_r
#if NY==1
        ! ddxm=dx_l_cell(i)
        ! ddxp=dx_r_cell(i)
        ddxm=dx(i,idim)
        ddxp=dx(i,idim)
#endif
#endif
        qm(i,irhod(idust),idim)   = q(i,irhod(idust))    + half*dt*sr0 + half*dq(i,irhod(idust),idim) * ddxm 
        qm(i,ivdx(idust),idim)    = q(i,ivdx(idust))     + half*dt*su0 + half*dq(i,ivdx(idust),idim)  * ddxm 
        qm(i,ivdy(idust),idim)    = q(i,ivdy(idust))     + half*dt*sv0 + half*dq(i,ivdy(idust),idim)  * ddxm 
        qm(i,ivdz(idust),idim)    = q(i,ivdz(idust))     + half*dt*sw0 + half*dq(i,ivdz(idust),idim)  * ddxm 
        
        qp(i,irhod(idust),idim)   = q(i,irhod(idust))    + half*dt*sr0 - half*dq(i,irhod(idust),idim)  *  ddxp
        qp(i,ivdx(idust),idim)    = q(i,ivdx(idust))     + half*dt*su0 - half*dq(i,ivdx(idust),idim)   *  ddxp
        qp(i,ivdy(idust),idim)    = q(i,ivdy(idust))     + half*dt*sv0 - half*dq(i,ivdy(idust),idim)   *  ddxp
        qp(i,ivdz(idust),idim)    = q(i,ivdz(idust))     + half*dt*sw0 - half*dq(i,ivdz(idust),idim)   *  ddxp

    end do


      end do
#endif

#if MHD==1
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
  use hydro_solvers

  implicit none
  integer :: i,idust,ivar,ix,iy,il,icell,idim
  integer :: ixx,iyy

  real(dp), dimension(:,:)  , allocatable  :: delta_U
  real(dp), dimension(:,:,:), allocatable  :: flux
  real(dp), dimension(1:nvar) :: qleft,qright,flx
  real(dp) :: csr,csl,barotrop,cs_eos


  if(static) return

  ! Initialisation of the flux,lambda_llf and delta U: they must be allocatable for 2D simus with h-res


  allocate(delta_U(1:ncells,1:nvar))
  allocate(flux(1:ncells,1:nvar,1:ndim))
  flux       = 0.0d0
  delta_U    = 0.0d0

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,ivar,ix,iy,il,idim,qleft,qright,flx,csr,csl)
  !$OMP DO
  do i = 1, ncells
    if(active_cell_predictor(i)==1) then
      ix=ixx(i)
      iy=iyy(i)
      do idim = 1,ndim 
        il = icell(ix-1,iy)
        if(idim==2) then
          il = icell(ix,iy-1)
        endif

        do ivar=1,nvar
            qleft(ivar)  = qm(il,ivar,idim)
            qright(ivar) = qp(i,ivar,idim)
            flx(ivar)    = 0.0d0
        end do

        csl=sqrt(gamma*qleft(iP)/qleft(irho))
        csr=sqrt(gamma*qright(iP)/qright(irho))
        if(iso_cs==1) then
            csl=cs(il)
            csr=cs(i)
        endif
        if(non_standard_eos==1) then
            csl=cs_eos(barotrop(qleft(irho)))
            csr=cs_eos(barotrop(qright(irho)))
        endif

#if SOLVER==0        
        call solver_llf(qleft,qright,flx,csl,csr,idim)
#endif
#if SOLVER==1       
        call solver_hll(qleft,qright,flx,csl,csr,idim)
#endif
#if SOLVER==2        
        call solver_hllc(qleft,qright,flx,csl,csr,idim)
#endif



#if NDUST>0

#if SOLVERDUST==0

    call solver_dust_Huang_Bai(qleft,qright,flx,idim)

#endif

#if SOLVERDUST==1
        call solver_dust_llf(qleft,qright,flx,idim)
#endif

#if SOLVERDUST==2
        call solver_dust_hll(qleft,qright,flx,idim)
#endif

#endif


#if MHD==1
     
#if SOLVERB==0

        call solver_induction_llf(qleft,qright,flx,csl,csr,idim)
#endif

#if SOLVERB==1

        call solver_induction_Huang_Bai(qleft,qright,flx,idim)
#endif

#if SOLVERB==2

        call solver_induction_hll(qleft,qright,flx,csl,csr,idim)
#endif

#endif


        do ivar=1,nvar
            flux(i,ivar,idim)=flx(ivar) 
        end do
      end do
    end if
  end do
  !$OMP END DO
 !stop
  !$OMP BARRIER

  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
        ix = ixx(i)
        iy = iyy(i)
        do ivar = 1, nvar
            delta_U(i,ivar)=(flux(i,ivar,1)*surf(i,1)-flux(icell(ix+1,iy),ivar,1)*surf(icell(ix+1,iy),1))/vol(i)*dt
            if(ndim==2) then
              delta_U(i,ivar)=delta_U(i,ivar)+(flux(i,ivar,2)*surf(i,2)-flux(icell(ix,iy+1),ivar,2)*surf(icell(ix,iy+1),2))/vol(i)*dt
            endif
        end do
      endif
    end do
  !end do
  !$OMP END DO
  !$OMP END PARALLEL
  !Update state vector 
  u_prim=u_prim+delta_U



  deallocate(delta_U)
  deallocate(flux)
end subroutine add_delta_u



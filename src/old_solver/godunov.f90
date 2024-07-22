! This routine performs the predictor operation ! 

! Reference for cylindrical geom :
! The Athena Astrophysical MHD Code in Cylindrical Geometry
! A Proper Discretization of Hydrodynamic Equations in Cylindrical Coordinates for Astrophysical Simulations
subroutine predictor
  use parameters
  use commons
  use units
  use OMP_LIB
  implicit none
  integer  :: i,ix,iy,icell,iymin,iymax,idust
  integer  :: ixx,iyy
  integer  :: ivar,idim,ix0,iy0
  real(dp) :: slope_lft,slope_rgt,slope_lim
  real(dp) :: drx,dry,dpx,dpy,dux,duy,dvx,dvy,dwx,dwy,r_rho,u,v,w,p,sr0,sp0,su0,sv0,sw0,dcen,dsgn,dlim,slop
  real(dp) :: slope_theta=1.5d0
  real(dp), dimension(:,:,:),allocatable :: dq
  real(dp), dimension(:,:)  ,allocatable :: qpred
  if(static) return


  ! Initialise to zero
  allocate(dq(1:ncells,1:nvar,1:ndim))
  allocate(qpred(1:ncells,1:nvar))

  dq    = 0.0d0
  qp    = 0.0d0
  qm    = 0.0d0
  ! Computes primitive variables
  qpred = q
  !call dust_drag_pred(qpred) ! Predictor step for the drag
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,ivar,ix,iy,idim,ix0,iy0,idust,slope_lft,slope_rgt,slop,dcen,dsgn,dlim,drx,dry,dpx,dpy,dux,duy,dvx,dvy,dwx,dwy,r_rho,u,v,w,p,sr0,sp0,su0,sv0,sw0)
  
  dry = 0.0d0 
  dpy = 0.0d0
  duy = 0.0d0
  dvx = 0.0d0
  dvy = 0.0d0
  dwx = 0.0d0
  dwy = 0.0d0
  v   = 0.0d0
  w   = 0.0d0
  !$OMP DO
  do i=1,ncells
    if(active_cell_predictor(i)==1) then
      ix=ixx(i)
      iy=iyy(i)

      ! We compute the slope with a slope limiter
      ! Note that in cylindrical coordinates all idim=2 derivative are divided by r (r=1 if GEOM=0)
      if(slope_type==1) then
      do ivar = 1,nvar
        slope_lft    = 2.0d0*(q(i,ivar)   - q(icell(ix-1,iy),ivar))/(dx(i,1)+dx(icell(ix-1,iy),1))
        slope_rgt    = 2.0d0*(q(icell(ix+1,iy),ivar) - q(i,ivar))/(dx(icell(ix+1,iy),1)+dx(i,1))
        dq(i,ivar,1) = slope_lft
        if(abs(slope_rgt)<abs(slope_lft)) dq(i,ivar,1) = slope_rgt
        if(slope_rgt*slope_lft<0.0d0) dq(i,ivar,1) = 0.0d0
#if NY>1
        slope_lft    = 2.0d0*(q(i,ivar) - q(icell(ix,iy-1),ivar))/(dx(i,2)+dx(icell(ix,iy-1),2))
        slope_rgt    = 2.0d0*(q(icell(ix,iy+1),ivar) - q(i,ivar))/(dx(icell(ix,iy+1),2)+dx(i,2))
        dq(i,ivar,2) = slope_lft
        if(abs(slope_rgt)<abs(slope_lft)) dq(i,ivar,2) = slope_rgt
        if(slope_rgt*slope_lft<0.0d0) dq(i,ivar,2) = 0.0d0
#endif        
    end do
      else if (slope_type==2) then ! VL
      do ivar=1,nvar
        slope_lft    = 2.0d0*(q(i,ivar)   - q(icell(ix-1,iy),ivar))/(dx(i,1)+dx(icell(ix-1,iy),1))
        slope_rgt    = 2.0d0*(q(icell(ix+1,iy),ivar) - q(i,ivar))/(dx(icell(ix+1,iy),1)+dx(i,1))
        dq(i,ivar,1) = 2.0d0*slope_lft*slope_rgt/(slope_lft+slope_rgt)
        if(slope_rgt*slope_lft<=0.0d0) dq(i,ivar,1) = 0.0d0
#if NY>1
        slope_lft    = 2.0d0*(q(i,ivar) - q(icell(ix,iy-1),ivar))/(dx(i,2)+dx(icell(ix,iy-1),2))
        slope_rgt    = 2.0d0*(q(icell(ix,iy+1),ivar) - q(i,ivar))/(dx(icell(ix,iy+1),2)+dx(i,2))
        dq(i,ivar,2) = 2.0d0*slope_lft*slope_rgt/(slope_lft+slope_rgt)
        if(slope_rgt*slope_lft<=0.0d0) dq(i,ivar,2) = 0.0d0
#endif    
    end do
    else if (slope_type==3) then ! Moncen
     do ivar=1,nvar
        slope_lft  = 2.0d0*(q(i,ivar)   - q(icell(ix-1,iy),ivar))/(dx(i,1)+dx(icell(ix-1,iy),1))
        slope_rgt  = 2.0d0*(q(icell(ix+1,iy),ivar) - q(i,ivar))/(dx(icell(ix+1,iy),1)+dx(i,1))
        dcen = half*(slope_lft+slope_rgt)
        dsgn = sign(1.0d0,dcen)
        slop = min(slope_theta*abs(slope_lft),slope_theta*abs(slope_rgt))
        dlim = slop
        if((slope_lft*slope_rgt)<=0.0d0)dlim=0.d0
        dq(i,ivar,1)  = dsgn*min(dlim,abs(dcen))
#if NY>1
        slope_lft  = 2.0d0*(q(i,ivar)   - q(icell(ix,iy-1),ivar))/(dx(i,2)+dx(icell(ix,iy-1),2))
        slope_rgt = 2.0d0*(q(icell(ix,iy+1),ivar) - q(i,ivar))/(dx(icell(ix,iy+1),2)+dx(i,2))
        dcen = half*(slope_lft+slope_rgt)
        dsgn = sign(1.0d0,dcen)
        slop = min(slope_theta*abs(slope_lft),slope_theta*abs(slope_rgt))
        dlim = slop
        if((slope_lft*slope_rgt)<=0.0d0)dlim=0.d0
        dq(i,ivar,2)  = dsgn*min(dlim,abs(dcen))
#endif    
     end do
    endif    
      r_rho = q(i,irho)
      u     = q(i,ivx)
      dux   = dq(i,ivx,1)
      drx   = dq(i,irho,1)
#if NY>1
      v     = q(i,ivy)
      dvx   = dq(i,ivy,1)
      duy   = dq(i,ivx,2)
      dry   = dq(i,irho,2)
      dvy   = dq(i,ivy,2)
      dpy   = dq(i,iP,2)

      w     = q(i,ivz)
      dwx   = dq(i,ivz,1)
      dwy   = dq(i,ivz,2)
#endif
      p     = q(i,iP)
      dPx   = dq(i,iP,1)

      sr0    = -u*drx-v*dry - (dux+dvy)*r_rho
      sp0    = -u*dpx-v*dpy - (dux+dvy)*gamma*p
      su0    = -u*dux-v*duy - (dpx        )/r_rho
#if NY>1      
      sv0    = -u*dvx-v*dvy - (dpy        )/r_rho
      sw0    = -u*dwx-v*dwy
#endif
      !First, we take care of the 1D terms 
      qpred(i,irho) = qpred(i,irho)   + half*dt*sr0
      qpred(i,ivx)   = qpred(i,ivx)   + half*dt*su0
      qpred(i,iP)   = qpred(i,iP)     + half*dt*sP0
      !We now add the terms that are specific to 2D problems
#if NY>1
      qpred(i,ivy)   = qpred(i,ivy)   + half*dt*sv0
      qpred(i,ivz)   = qpred(i,ivz)   + half*dt*sw0
#endif 
    ! Dust terms: same remark as for the gas
#if NDUST>0
      do idust=1,ndust
        r_rho = q(i,irhod(idust))
        u     = q(i,ivdx(idust))
        dux   = dq(i,ivdx(idust),1)
        drx   = dq(i,irhod(idust),1)
#if NY>1
        v     = q(i,ivdy(idust))
        dvx   = dq(i,ivdy(idust),1)
        dvy   = dq(i,ivdy(idust),2)
        duy   = dq(i,ivdx(idust),2)
        dry   = dq(i,irhod(idust),2)

        w     = q(i,ivdz(idust))
        dwx   = dq(i,ivdz(idust),1)
        dwy   = dq(i,ivdz(idust),2)
#endif
        sr0    = -u*drx-v*dry - (dux+dvy)*r_rho
        su0    = -u*dux-v*duy 
#if NY>1      
        sv0    = -u*dvx-v*dvy 
        sw0    = -u*dwx-v*dwy
#endif
        qpred(i,irhod(idust)) = qpred(i,irhod(idust))    + half*dt*sr0
        qpred(i,ivdx(idust))   = qpred(i,ivdx(idust))    + half*dt*su0
#if NY>1
        qpred(i,ivdy(idust))  = qpred(i,ivdy(idust))   + half*dt*sv0
        qpred(i,ivdz(idust))  = qpred(i,ivdz(idust))   + half*dt*sw0
#endif 

      end do
#endif
  if(iso_cs==1) qpred(i,iP)   = cs(i)**2*qpred(i,irho)

  end if
  end do
  !$OMP END DO

  !$OMP BARRIER

  !$OMP DO
      do i=1,ncells
        do ivar = 1,nvar
          qm(i,ivar,1)=qpred(i,ivar) - half*dq(i,ivar,1)*dx(i,1)
          qp(i,ivar,1)=qpred(i,ivar) + half*dq(i,ivar,1)*dx(i,1)
#if NY>1    
          qm(i,ivar,2)=qpred(i,ivar) - half*dq(i,ivar,2)*dx(i,2) 
          qp(i,ivar,2)=qpred(i,ivar) + half*dq(i,ivar,2)*dx(i,2)
#endif          
      end do

    end do
  !$OMP END DO
  !$OMP END PARALLEL
  deallocate(dq)
  deallocate(qpred)
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
  real(dp), dimension(1:nvar) :: qleft,qright,flx
  real(dp) :: csr,csl


  if(static) return

  !Initialisation of the flux,lambda_llf and delta U: they must be allocatable for 2D simus with h-res


  allocate(delta_U(1:ncells,1:nvar))
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
        il = icell(ix+1,iy)
        if(idim==2) then
          il = icell(ix,iy+1)
        endif

        do ivar=1,nvar
            qleft(ivar)  = qp(i,ivar,idim)
            qright(ivar) = qm(il,ivar,idim)
            flx(ivar)=0.0d0
        end do

        csl=sqrt(gamma*qleft(iP)/qleft(irho))
        csr=sqrt(gamma*qright(iP)/qright(irho))

        if(iso_cs==1) then
            csl=cs(i)
            csr=cs(il)
        endif

        call solver_hllc(qleft,qright,flx,csl,csr,idim)
#if NDUST>0
        call solver_dust(qleft,qright,flx,idim)
#endif
        do ivar=1,nvar
            flux(i,ivar,idim)=flx(ivar)
        end do
      end do
    end if
  end do
  !$OMP END DO

  !$OMP BARRIER

  !$OMP DO
  do i=1,ncells
    if(active_cell(i)==1) then
        ix = ixx(i)
        iy = iyy(i)
        do ivar = 1, nvar
            il = icell(ix-1,iy)
            delta_U(i,ivar)=(flux(il,ivar,1)*surf(il,1)-flux(i,ivar,1)*surf(i,1))/vol(i)*dt
            if(ndim==2) then
              il = icell(ix,iy-1)
              delta_U(i,ivar)=delta_U(i,ivar)+(flux(il,ivar,2)*surf(il,2)-flux(i,ivar,2)*surf(i,2))/vol(i)*dt
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
end subroutine add_delta_u



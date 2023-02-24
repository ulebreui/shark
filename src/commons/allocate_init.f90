!Allocation of all the quantities
subroutine allocate_init
  use parameters
  use commons
  use units
  implicit none

  integer :: idust

  call get_active_cells

  ! Grid related quantities
  allocate(vol(1:ncells))
  allocate(dx(1:ncells,1:ndim))
  allocate(surf(1:ncells,1:ndim))
  allocate(position(1:ncells,1:ndim))

  vol  = 0.d0
  dx   = 0.d0
  Surf = 0.d0
  position = 0.0d0
#if GEOM==1
  allocate(polar_radii(1:ncells))
  polar_radii=0.0d0 
  allocate(theta(1:ncells))
  theta=0.0d0 
#endif
  ! Variable related quantities
  allocate(u_prim(1:ncells,1:nvar))
  allocate(q(1:ncells,1:nvar))
  allocate(qm(1:ncells,1:nvar,1:ndim))
  allocate(qp(1:ncells,1:nvar,1:ndim))
  allocate(unit_var(1:nvar))
  unit_var = 1.0d0
  ! Force on the gas
  allocate(force(1:ncells,1:ndim+nivz))
#if TURB==1
  allocate(f_turb(1:ncells,1:ndim))
  f_turb=0.0d0
  allocate(f_turb_old(1:ncells,1:ndim))
  f_turb_old=0.0d0
#endif
  allocate(cs(1:ncells))

#if NDUST>0
  ! Dust variables when included
  call allocate_dust
#endif

  u_prim   = 0.0d0
  q      = 0.0d0
  qm     = 0.0d0
  qp     = 0.0d0
  force  = 0.0d0
  cs     = 0.0d0
  !Indexation of variables /!\ Every index must be unique but the order does not matter
  irho= 1
  iv  = 2
  if(ndim==2) ivy=3  
#if IVZ==1  
   ivz = 4
#endif
  iP  = 2+nivz+ndim ! Energy equation
  print *,'Are the indices in order ? '

  print *,'irho =', irho
  print *,'iv   =', iv
  if(ndim==2)print *,'ivy  =', ivy
#if IVZ==1
  if(ndim==2)print *,'ivz  =', ivz
#endif 
  print *,'iP   =', iP

#if NDUST>0  
  do idust=1,ndust
    print *,'idust   =', idust, 'irhod = ',2+ndim+idust, ' ivd = ', 2+ndim+ndust+idust,' ivdy = ',  2+ndim+2*ndust+idust
    irhod(idust)= iP+idust
    ivd(idust)  = iP+ndust+idust
    ivdy(idust) = iP+2*ndust+idust
#if IVZ==1    
    ivdz(idust) = iP+3*ndust+idust
#endif
  end do
#endif
  print *, 'Nvar =', nvar


end subroutine allocate_init

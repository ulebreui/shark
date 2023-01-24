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

  ! Variable related quantities
  allocate(unew(1:ncells,1:nvar))
  allocate(uold(1:ncells,1:nvar))
  allocate(q(1:ncells,1:nvar))
  allocate(ql(1:ncells,1:nvar,1:ndim))
  allocate(qr(1:ncells,1:nvar,1:ndim))

  ! Force on the gas
  allocate(force(1:ncells))
  allocate(cs(1:ncells))

#if NDUST>0
  ! Dust variables when included
  call allocate_dust
#endif

  unew   = 0.0d0
  uold   = 0.0d0
  q      = 0.0d0
  ql     = 0.0d0
  qr     = 0.0d0
  force  = 0.0d0
  cs     = 0.0d0
  !Indexation of variables /!\ Every index must be unique but the order does not matter
  irho= 1
  iv  = 2
  if(ndim==2) ivy=3  
  iP  = 2+ndim ! Energy equation
  print *,'Are the indices in order ? '

  print *,'irho =', irho
  print *,'iv   =', iv
  if(ndim==2)print *,'ivy  =', ivy
  print *,'iP   =', iP
  print *, 'Nvar =', nvar

#if NDUST>0  
  do idust=1,ndust
     irhod(idust)= 2+ndim+idust
     ivd(idust)  = 2+ndim+ndust+idust
     ivdy(idust) = 2+ndim+2*ndust+idust
  end do
#endif


end subroutine allocate_init

!Allocation of all the quantities
subroutine allocate_init
  use parameters
  use commons
  use units
  implicit none

  integer :: idust


  ! Grid related quantities
  allocate(r(1:ncells))
  allocate(r_c(1:ncells))
  allocate(vol_cell(1:ncells))
  allocate(dvol(1:ncells))
  allocate(dx(1:ncells))
  allocate(dx_l(1:ncells))
  allocate(dx_r(1:ncells))
  allocate(dx_l_cell(1:ncells))
  allocate(dx_r_cell(1:ncells))
  allocate(dx_c(1:ncells))
  allocate(Surf_p(1:ncells))
  allocate(Surf_m(1:ncells))

  r       = 0.0d0
  r_c     = 0.d0
  vol_cell= 0.d0

  dx   = 0.d0
  dx_l = 0.d0
  dx_r = 0.d0
  dx_l_cell=0.d0
  dx_r_cell=0.d0
  dx_c=0.d0
  Surf_p=0.d0
  Surf_m=0.0d0

  ! Variable related quantities
  allocate(unew(1:ncells,1:nvar))
  allocate(uold(1:ncells,1:nvar))
  allocate(q(1:ncells,1:nvar))
  allocate(ql(1:ncells,1:nvar))
  allocate(qr(1:ncells,1:nvar))

  !Sound speed
  allocate(cs(1:ncells))
  allocate(csl(1:ncells))
  allocate(csr(1:ncells))

  ! Magnetic field
  allocate(B_cell(1:ncells))
  ! Force on the gas
  allocate(force(1:ncells))

#if NDUST>0
  ! Dust variables when included
  call allocate_dust
#endif

  unew  = 0.0d0
  uold  = 0.0d0
  q     = 0.0d0
  ql    = 0.0d0
  qr    = 0.0d0
  cs    = 0.0d0
  csl   = 0.0d0
  csr   = 0.0d0

  B_cell = 0.0d0
  force  = 0.0d0

#if GRAVITY==1  
  allocate(Mc(1:ncells))
  Mc= 0.0d0
#endif  

  !Indexation of variables /!\ Every index must be unique but the order does not matter
  irho=1
  iv=2
  iP=0
#if ENERGY>0
  iP=3 ! Energy equation not coded. It will be done soon by the next release
#endif
  
#if NDUST>0  
  do idust=1,ndust
     irhod(idust)= 2+idust+nenergy
     ivd(idust)  = 2+ndust+idust+nenergy
  end do
#endif


end subroutine allocate_init

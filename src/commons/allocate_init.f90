
!Allocation of all the quantities
subroutine allocate_init
  use parameters
  use commons
  use units
  implicit none

  integer :: idust,i,ipscal,icountpscal

  call get_active_cells

  ! Grid related quantities
  allocate(vol(1:nx_max,1:ny_max))
  allocate(surf(1:nx_max,1:ny_max,1:ndim))
  allocate(dx(1:nx_max,1:ny_max,1:ndim))



  allocate(position(1:nx_max,1:ny_max,1:ndim))
  allocate(radii_c(1:ncells))

  vol      = 0.d0
  Surf     = 0.d0
  dx       = 0.d0
  
  position = 0.0d0

  do i=1,ncells
    radii_c(i) = 1.0d0 ! We set it to 1 for cartesian geometry to avoid code duplication
  end do

#if GEOM==2
    allocate(phi(1:nx_max,1:ny_max))
    phi=1.0d0
#endif

  ! Variable related quantities
  allocate(q(1:nx_max,1:ny_max,1:nvar))
  allocate(u_prim(1:nx_max,1:ny_max,1:nvar))
  allocate(qm(1:nx_max,1:ny_max,1:nvar,1:ndim))
  allocate(qp(1:nx_max,1:ny_max,1:nvar,1:ndim))



  ! Force on the gas
  allocate(force(1:ncells,1:3))
  allocate(flux(1:ncells,1:nvar,1:ndim))
  flux=0.0d0

  allocate(cs(1:ncells))

#if NDUST>0
  ! Dust variables when included
  call allocate_dust
#endif

  u_prim = 0.0d0
  q      = 0.0d0
  qm     = 0.0d0
  qp     = 0.0d0
  force  = 0.0d0
  cs     = 0.0d0

  !Indexation of variables /!\ Every index must be unique but the order does not matter
  irho = 1
  ivx  = 2
  ivy  = 3  
  ivz  = 4
  iP   = 5 ! Energy equation
  print *,'Are the indices in order ? '
  print *,'irho =', irho
  print *,'ivx   =', ivx
  print *,'ivy  =', ivy
  print *,'ivz  =', ivz
  print *,'iP   =', iP
  index_vn(1) = ivx  
  index_vt(1) = ivy
  index_vn(2) = ivy
  index_vt(2) = ivx

  
#if NDUST>0  
  icountpscal=1
  do idust=1,ndust
    print *,'idust   =', idust, 'irhod = ',iP+idust, ' ivdx = ', iP+ndust+idust,' ivdy = ',  iP+2*ndust+idust,' ivdz = ',  iP+3*ndust+idust
    irhod(idust)= iP+idust
    ivdx(idust) = iP+ndust+idust
    ivdy(idust) = iP+2*ndust+idust
    ivdz(idust) = iP+3*ndust+idust
#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
      idust_pscal(idust,ipscal) = iP +4*ndust+icountpscal
      print *,'idustpscal   =', idust_pscal(idust,ipscal)
      icountpscal= icountpscal +1
    end do
#endif     
    index_vdn(idust,1) = ivdx(idust)
    index_vdt(idust,1) = ivdy(idust)
    index_vdn(idust,2) = ivdy(idust)
    index_vdt(idust,2) = ivdx(idust)
  end do
#endif
  smallr=smallr/unit_d
  smallc=smallc/unit_v
  smallp=smallr*smallc**2
  print *, 'Nvar =', nvar


end subroutine allocate_init

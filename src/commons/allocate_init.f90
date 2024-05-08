
!Allocation of all the quantities
subroutine allocate_init
  use parameters
  use commons
  use units
  implicit none

  integer :: idust,i,ipscal,icountpscal

  call get_active_cells

  ! Grid related quantities
  allocate(vol(1:ncells))
  allocate(dvol(1:ncells))
  allocate(dx(1:ncells,1:ndim))
  allocate(surf(1:ncells,1:ndim))
  allocate(position(1:ncells,1:ndim))
  allocate(radii(1:ncells))
  allocate(radii_c(1:ncells))

  vol      = 0.d0
  dvol     = 0.d0
  dx       = 0.d0
  Surf     = 0.d0
  position = 0.0d0

  do i=1,ncells
    radii(i)   = 1.0d0 ! We set it to 1 for cartesian geometry to avoid code duplication
    radii_c(i) = 1.0d0 ! We set it to 1 for cartesian geometry to avoid code duplication
  end do

#if GEOM==1
  radii=0.0d0 
  radii_c=0.0d0
  allocate(theta(1:ncells))
  theta=0.0d0 
  allocate(dx_l(1:ncells))
  allocate(dx_r(1:ncells))
  allocate(dx_l_cell(1:ncells))
  allocate(dx_r_cell(1:ncells))
  allocate(dx_c(1:ncells))
  dx_l = 0.d0
  dx_r = 0.d0
  dx_l_cell=0.d0
  dx_r_cell=0.d0
  dx_c=0.d0
#endif
#if GEOM==2
    allocate(phi(1:ncells))
    phi=1.0d0
#if GRIDSPACE==1
    allocate(dx_l_cell(1:nx_max))
    allocate(dx_r_cell(1:nx_max))
    dx_l_cell=0.d0
    dx_r_cell=0.d0
#endif
#endif
  ! Variable related quantities
  allocate(u_prim(1:ncells,1:nvar))
  allocate(q(1:ncells,1:nvar))
  allocate(qm(1:ncells,1:nvar,1:ndim))
  allocate(qp(1:ncells,1:nvar,1:ndim))
  ! Force on the gas
  allocate(force(1:ncells,1:3))
  allocate(flux(1:ncells,1:nvar,1:ndim))
  flux=0.0d0
  !Gravitational potential
  allocate(phi_sg(1:ncells))
  allocate(grad_phi_sg(1:ncells,1:2))
  phi_sg=0.0d0
  grad_phi_sg=0.0d0
#if TURB==1
  allocate(f_turb(1:ncells,1:3))
  f_turb=0.0d0
  allocate(f_turb_old(1:ncells,1:3))
  f_turb_old=0.0d0
#endif
  allocate(cs(1:ncells))

#if NDUST>0
  ! Dust variables when included
  call allocate_dust
#endif

#if GRAVITY==1  
  allocate(Mc(1:ncells))
  Mc= 0.0d0
#endif  

  if(fargo) then
    allocate(fargo_velocity(1:ncells))
    fargo_velocity = 0.0d0
  endif
  if(charging) then
    allocate(eta_a(1:ncells))
    allocate(eta_h(1:ncells))
    allocate(eta_o(1:ncells))
    allocate(sigma_o(1:ncells))
    allocate(sigma_p(1:ncells))
    allocate(sigma_h(1:ncells))
    allocate(ni(1:ncells))
    allocate(ne(1:ncells))
    allocate(psi_old(1:ncells))

    eta_a  =0.0d0
    eta_h  =0.0d0
    eta_o  =0.0d0
    sigma_o=0.0d0
    sigma_p=0.0d0
    sigma_h=0.0d0
    ni     =0.0d0
    ne     =0.0d0
    psi_old=0.0d0

  endif

  u_prim = 0.0d0
  q      = 0.0d0
  qm     = 0.0d0
  qp     = 0.0d0
  force  = 0.0d0
  cs     = 0.0d0
  allocate(eta_visc(1:ncells))
  eta_visc=0.0d0
  !Indexation of variables /!\ Every index must be unique but the order does not matter
  irho = 1
  ivx = 2
  ivy = 3  
  ivz = 4
#if MHD==1
  iBx=5
  iBy=6
  iBz=7
#endif  
  iP  = 5 + Nmhd! Energy equation
  print *,'Are the indices in order ? '

!To display B TODO
  print *,'irho =', irho
  print *,'ivx   =', ivx
  print *,'ivy  =', ivy
  print *,'ivz  =', ivz
#if MHD==1
  print *,'iBx  =', iBx
  print *,'iBy  =', iBy
  print *,'iBz  =', iBz
#endif

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

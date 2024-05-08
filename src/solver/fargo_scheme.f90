!This is the fargo orbital scheme : the solver is largely adapted from the DUMSES-hybrid solver.
subroutine fargo_scheme
  use parameters
  use commons
  use units
  use OMP_LIB
  use slope_limiter
  integer :: i,ix,iy,il,ir,nshft,jshft,jshftp,iprim,iprimprim,ivar
  real(dp),dimension(:,:), allocatable :: uprim_old,duprim
  real(dp) :: radius_polar,vk,dy_loc,forth
  allocate(uprim_old(1:ncells,1:nvar))

  allocate(duprim(1:ncells,1:nvar))

  call apply_boundaries

  uprim_old = u_prim
  forth     = 0.25d0
  !To do
  do ix = 1, nx
    do iy = 1, ny
      i = icell(ix,iy)
      if(active_cell_predictor(i)==1) then
        il = icell(ix,iy-1)
        ir = icell(ix,iy+1)
        do ivar=1,nvar
          duprim(i,ivar) = slope_limit((u_prim(i,ivar) - u_prim(il,ivar)),(u_prim(ir,ivar) - u_prim(i,ivar)))
        end do
      endif 
    end do
  end do 

    !To do
  do ix = 1, nx_max
    do iy = 1, ny_max
      i = icell(ix,iy)
      if(active_cell(i)==1) then
#if GEOM<2
      radius_polar=1.0d0
#endif      
#if GEOM==2
      radius_polar=radii_c(i)
#endif
      vk     = fargo_velocity(i)
      dy_loc = dx(i,2)*radius_polar
      nshft  = int(abs(vK*dt)/dy_loc)
      eps    = mod(abs(vK*dt),dy_loc)/dy_loc
      if (vk > 0.0d0) then
         jshft = iy - nshft - 1
         if (jshft < nghost+1) jshft  = jshft  + ny
         jshftp = jshft + 1
         if (jshftp< nghost+1) jshftp = jshftp + ny
         iprim     = icell(ix,jshft)
         iprimprim = icell(ix,jshftp)

         u_prim(i,1:nvar) = uprim_old(iprim,1:nvar)*eps &
                              & + uprim_old(iprimprim,1:nvar)*(one - eps)
         u_prim(i,1:nvar) = u_prim(i,1:nvar) &
             & + half*duprim(iprim,1:nvar)*(forth - (half - eps)**2)&
             & + half*duprim(iprimprim,1:nvar)*((half - eps)**2 - forth)
      else
         jshft  = iy + nshft
         if (jshft > ny + nghost) jshft   = jshft - ny
         jshftp = jshft + 1
         if (jshftp > ny + nghost) jshftp = jshftp - ny
         iprim     = icell(ix,jshft)
         iprimprim = icell(ix,jshftp)
         !print *, ix, iy, jshft,jshftp,iprim,iprimprim
         u_prim(i,1:nvar) = uprim_old(iprim,1:nvar)*(one - eps) &
                              & + uprim_old(iprimprim,1:nvar)*eps
         u_prim(i,1:nvar) = u_prim(i,1:nvar) &
             & + half*duprim(iprim,1:nvar)*(forth - (half - eps)**2) &
             & + half*duprim(iprimprim,1:nvar)*((half - eps)**2 - forth)
      endif
      endif 
    end do
  end do 
  
  !u_prim = uprim_old
  call apply_boundaries
end subroutine fargo_scheme


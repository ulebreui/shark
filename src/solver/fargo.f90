#if GEOM==2
! Mean azimuthal gas velocity of each ring (from q, after ctoprim)
subroutine compute_vbar
  use parameters
  use commons
  implicit none
  integer  :: ix, iy
  real(dp) :: s
  !$omp parallel do default(shared) schedule(static) private(ix,iy,s)
  do ix = 1, nx_max
     s = 0.0d0
     do iy = first_active_y, last_active_y
        s = s + q(ivy,ix,iy)
     end do
     vbar(ix) = s/dble(ny)
  end do
  !$omp end parallel do
end subroutine compute_vbar

! Conservative shift of each ring by vbar*dt: integer part + 2nd-order upwind fractional part
subroutine fargo_shift
  use parameters
  use commons
  use slope_limiter
  implicit none
  integer  :: ix, ivar, nshift, j
  real(dp) :: delta, f
  real(dp), dimension(0:ny-1) :: ring, flx

  !$omp parallel do default(shared) schedule(static) private(ix,ivar,nshift,j,delta,f,ring,flx)
  do ix = first_active, last_active
     delta  = vbar(ix)*dt/(radii(ix,first_active_y)*dx(ix,first_active_y,2)) ! shift in cells
     nshift = floor(delta)
     f      = delta - dble(nshift)                                          ! 0 <= f < 1
     do ivar = 1, nvar
        do j = 0, ny-1                                  ! integer shift (periodic)
           ring(j) = u_prim(ivar,ix,first_active_y + modulo(j - nshift, ny))
        end do
        do j = 0, ny-1                                  ! flux through the +phi face of cell j
           flx(j) = f*(ring(j) + half*(1.0d0 - f)* &
                &   slope_limit(ring(j) - ring(modulo(j-1,ny)), ring(modulo(j+1,ny)) - ring(j)))
        end do
        do j = 0, ny-1
           u_prim(ivar,ix,first_active_y + j) = ring(j) - flx(j) + flx(modulo(j-1,ny))
        end do
     end do
  end do
  !$omp end parallel do
end subroutine fargo_shift
#endif
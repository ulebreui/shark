! Disk on a spherical (r,theta) wedge: all ghost cells fixed to the initial state
subroutine apply_boundaries
  use parameters
  use commons
  use units
  implicit none
  integer :: ix, iy, k, ivar

  ! theta ghosts: density profile of the IC scaled to the last active cell, v_theta reflecting
  do ix = 1, nx_max
     do k = 1, nghost
        call th_ghost(ix, first_active_y - k, first_active_y, first_active_y + k - 1)
        call th_ghost(ix, last_active_y  + k, last_active_y,  last_active_y  - k + 1)
     end do
  end do

  ! r ghosts, active theta rows
  do iy = first_active_y, last_active_y
     do k = 1, nghost
        do ivar = 1, nvar
           u_prim(ivar,first_active-k,iy) = uprim_condinit(ivar,first_active-k,iy)
           u_prim(ivar,last_active +k,iy) = uprim_condinit(ivar,last_active +k,iy)
        end do
     end do
  end do
end subroutine apply_boundaries

! ghost (ix,ig): rho = rho_IC(ig) * rho(ia)/rho_IC(ia), v_r and v_phi from the IC, v_theta = -v_theta(im)
! ia = first/last active row, im = mirror active row of the ghost
subroutine th_ghost(ix, ig, ia, im)
  use parameters
  use commons
  implicit none
  integer, intent(in) :: ix, ig, ia, im
  real(dp) :: f, rho_g
  f     = u_prim(irho,ix,ia)/uprim_condinit(irho,ix,ia)
  rho_g = uprim_condinit(irho,ix,ig)*f
  u_prim(irho,ix,ig) = rho_g
  u_prim(ivx ,ix,ig) = uprim_condinit(ivx,ix,ig)*f                       ! rho*v_r  (v_r from IC)
  u_prim(ivz ,ix,ig) = uprim_condinit(ivz,ix,ig)*f                       ! rho*v_phi (v_phi from IC)
  u_prim(ivy ,ix,ig) = -u_prim(ivy,ix,im)/u_prim(irho,ix,im)*rho_g        ! reflecting v_theta
  u_prim(iP  ,ix,ig) = uprim_condinit(iP,ix,ig)*f                        ! isothermal: unused
end subroutine th_ghost
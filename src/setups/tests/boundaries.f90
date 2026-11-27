! Boundaries of the 1D validation tests
!   x: zero-gradient (sod), fixed to the initial state (dustyshock), periodic (others)
!   y: periodic
subroutine apply_boundaries
  use parameters
  use commons
  use units
  implicit none
  integer :: ix, iy, ivar

  do iy = 1, ny_max
     do ix = 1, nghost
        do ivar = 1, nvar
           select case (trim(test_case))
           case ('sod')
              u_prim(ivar,ix,iy)          = u_prim(ivar,first_active,iy)
              u_prim(ivar,nx_max+1-ix,iy) = u_prim(ivar,last_active,iy)
           case ('dustyshock')
              u_prim(ivar,ix,iy)          = uprim_condinit(ivar,ix,iy)
              u_prim(ivar,nx_max+1-ix,iy) = uprim_condinit(ivar,nx_max+1-ix,iy)
           case default
              u_prim(ivar,ix,iy)          = u_prim(ivar,last_active-nghost+ix,iy)
              u_prim(ivar,nx_max+1-ix,iy) = u_prim(ivar,first_active+nghost-ix,iy)
           end select
        end do
     end do
  end do

  do ix = 1, nx_max
     do iy = 1, nghost
        do ivar = 1, nvar
           u_prim(ivar,ix,iy)          = u_prim(ivar,ix,last_active_y-nghost+iy)
           u_prim(ivar,ix,ny_max+1-iy) = u_prim(ivar,ix,first_active_y+nghost-iy)
        end do
     end do
  end do

end subroutine apply_boundaries

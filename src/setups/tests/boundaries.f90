! Boundaries of the 1D validation tests. Along the direction of the problem (direction = 1: x, 2: y):
!   zero-gradient (sod), fixed to the initial state (dustyshock), periodic (others).
! Along the other direction: periodic.
subroutine apply_boundaries
  use parameters
  use commons
  use units
  implicit none
  integer :: ix, iy, ivar
  character(len=10) :: bcx, bcy, bc

  select case (trim(test_case))
  case ('sod');        bc = 'zerograd'
  case ('dustyshock'); bc = 'fixed'
  case default;        bc = 'periodic'
  end select
  bcx = 'periodic'; bcy = 'periodic'
  if (direction == 1) bcx = bc
  if (direction == 2) bcy = bc

  ! x boundaries
  do iy = 1, ny_max
     do ix = 1, nghost
        do ivar = 1, nvar
           select case (trim(bcx))
           case ('zerograd')
              u_prim(ivar,ix,iy)          = u_prim(ivar,first_active,iy)
              u_prim(ivar,nx_max+1-ix,iy) = u_prim(ivar,last_active,iy)
           case ('fixed')
              u_prim(ivar,ix,iy)          = uprim_condinit(ivar,ix,iy)
              u_prim(ivar,nx_max+1-ix,iy) = uprim_condinit(ivar,nx_max+1-ix,iy)
           case default
              u_prim(ivar,ix,iy)          = u_prim(ivar,last_active-nghost+ix,iy)
              u_prim(ivar,nx_max+1-ix,iy) = u_prim(ivar,first_active+nghost-ix,iy)
           end select
        end do
     end do
  end do

  ! y boundaries
  do ix = 1, nx_max
     do iy = 1, nghost
        do ivar = 1, nvar
           select case (trim(bcy))
           case ('zerograd')
              u_prim(ivar,ix,iy)          = u_prim(ivar,ix,first_active_y)
              u_prim(ivar,ix,ny_max+1-iy) = u_prim(ivar,ix,last_active_y)
           case ('fixed')
              u_prim(ivar,ix,iy)          = uprim_condinit(ivar,ix,iy)
              u_prim(ivar,ix,ny_max+1-iy) = uprim_condinit(ivar,ix,ny_max+1-iy)
           case default
              u_prim(ivar,ix,iy)          = u_prim(ivar,ix,last_active_y-nghost+iy)
              u_prim(ivar,ix,ny_max+1-iy) = u_prim(ivar,ix,first_active_y+nghost-iy)
           end select
        end do
     end do
  end do

end subroutine apply_boundaries

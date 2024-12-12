! Force kick (if applicable)
subroutine kick
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: idust,ix,iy
  real(dp) :: energy_old,energy_new,coeffdt

  if(static) return
  if(iso_cs.ne.1 .or. non_standard_eos .ne. 1) then
    !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy,energy_old,energy_new)
    do iy=first_active_y,last_active_y
      do ix=first_active,last_active
          energy_old         = half*u_prim(ivx,ix,iy)**2/u_prim(irho,ix,iy) +  half*u_prim(ivy,ix,iy)**2/u_prim(irho,ix,iy)  + half*u_prim(ivz,ix,iy)**2/u_prim(irho,ix,iy)
          u_prim(ivx,ix,iy)  = u_prim(ivx,ix,iy)  + u_prim(irho,ix,iy)*dt*force_x(ix,iy)
          u_prim(ivy,ix,iy)  = u_prim(ivy,ix,iy)  + u_prim(irho,ix,iy)*dt*force_y(ix,iy)
          u_prim(ivz,ix,iy)  = u_prim(ivz,ix,iy)  + u_prim(irho,ix,iy)*dt*force_z(ix,iy)

          energy_new         = half*u_prim(ivx,ix,iy)**2/u_prim(irho,ix,iy) +  half*u_prim(ivy,ix,iy)**2/u_prim(irho,ix,iy) + half*u_prim(ivz,ix,iy)**2/u_prim(irho,ix,iy)
          u_prim(iP,ix,iy)   = u_prim(iP,ix,iy) +  (energy_new - energy_old)
#if NDUST>0
        do idust=1,ndust
          u_prim(ivdx(idust),ix,iy)            = u_prim(ivdx(idust),ix,iy) + u_prim(irhod(idust),ix,iy)*dt*force_dust_x(idust,ix,iy)
          u_prim(ivdy(idust),ix,iy)            = u_prim(ivdy(idust),ix,iy) + u_prim(irhod(idust),ix,iy)*dt*force_dust_y(idust,ix,iy)
          u_prim(ivdz(idust),ix,iy)            = u_prim(ivdz(idust),ix,iy) + u_prim(irhod(idust),ix,iy)*dt*force_dust_z(idust,ix,iy)    
        end do
#endif
    end do
  end do

else

  !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ix,iy)
  do iy=first_active_y,last_active_y
    do ix=first_active,last_active
        u_prim(ivx,ix,iy)  = u_prim(ivx,ix,iy)  + u_prim(irho,ix,iy)*dt*force_x(ix,iy)
        u_prim(ivy,ix,iy)  = u_prim(ivy,ix,iy)  + u_prim(irho,ix,iy)*dt*force_y(ix,iy)
        u_prim(ivz,ix,iy)  = u_prim(ivz,ix,iy)  + u_prim(irho,ix,iy)*dt*force_z(ix,iy)
#if NDUST>0
      do idust=1,ndust
        u_prim(ivdx(idust),ix,iy)            = u_prim(ivdx(idust),ix,iy) + u_prim(irhod(idust),ix,iy)*dt*force_dust_x(idust,ix,iy)
        u_prim(ivdy(idust),ix,iy)            = u_prim(ivdy(idust),ix,iy) + u_prim(irhod(idust),ix,iy)*dt*force_dust_y(idust,ix,iy)
        u_prim(ivdz(idust),ix,iy)            = u_prim(ivdz(idust),ix,iy) + u_prim(irhod(idust),ix,iy)*dt*force_dust_z(idust,ix,iy)    
      end do
#endif
  end do
end do
endif
end subroutine kick


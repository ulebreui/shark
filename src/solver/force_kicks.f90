! Force kick (if applicable)
subroutine kick(coeffdt)
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: idust,ix,iy,icell
  real(dp) :: energy_old,energy_new,coeffdt

  if(static) return

  do iy=first_active_y,last_active_y
   do ix=first_active,last_active
      if(iso_cs.ne.1 .or. non_standard_eos .ne. 1) then
        energy_old         = half*u_prim(ix,iy,ivx)**2/u_prim(ix,iy,irho) +  half*u_prim(ix,iy,ivy)**2/u_prim(ix,iy,irho)  + half*u_prim(ix,iy,ivz)**2/u_prim(ix,iy,irho)
      endif
        u_prim(ix,iy,ivx)  = u_prim(ix,iy,ivx)  + u_prim(ix,iy,irho)*coeffdt*dt*force_x(ix,iy)
        u_prim(ix,iy,ivy)  = u_prim(ix,iy,ivy)  + u_prim(ix,iy,irho)*coeffdt*dt*force_y(ix,iy)
        u_prim(ix,iy,ivz)  = u_prim(ix,iy,ivz)  + u_prim(ix,iy,irho)*coeffdt*dt*force_z(ix,iy)

      if(iso_cs.ne.1 .or. non_standard_eos .ne. 1) then
        energy_new         = half*u_prim(ix,iy,ivx)**2/u_prim(ix,iy,irho) +  half*u_prim(ix,iy,ivy)**2/u_prim(ix,iy,irho) + half*u_prim(ix,iy,ivz)**2/u_prim(ix,iy,irho)
        u_prim(ix,iy,iP)   = u_prim(ix,iy,iP) +  (energy_new - energy_old)
      endif
#if NDUST>0
      do idust=1,ndust
        u_prim(ix,iy,ivdx(idust))            = u_prim(ix,iy,ivdx(idust)) + u_prim(ix,iy,irhod(idust))*coeffdt*dt*force_dust_x(ix,iy,idust)
        u_prim(ix,iy,ivdy(idust))            = u_prim(ix,iy,ivdy(idust)) + u_prim(ix,iy,irhod(idust))*coeffdt*dt*force_dust_y(ix,iy,idust)
        u_prim(ix,iy,ivdz(idust))            = u_prim(ix,iy,ivdz(idust)) + u_prim(ix,iy,irhod(idust))*coeffdt*dt*force_dust_z(ix,iy,idust)    
      end do
#endif
  end do
end do


end subroutine kick


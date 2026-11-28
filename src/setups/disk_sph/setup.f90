! Axisymmetric disk on a spherical (r,theta) grid (GEOM=3). Empty skeleton: uniform gas at rest.
subroutine setup
  use parameters
  use commons
  use units
  implicit none
  integer :: ix, iy
  real(dp) :: rr, RR_cyl, cs2, hr2
  real(dp) :: rd_in, rd_out
  call allocate_init
  call gridinit_sph(r_in, r_out, th_min, th_max)
  allocate(uprim_condinit(1:nvar,1:nx_max,1:ny_max))

  q = 0.0d0
  iso_cs = 1
  do iy = 1, ny_max
     do ix = 1, nx_max
        rr  = radii(ix,iy)
        RR_cyl = rr*sin(phi(ix,iy))
        cs2 = (HoverR**2*Mstar/r0)*(RR_cyl/r0)**q_temp
        hr2 = cs2*RR_cyl/Mstar                                  ! (H/R)^2 = cs^2/vK^2
        cs(ix,iy)     = sqrt(cs2)
        q(irho,ix,iy) = rho0/unit_d*(RR_cyl/r0)**p_dens*exp(Mstar/cs2*(1.0d0/rr - 1.0d0/RR_cyl))
        q(ivx ,ix,iy) = 0.0d0
        q(ivy ,ix,iy) = 0.0d0
        q(ivz ,ix,iy) = sqrt(Mstar/RR_cyl)*sqrt((p_dens+q_temp)*hr2 + 1.0d0 + q_temp - q_temp*RR_cyl/rr)
        q(iP  ,ix,iy) = q(irho,ix,iy)*cs2
     end do
  end do

  call primtoc
  uprim_condinit = u_prim
  allocate(damp_rate(1:nx_max,1:ny_max)); damp_rate = 0.0d0
  if (damping) then
     rd_in  = damp_in_frac*r_in
     rd_out = damp_out_frac*r_out
     do iy = 1, ny_max
        do ix = 1, nx_max
           rr = radii(ix,iy)
           if (rr < rd_in)  damp_rate(ix,iy) = ((rd_in - rr)/(rd_in - r_in))**2 &
                &                              /(damp_tau*2.0d0*pi*sqrt(r_in**3/Mstar))
           if (rr > rd_out) damp_rate(ix,iy) = ((rr - rd_out)/(r_out - rd_out))**2 &
                &                              /(damp_tau*2.0d0*pi*sqrt(r_out**3/Mstar))
        end do
     end do
  end if
#if NDUST>0
  call compute_tstop
#endif
end subroutine setup


subroutine get_rhoturb(pert,del)
  use random
  use precision
  implicit none
  integer :: iseed = 0
  integer, dimension(1,1:IRandNumSize) :: allseed
  integer, dimension(IRandNumSize)     :: localseed = -1
  real(dp) :: pert, del, randno
  if (localseed(1) == -1) then
     call rans(1,iseed,allseed)
     localseed = allseed(1,1:IRandNumSize)
  end if
  call ranf(localseed,randno)
  del = randno*pert - 0.5d0*pert
end subroutine get_rhoturb


subroutine write_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer :: ilun
  write(ilun,*) time
end subroutine write_setup_info


subroutine read_setup_info(ilun)
  use parameters
  use commons
  use units
  implicit none
  integer  :: ilun
  real(dp) :: info1
  read(ilun,*) info1
  time = info1
  close(ilun)
end subroutine read_setup_info


subroutine read_setup_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70) :: nmlfile
  character(len=200) :: msg
  integer :: io, ilun
  namelist/setup_params/r_in,r_out,th_min,th_max,Mstar,HoverR,r0,rho0,p_dens,q_temp,damping,damp_in_frac,damp_out_frac,damp_tau
  msg = ''
  read(13,setup_params,IOSTAT=io,IOMSG=msg)
  rewind(13)
  if (io /= 0) then
     write(*,*) 'Invalid line in the setup namelist: io =', io, ' ', trim(msg)
     stop
  end if
end subroutine read_setup_params


subroutine flag_continue(continue_sim)
  use parameters
  use commons
  use units
  implicit none
  logical :: continue_sim
  if (time >= tend) continue_sim = .false.
end subroutine flag_continue


subroutine check_output(icount,iout,outputing,verbose)
  use parameters
  use commons
  use units
  implicit none
  integer :: icount, iout
  logical :: outputing, verbose
  if (icount .eq. freq_out) then
     verbose    = .true.
     outputing  = .true.
     icount     = 0
  end if
  if (time .eq. 0) outputing = .true.
  if (outputing) then
     call output(iout)
     iout = iout + 1
  end if
  outputing = .false.
end subroutine check_output


subroutine setup_preloop
  use parameters
  use commons
  use units
  implicit none
end subroutine setup_preloop


subroutine setup_inloop
  use parameters
  use commons
  use units
  implicit none
  integer  :: ix, iy, ivar
  real(dp) :: f1
  if (.not. damping) return
  !$omp parallel do default(shared) schedule(static) private(ix,iy,ivar,f1)
  do iy = first_active_y, last_active_y
     do ix = first_active, last_active
        if (damp_rate(ix,iy) > 0.0d0) then
           f1 = dt*damp_rate(ix,iy)/(1.0d0 + dt*damp_rate(ix,iy))
           do ivar = 1, nvar
              u_prim(ivar,ix,iy) = u_prim(ivar,ix,iy) - f1*(u_prim(ivar,ix,iy) - uprim_condinit(ivar,ix,iy))
           end do
        end if
     end do
  end do
  !$omp end parallel do
end subroutine setup_inloop


subroutine update_force_setup
  use parameters
  use commons
  use units
  implicit none
  integer :: ix, iy, idust
  !$omp parallel do default(shared) schedule(static) private(ix,iy,idust)
  do iy = 1, ny_max
     do ix = 1, nx_max
        force_x(ix,iy) = -Mstar/radii(ix,iy)**2
        force_y(ix,iy) = 0.0d0
        force_z(ix,iy) = 0.0d0
#if NDUST>0
        do idust = 1, ndust
           force_dust_x(idust,ix,iy) = force_x(ix,iy)
           force_dust_y(idust,ix,iy) = 0.0d0
           force_dust_z(idust,ix,iy) = 0.0d0
        end do
#endif
     end do
  end do
  !$omp end parallel do
end subroutine update_force_setup


#if NDUST>0
! TODO: Epstein stopping time
subroutine compute_tstop
  use parameters
  use commons
  use units
  implicit none
  tstop = 1.0d30
end subroutine compute_tstop
#endif


! Gas scale height (only used by the vertically integrated kernel, GEOM=2)
function gas_scale_height(ix,iy)
  use parameters
  use commons
  implicit none
  integer  :: ix, iy
  real(dp) :: gas_scale_height
  gas_scale_height = cs(ix,iy)*(radii(ix,iy)*sin(phi(ix,iy)))**1.5d0/sqrt(Mstar)
end function gas_scale_height

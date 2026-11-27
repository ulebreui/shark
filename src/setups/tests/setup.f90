! Standard 1D validation tests (Cartesian, GEOM=0, NY=1). See setup_commons.f90 for the list.
subroutine setup
  use parameters
  use commons
  use units
  implicit none
  integer  :: ix, iy, idust
  real(dp) :: xx, kw, sn, mg, md, pitot, vdn

  call allocate_init
  call gridinit(box_l, box_l_y)
  allocate(uprim_condinit(1:nvar,1:nx_max,1:ny_max))

  q = 0.0d0
  iso_cs = 1
  if (trim(test_case) == 'sod') iso_cs = 0
  kw = 2.0d0*pi*nwave/box_l

  ! Dusty shock: downstream state from mass and total momentum conservation
  mg  = v_shock; md = eps_shock*v_shock
  pitot = mg*v_shock + cs0**2 + md*v_shock
  vdn = (pitot - sqrt(pitot**2 - 4.0d0*(mg + md)*cs0**2*mg))/(2.0d0*(mg + md))

  do iy = 1, ny_max
     do ix = 1, nx_max
        xx = position(ix,iy,1)
        sn = sin(kw*xx)
        cs(ix,iy) = cs0
        select case (trim(test_case))
        case ('sod')
           if (xx < x_disc) then
              q(irho,ix,iy) = sod_rho_l; q(iP,ix,iy) = sod_p_l
           else
              q(irho,ix,iy) = sod_rho_r; q(iP,ix,iy) = sod_p_r
           end if
           cs(ix,iy) = sqrt(gamma*q(iP,ix,iy)/q(irho,ix,iy))
        case ('soundwave')                       ! right-moving wave
           q(irho,ix,iy) = 1.0d0 + amp*sn
           q(ivx,ix,iy)  = cs0*amp*sn
        case ('dustybox')
           q(irho,ix,iy) = 1.0d0
           q(ivx,ix,iy)  = vg0
        case ('dustywave')
           q(irho,ix,iy) = 1.0d0 + amp*sn
           q(ivx,ix,iy)  = cs0*amp*sn
        case ('dustyshock')
           if (xx < x_disc) then
              q(irho,ix,iy) = 1.0d0;  q(ivx,ix,iy) = v_shock
           else
              q(irho,ix,iy) = mg/vdn; q(ivx,ix,iy) = vdn
           end if
        case default
           print *, 'Unknown test_case: ', trim(test_case)
           stop
        end select
        if (iso_cs == 1) q(iP,ix,iy) = q(irho,ix,iy)*cs0**2

#if NDUST>0
        do idust = 1, ndust
           select case (trim(test_case))
           case ('dustybox')
              q(irhod(idust),ix,iy) = eps_dust(idust)
              q(ivdx(idust),ix,iy)  = vd0(idust)
           case ('dustywave')
              q(irhod(idust),ix,iy) = eps_dust(idust)*(1.0d0 + amp*sn)
              q(ivdx(idust),ix,iy)  = cs0*amp*sn
           case ('dustyshock')
              if (xx < x_disc) then
                 q(irhod(idust),ix,iy) = eps_shock;  q(ivdx(idust),ix,iy) = v_shock
              else
                 q(irhod(idust),ix,iy) = md/vdn;     q(ivdx(idust),ix,iy) = vdn
              end if
           case default
              print *, 'test_case ', trim(test_case), ' must be compiled with NDUST=0'
              stop
           end select
        end do
#endif
     end do
  end do

  call primtoc
  uprim_condinit = u_prim
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
  integer :: io, ilun
  namelist/setup_params/test_case,box_l,box_l_y,cs0,amp,nwave,sod_rho_l,sod_p_l,sod_rho_r,sod_p_r,x_disc, &
       & eps_dust,ts_dust,vd0,vg0,v_shock,eps_shock,K_drag
  read(13,setup_params,IOSTAT=io)
  rewind(13)
  if (io /= 0) then
     write(*,*) 'Invalid line in the setup namelist'
     stop
  end if
  print *, 'Validation test: ', trim(test_case)
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


! Adiabatic tests: sound speed updated from the new state (used by courant)
subroutine setup_inloop
  use parameters
  use commons
  use units
  implicit none
  integer  :: ix, iy
  real(dp) :: rho, ekin, pres
  if (iso_cs == 1) return
  !$omp parallel do default(shared) schedule(static) private(ix,iy,rho,ekin,pres)
  do iy = first_active_y, last_active_y
     do ix = first_active, last_active
        rho  = u_prim(irho,ix,iy)
        ekin = half*(u_prim(ivx,ix,iy)**2 + u_prim(ivy,ix,iy)**2 + u_prim(ivz,ix,iy)**2)/rho
        pres = max((gamma - 1.0d0)*(u_prim(iP,ix,iy) - ekin), smallp)
        cs(ix,iy) = sqrt(gamma*pres/rho)
     end do
  end do
  !$omp end parallel do
end subroutine setup_inloop


! No external force in these tests
subroutine update_force_setup
  use parameters
  use commons
  use units
  implicit none
  force_x = 0.0d0; force_y = 0.0d0; force_z = 0.0d0
#if NDUST>0
  force_dust_x = 0.0d0; force_dust_y = 0.0d0; force_dust_z = 0.0d0
#endif
end subroutine update_force_setup


#if NDUST>0
! Stopping time: constant (dustybox, dustywave) or constant drag coefficient K (dustyshock: t_s = rho_d/K)
subroutine compute_tstop
  use parameters
  use commons
  use units
  implicit none
  integer :: idust, ix, iy
  !$omp parallel do default(shared) schedule(static) private(idust,ix,iy)
  do iy = 1, ny_max
     do ix = 1, nx_max
        do idust = 1, ndust
           if (trim(test_case) == 'dustyshock') then
              tstop(idust,ix,iy) = max(q(irhod(idust),ix,iy), smallr)/K_drag
           else
              tstop(idust,ix,iy) = ts_dust(idust)
           end if
        end do
     end do
  end do
  !$omp end parallel do
end subroutine compute_tstop
#endif

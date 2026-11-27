! Standard 1D validation tests (Cartesian, GEOM=0), along x (NY=1, direction=1) or y (NX=1, direction=2).
! See setup_commons.f90 for the list.
subroutine setup
  use parameters
  use commons
  use units
  implicit none
  integer  :: ix, iy, idust, iv, ivdd
  real(dp) :: xx, kw, sn, mg, md, pitot, vdn, lbox

  call allocate_init
#if GEOM==2
  call gridinit(box_l, box_l_y, r_in_test)       ! polar grid (disk_kernel test)
#else
  call gridinit(box_l, box_l_y)
#endif
  allocate(uprim_condinit(1:nvar,1:nx_max,1:ny_max))

  q = 0.0d0
  iso_cs = 1
  if (trim(test_case) == 'sod') iso_cs = 0
  ! Direction of the 1D problem: coordinate, length and velocity component
  lbox = box_l; iv = ivx
  if (direction == 2) then
     lbox = box_l_y; iv = ivy
  end if
  kw = 2.0d0*pi*nwave/lbox
#if NDUST>0
  if (trim(test_case) == 'smoluchowski' .or. trim(test_case) == 'fragmentation' .or. trim(test_case) == 'disk_kernel') &
       & call distribution_dust   ! bins + IC
#endif

  ! Dusty shock: downstream state from mass and total momentum conservation
  mg  = v_shock; md = eps_shock*v_shock
  pitot = mg*v_shock + cs0**2 + md*v_shock
  vdn = (pitot - sqrt(pitot**2 - 4.0d0*(mg + md)*cs0**2*mg))/(2.0d0*(mg + md))

  do iy = 1, ny_max
     do ix = 1, nx_max
        xx = position(ix,iy,direction)
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
           q(iv,ix,iy)   = cs0*amp*sn
        case ('dustybox')
           q(irho,ix,iy) = 1.0d0
           q(iv,ix,iy)   = vg0
        case ('dustywave')
           q(irho,ix,iy) = 1.0d0 + amp*sn
           q(iv,ix,iy)   = cs0*amp*sn
        case ('smoluchowski', 'disk_kernel', 'pscal_advect')   ! uniform gas at rest
           q(irho,ix,iy) = 1.0d0
        case ('fragmentation')                   ! uniform gas at rest (cgs units)
           q(irho,ix,iy) = rho_frag
        case ('dustyshock')
           if (xx < x_disc) then
              q(irho,ix,iy) = 1.0d0;  q(iv,ix,iy) = v_shock
           else
              q(irho,ix,iy) = mg/vdn; q(iv,ix,iy) = vdn
           end if
        case default
           print *, 'Unknown test_case: ', trim(test_case)
           stop
        end select
        if (iso_cs == 1) q(iP,ix,iy) = q(irho,ix,iy)*cs0**2

#if NDUST>0
        do idust = 1, ndust
           ivdd = ivdx(idust)
           if (direction == 2) ivdd = ivdy(idust)
           select case (trim(test_case))
           case ('dustybox')
              q(irhod(idust),ix,iy) = eps_dust(idust)
              q(ivdd,ix,iy)         = vd0(idust)
           case ('dustywave')
              q(irhod(idust),ix,iy) = eps_dust(idust)*(1.0d0 + amp*sn)
              q(ivdd,ix,iy)         = cs0*amp*sn
           case ('smoluchowski', 'disk_kernel')  ! n(m,0) = exp(-m): rho_i = m_i exp(-m_i) dm_i
              q(irhod(idust),ix,iy) = max(epsilondust(1,idust), dust_ratio_min*q(irho,ix,iy))   ! floor: empty bins would give 1/0 in ctoprim
           case ('pscal_advect')                 ! dust advected at constant velocity vd0, no drag
              q(irhod(idust),ix,iy) = eps_dust(idust)*(1.0d0 + 0.5d0*sn)
              q(ivdd,ix,iy)         = vd0(idust)
           case ('fragmentation')                ! MRN between scutmin and scut (dust_distribution = 'mrn')
              q(irhod(idust),ix,iy) = epsilondust(1,idust)*q(irho,ix,iy)
           case ('dustyshock')
              if (xx < x_disc) then
                 q(irhod(idust),ix,iy) = eps_shock;  q(ivdd,ix,iy) = v_shock
              else
                 q(irhod(idust),ix,iy) = md/vdn;     q(ivdd,ix,iy) = vdn
              end if
           case default
              print *, 'test_case ', trim(test_case), ' must be compiled with NDUST=0'
              stop
           end select
#if NDUSTPSCAL>0
           ! Passive scalar (first one): profile for pscal_advect, uniform 1 otherwise (must stay exactly 1)
           q(idust_pscal(idust,1),ix,iy) = 1.0d0
           if (trim(test_case) == 'pscal_advect') q(idust_pscal(idust,1),ix,iy) = 1.0d0 + ps_amp*sin(kw*xx + 1.0d0)
#endif
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
       & eps_dust,ts_dust,vd0,vg0,v_shock,eps_shock,K_drag,rho_frag,h_test,r_in_test,direction,ps_amp
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
           else if (trim(test_case) == 'disk_kernel') then       ! t_s proportional to size
              tstop(idust,ix,iy) = ts_dust(1)*sdust(idust)/sdust(1)
           else if (trim(test_case) == 'fragmentation') then     ! Epstein, isothermal (cgs units)
              tstop(idust,ix,iy) = sqrt(pi/8.0d0)*rhograin*sdust(idust)/(q(irho,ix,iy)*cs(ix,iy))
           else
              tstop(idust,ix,iy) = ts_dust(idust)
           end if
        end do
     end do
  end do
  !$omp end parallel do
end subroutine compute_tstop
#endif


! Gas scale height, used by the vertically integrated coagulation kernel (GEOM=2)
function gas_scale_height(ix,iy)
  use parameters
  use commons
  implicit none
  integer  :: ix, iy
  real(dp) :: gas_scale_height
  gas_scale_height = h_test
end function gas_scale_height

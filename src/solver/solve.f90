!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This solves the equations of hydro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solve(verbose,outputing)
  use parameters
  use commons
  use units
  implicit none
  logical :: verbose,outputing
  integer::clock_rate, clock_max,t1,t2,t3,t4,t5,t6,t7,t8,t9, t10,i
  real(dp):: tall
  call system_clock ( t1, clock_rate, clock_max )

  call apply_boundaries !Boundaries are applied here.
  if(force_kick) call kick(0.5d0)
  call ctoprim
  call system_clock ( t2, clock_rate, clock_max )


  if(force_kick) call update_force_setup

  call system_clock ( t3, clock_rate, clock_max )
#if NDUST>0
  ! Re-calc distribution
  call distribution_dust(.false.)
  call compute_tstop  !Re-calc distribution
#endif


  if(charging) then
    if (analytical_charging .eqv. .false.) call charge !Set res_Marchand=True to compute charges AND res

#if NDUST>0
    if(analytical_charging) call analytical_charge
#endif
   endif

#if MHD==1
#if NDUST>0
  if(dust_inertia) then

    if(dusty_nonideal_MHD_no_electron) then
    
        call Hall_factor
        if (hyper_diffusion) call effective_diffusion_coef_induction 
        if (call_electric_field) call electric_field
        if (apply_Lorentz_force) call Lorentz_force


    endif

  endif
#endif
#endif






  call system_clock ( t4, clock_rate, clock_max )

  ! We compute the stability timestep
  call courant


  call system_clock ( t5, clock_rate, clock_max )
  
  ! Predictor step. Variables are estimated at cell interfaces and half dt
  call predictor
  call system_clock ( t6, clock_rate, clock_max )

  ! Flux are computed and added to u_prim
  call add_delta_u
  !call ctoprim
  !call apply_boundaries


  call system_clock ( t7, clock_rate, clock_max )


  ! Source terms are computed and added to u_prim

  call apply_boundaries
  call source_terms
     ! print*, 'dxBy=', dxBy
     ! print*, 'dxBz=', dxBz
    ! print *, 'Ex=', E_x(:)
    ! print *, 'Ey=', E_y(:)
  ! print *, 'Ez=', E_z(:)
    ! print*, 'By=', q(:,iBy)
    ! print*, 'Bz=', q(:,iBz)
     ! print*, 'FL=', FLor_x_d(:,:)
     ! print*, 'FL=', FLor_x_d(:,:)


  if(fargo) call fargo_scheme

  call system_clock ( t8, clock_rate, clock_max )

#if NDUST>0
  ! Dust step (dynamics, growth, charging)
  if(drag)   call dust_drag(1.0d0) ! Second half kick
  if(growth) then
     call ctoprim
     call dust_growth(verbose)
  endif

#if NDUSTPSCAL > 0 
  if(growth_step) call dust_growth_stepinski(verbose) ! Dust growth with Stepinski /!\ dust size is in the first pscal
#endif
#endif
  call system_clock ( t9, clock_rate, clock_max )

  ! if(force_kick) call kick(1.0d0)
  if(force_kick) call kick(0.5d0)

#if TURB>0
!Maybe to put in force_kick

  call compute_rms_velocity

  if(driven_turb) then
    !!Initialization of velocity/acceleration modes to be done within setup
#if NY==1
     !call driven_turbulence !works only in 1D - 1.5D
     if (phase_drift) then
        call kick_phase_drift
        iseed_phase_drift = iseed_phase_drift + 1
     end if 

     call add_driven_turb_kick
#endif
  end if


#endif 
  ! Setup related modifs
  call setup_inloop
  call system_clock ( t10,  clock_rate, clock_max )
  
  t21  = t21   + real ( t2 - t1 )  / real ( clock_rate )
  t32  = t32   + real ( t3 - t2 )  / real ( clock_rate )
  t43  = t43   + real ( t4 - t3 )  / real ( clock_rate )
  t54  = t54   + real ( t5 - t4 )  / real ( clock_rate )
  t65  = t65   + real ( t6 - t5 )  / real ( clock_rate )
  t76  = t76   + real ( t7 - t6 )  / real ( clock_rate )
  t87  = t87   + real ( t8 - t7 )  / real ( clock_rate )
  t98  = t98   + real ( t9 - t8 )  / real ( clock_rate )
  t109 = t109  + real ( t10 - t9 ) / real ( clock_rate )

  if(verbose) then
     write(*,*) "Time spent in each routines: cumulative percentage & real time of current timestep & real time cumulative  "
     tall = t21 + t32 + t43 + t54 + t65 + t76 + t87 + t98 + t109
     write ( *, * ) 'Boundaries + ctoprim      ', t21/tall*100. ,'(%) ,',t21,' seconds'
     write ( *, * ) 'Forces                    ', t32/tall*100. ,'(%) ,',t32,' seconds'
     write ( *, * ) 'Charging + stoping time   ', t43/tall*100. ,'(%) ,',t43,' seconds'
     write ( *, * ) 'Courant                   ', t54/tall*100. ,'(%) ,',t54,' seconds'
     write ( *, * ) 'Predictor step            ', t65/tall*100. ,'(%) ,',t65,' seconds'
     write ( *, * ) 'Corrector step            ', t76/tall*100. ,'(%) ,',t76,' seconds'
     write ( *, * ) 'Source terms              ', t87/tall*100. ,'(%) ,',t87,' seconds'
     write ( *, * ) 'Dust                      ', t98/tall*100. ,'(%) ,',t98,' seconds'
     write ( *, * ) 'Force kick + setup action ', t109/tall*100. ,'(%) ,',t109,' seconds'
     write ( *, * ) 'Running time is ', tall,' seconds'

  endif
  
end subroutine solve



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This routine computes the primitive variables from the conservative ones
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ctoprim
  use parameters
  use commons
  use units
  use OMP_LIB

  implicit none
  integer :: i,idust,ipscal
  real(dp):: ekin,cs_eos,barotrop,emag

  ! Gas related primitive quantities
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,emag,ekin)
  !$OMP DO
  do i = 1,ncells
    q(i,irho)             = max(u_prim(i,irho),smallr)
    q(i,ivx)              = u_prim(i,ivx)/u_prim(i,irho)
    q(i,ivy)              = u_prim(i,ivy)/u_prim(i,irho)
    q(i,ivz)              = u_prim(i,ivz)/u_prim(i,irho)
    
    ekin = half*u_prim(i,irho)*((u_prim(i,ivx)/u_prim(i,irho))**2.0) + half*u_prim(i,irho)*((u_prim(i,ivy)/u_prim(i,irho))**2.0) + half*u_prim(i,irho)*((u_prim(i,ivz)/u_prim(i,irho))**2.0)

    emag=0.0d0
#if MHD==1   
    emag                  = half*(u_prim(i,iBx)**2.0+u_prim(i,iBy)**2.0+u_prim(i,iBz)**2.0)
#endif    
    q(i,iP)               = max((gamma-1.0d0)*(u_prim(i,iP)-ekin-emag),smallp) !TODO : substract magnetic nrj

    if(iso_cs<1)              cs(i) = sqrt(gamma*q(i,iP)/q(i,irho))
    if(non_standard_eos == 1) cs(i) = cs_eos(barotrop(q(i,irho)))
    if(iso_cs==1 .or. non_standard_eos == 1) q(i,iP)  =  u_prim(i,irho)*cs(i)**2
#if NDUST>0
    do idust = 1,ndust
        q(i,irhod(idust)) = u_prim(i,irhod(idust))
        q(i,ivdx(idust))  = u_prim(i,ivdx(idust))/u_prim(i,irhod(idust))
        q(i,ivdy(idust))  = u_prim(i,ivdy(idust))/u_prim(i,irhod(idust))
        q(i,ivdz(idust))  = u_prim(i,ivdz(idust))/u_prim(i,irhod(idust))
#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
        q(i,idust_pscal(idust,ipscal))  = u_prim(i,idust_pscal(idust,ipscal))/u_prim(i,irhod(idust))
        if(growth_step)  sdust(i,idust) = q(i,idust_pscal(idust,ipscal))
    end do
#endif             
     end do
#endif
#if MHD==1
      q(i,iBx)              = u_prim(i,iBx)
      q(i,iBy)              = u_prim(i,iBy)
      q(i,iBz)              = u_prim(i,iBz)
#endif
  end do
  !$OMP END DO
  !$OMP END PARALLEL
#if GRAVITY==1  
  call mtot
#endif

end subroutine ctoprim


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This routine computes the conservative variables from the primitive ones
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine primtoc
  use parameters
  use commons
  use units 
  use OMP_LIB

  implicit none
  integer :: i,idust,ipscal

  do i = 1 ,ncells
     u_prim(i,irho)              = q(i,irho)
     u_prim(i,ivx)               = q(i,irho) * q(i,ivx)
     u_prim(i,ivy)               = q(i,irho) * q(i,ivy)
     u_prim(i,ivz)               = q(i,irho) * q(i,ivz)
     u_prim(i,iP)                = q(i,iP)/(gamma-1.0d0)+half* q(i,irho)*q(i,ivx)**2 + half* q(i,irho)*q(i,ivy)**2 + half* q(i,irho)*q(i,ivz)**2
#if NDUST>0     
     do idust = 1,ndust
        u_prim(i,irhod(idust))  = q(i,irhod(idust))
        u_prim(i,ivdx(idust))   = q(i,irhod(idust)) * q(i,ivdx(idust))
        u_prim(i,ivdy(idust))   = q(i,irhod(idust)) * q(i,ivdy(idust))
        u_prim(i,ivdz(idust))   = q(i,irhod(idust)) * q(i,ivdz(idust))
#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
        u_prim(i,idust_pscal(idust,ipscal))  = q(i,irhod(idust))*q(i,idust_pscal(idust,ipscal))
    end do
#endif
     end do
#endif    
#if MHD==1
      u_prim(i,iBx)              = q(i,iBx)
      u_prim(i,iBy)              = q(i,iBy)
      u_prim(i,iBz)              = q(i,iBz)
      u_prim(i,iP)               = u_prim(i,iP) + half*(q(i,iBx)**2.0+q(i,iBy)**2.0+q(i,iBz)**2.0)
#endif 
  end do


end subroutine primtoc



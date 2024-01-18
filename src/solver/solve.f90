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
  integer::clock_rate, clock_max,t1,t2,t3,t4,t5,t6,t7,t8
  real(dp):: tall
  call system_clock ( t1, clock_rate, clock_max )

  call apply_boundaries !Boundaries are applied here.
  call ctoprim

  if(force_kick) call update_force_setup

#if NDUST>0
  ! Re-calc distribution
  call distribution_dust(.false.)
  call compute_tstop  !Re-calc distribution
#endif
  if(charging)   call charge
 
  ! We compute the stability timestep
  call courant

  call system_clock ( t2, clock_rate, clock_max )
  
  ! Predictor step. Variables are estimated at cell interfaces and half dt
  call predictor
  call system_clock ( t3, clock_rate, clock_max )

  ! Flux are computed and added to u_prim
  call add_delta_u


  ! Source terms are computed and added to u_prim
  call source_terms

  call system_clock ( t4, clock_rate, clock_max )

#if NDUST>0
  ! Dust step (dynamics, growth, charging)
  if(drag)   call dust_drag(1.0d0) ! Second half kick
  if(growth) call dust_growth(verbose)
#endif
  call system_clock ( t5, clock_rate, clock_max )

  if(force_kick) call kick(1.0d0)
  ! Setup related modifs
  call setup_inloop
  call system_clock ( t6,  clock_rate, clock_max )
  
  t21 = t21 + real ( t2 - t1 ) / real ( clock_rate )
  t32 = t32 + real ( t3 - t2 ) / real ( clock_rate )
  t43 = t43 + real ( t4 - t3 ) / real ( clock_rate )
  t54 = t54 + real ( t5 - t4 ) / real ( clock_rate )
  t65 = t65 + real ( t6 - t5 ) / real ( clock_rate )

  if(verbose) then
     write(*,*) "Time spent in each routines: cumulative percentage & real time of current timestep & real time cumulative  "
     tall = t21 + t32 + t43 + t54
     write ( *, * ) 'Courant      ', t21/tall*100. ,'(%) &',real ( t2- t1 ) / real ( clock_rate )*1e6," (mus)", t21,"(s)"
     write ( *, * ) 'Predictor    ', t32/tall*100. ,'(%) &',real ( t3- t2 ) / real ( clock_rate )*1e6," (mus)", t32," (s)"
     write ( *, * ) 'Delta U  + Source    ', t43/tall*100. ,'(%) &',real ( t4- t3 ) / real ( clock_rate )*1e6," (mus)", t43," (s)"
#if NDUST>0     
     write ( *, * ) 'Dust         ', t54/tall*100. ,'(%) &',real ( t5- t4 ) / real ( clock_rate )*1e6," (mus)", t54," (s)"
#endif
     write ( *, * ) 'Others       ', t65/tall*100. ,'(%) &',real ( t6- t5 ) / real ( clock_rate )*1e6," (mus)", t65," (s)"
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
  integer :: i,idust
  real(dp):: ekin,cs_eos,barotrop,emag

  ! Gas related primitive quantities
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust,emag,ekin)
  !$OMP DO
  do i = 1,ncells
    q(i,irho)             = max(u_prim(i,irho),smallr)
    q(i,ivx)              = u_prim(i,ivx)/u_prim(i,irho)

    ekin = half*u_prim(i,irho)*((u_prim(i,ivx)/u_prim(i,irho))**2.0) + half*u_prim(i,irho)*((u_prim(i,ivy)/u_prim(i,irho))**2.0) + half*u_prim(i,irho)*((u_prim(i,ivz)/u_prim(i,irho))**2.0)
    q(i,ivz)              = u_prim(i,ivz)/u_prim(i,irho)
    q(i,ivy)              = u_prim(i,ivy)/u_prim(i,irho)
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
  integer :: i,idust

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



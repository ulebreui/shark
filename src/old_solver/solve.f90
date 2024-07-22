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

  if(force_kick) call kick(half)
#if NDUST>0
  !Dust step  - half kick
  if(drag)   call dust_drag(half)
  call system_clock ( t7, clock_rate, clock_max )
#endif
  call apply_boundaries !Boundaries are applied here.
  call ctoprim

  ! We compute the stability timestep
  call courant

  call system_clock ( t2, clock_rate, clock_max )

  !Predictor step. Variables are estimated at cell interfaces and half dt
  call predictor
  call system_clock ( t3, clock_rate, clock_max )

  !Flux are computed and added to u_prim
  !call add_delta_u_old
  call add_delta_u

  call system_clock ( t4, clock_rate, clock_max )

  call set_u_prim

  call system_clock ( t5, clock_rate, clock_max )
  

  call system_clock ( t6, clock_rate, clock_max )

  
#if NDUST>0
  !Dust step (dynamics, growth, charging)
  if(drag)   call dust_drag(half) ! Second half kick
  if(growth) call dust_growth(verbose)
  call system_clock ( t7, clock_rate, clock_max )
  t7=t6
#endif

  !Source terms are computed and added to u_prim
  call source_terms
  if(force_kick) call kick(half)
  !Setup related modifs
  call setup_inloop
  
  t21=t21+real ( t2- t1 ) / real ( clock_rate )
  t32=t32+real ( t3- t2 ) / real ( clock_rate )
  t43=t43+real ( t4- t3 ) / real ( clock_rate )
  t54=t54+real ( t5- t4 ) / real ( clock_rate )
  t65=t65+real ( t6- t5 ) / real ( clock_rate )
  t76=t76+real ( t7- t6 ) / real ( clock_rate )
  
  if(verbose) then
     write(*,*) "Time spent in each routines: cumulative percentage & real time of current timestep & real time cumulative  "
     tall= t21+t32+t43+t54+t65+t76
     write ( *, * ) 'Courant      ',t21/tall*100. ,'(%) &',real ( t2- t1 ) / real ( clock_rate )*1e6," (mus)", t21,"(s)"
     write ( *, * ) 'Predictor    ',t32/tall*100. ,'(%) &',real ( t3- t2 ) / real ( clock_rate )*1e6," (mus)",t32," (s)"
     write ( *, * ) 'Delta U      ',t43/tall*100. ,'(%) &',real ( t4- t3 ) / real ( clock_rate )*1e6," (mus)",t43," (s)"
     write ( *, * ) 'Source terms ',t54/tall*100. ,'(%) &',real ( t5- t4 ) / real ( clock_rate )*1e6," (mus)",t54," (s)"
     write ( *, * ) 'Set Uold     ',t65/tall*100. ,'(%) &',real ( t6- t5 ) / real ( clock_rate )*1e6," (mus)",t65," (s)"
#if NDUST>0     
     write ( *, * ) 'Dust         ',t76/tall*100. ,'(%) &',real ( t7- t6 ) / real ( clock_rate )*1e6," (mus)",t76," (s)"
#endif
  endif
  
end subroutine solve


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine sets u_prim and does some regularisation for the density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_u_prim
  use parameters
  use commons
  use units
  use OMP_LIB
  implicit none
  integer ::i


   !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i)
  !$OMP DO
  do i = 1,ncells
     u_prim(i,irho) = max(u_prim(i,irho),1d-27/unit_d) !smallr
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
end subroutine set_u_prim

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

  ! Gas related primitive quantities
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do i = 1,ncells
     q(i,irho)             = u_prim(i,irho)
     q(i,ivx)               = u_prim(i,ivx)/u_prim(i,irho)
#if NY==1
     q(i,iP)               = ((gamma-1.0d0)*u_prim(i,iP)-half*u_prim(i,irho)*((u_prim(i,ivx)/u_prim(i,irho))**2.0))
#endif
#if NY>1  
      q(i,iP)              = ((gamma-1.0d0)*u_prim(i,iP)-half*u_prim(i,irho)*((u_prim(i,ivx)/u_prim(i,irho))**2.0+(u_prim(i,ivy)/u_prim(i,irho))**2.0+(u_prim(i,ivz)/u_prim(i,irho))**2.0))
      q(i,ivz)             = u_prim(i,ivz)/u_prim(i,irho)
      q(i,ivy)             = u_prim(i,ivy)/u_prim(i,irho)
#endif
    if(iso_cs<0) cs(i)       = sqrt(gamma*q(i,iP)/q(i,irho))
    if(iso_cs.ge.0) q(i,iP)  =  u_prim(i,irho)*cs(i)**2
#if NDUST>0
    do idust = 1,ndust
        q(i,irhod(idust)) = u_prim(i,irhod(idust))
        q(i,ivdx(idust))   = u_prim(i,ivdx(idust))/u_prim(i,irhod(idust))
#if NY>1       
        q(i,ivdy(idust))  = u_prim(i,ivdy(idust))/u_prim(i,irhod(idust))
        q(i,ivdz(idust))  = u_prim(i,ivdz(idust))/u_prim(i,irhod(idust))
#endif       
     end do
#endif
  end do
  !$OMP END DO
  !$OMP END PARALLEL


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
     u_prim(i,irho)             = q(i,irho)
     u_prim(i,ivx)               = q(i,irho)*q(i,ivx)
#if NY>1
     u_prim(i,ivy)              = q(i,irho)*q(i,ivy)
     u_prim(i,ivz)              = q(i,irho)*q(i,ivz)
#endif 
     u_prim(i,iP)   = q(i,iP)/(gamma-1.0d0)+half* q(i,irho)*q(i,ivx)**2
#if NY>1
     u_prim(i,iP)   =  u_prim(i,iP) + half* q(i,irho)*q(i,ivy)**2
     u_prim(i,iP)   =  u_prim(i,iP) + half* q(i,irho)*q(i,ivz)**2
#endif
#if NDUST>0     
     do idust = 1,ndust
        u_prim(i,irhod(idust)) = q(i,irhod(idust))
        u_prim(i,ivdx(idust))   = q(i,irhod(idust))*q(i,ivdx(idust))
#if NY>1
        u_prim(i,ivdy(idust))  = q(i,irhod(idust))*q(i,ivdy(idust))
        u_prim(i,ivdz(idust))  = q(i,irhod(idust))*q(i,ivdz(idust))
#endif        
     end do
#endif     
  end do


end subroutine primtoc



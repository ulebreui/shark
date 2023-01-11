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
 
  ! We compute the stability timestep
  
  call courant
  call system_clock ( t2, clock_rate, clock_max )

  !Predictor step. Variables are estimated at cell interfaces and half dt
  
  call set_unew
  call predictor
  call system_clock ( t3, clock_rate, clock_max )

  !Flux are computed
  
  call add_delta_u
  call system_clock ( t4, clock_rate, clock_max )

  !Source terms are computed
  
  call source_terms
  call system_clock ( t5, clock_rate, clock_max )
    
  call set_uold
  call system_clock ( t6, clock_rate, clock_max )
  
  
#if NDUST>0
  !Dust step (dynamics, growth, charging)
  call dust(verbose,outputing)
  call system_clock ( t7, clock_rate, clock_max )
#else
  t7=t6
#endif

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
! This routine sets unew and computes the new primitive variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_unew
  use parameters
  use commons
  use units
  implicit none
  call apply_boundaries(1,uold,ncells,nvar) !Boundaries are applied to uold so it must be done as early as that
  !Unew = Uold, Unew is the temporary state vector
  unew = uold
  !We compute the new primitive variables.
  !Necessary for the Riemmann solver and for the predictor step
  call ctoprim

  
end subroutine set_unew

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine sets uold and does some regularisation for the density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_uold
  use parameters
  use commons
  use units
  implicit none
  integer ::i

  !We update the state vector (and impose a lower density to avoid negative values)
  uold = unew
  do i = 1,ncells
     uold(i,irho) = max(uold(i,irho),1d-27/unit_d) !smallr
  end do
  call apply_boundaries(1,uold,ncells,nvar)

  
end subroutine set_uold

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
  real(dp):: barotrop,cs_eos
  real(dp):: T
  ! Gas related primitive quantities
  do i = 1,ncells
     q(i,irho) = uold(i,irho)
     q(i,iv)   = uold(i,iv)/uold(i,irho)
     cs(i)     = cs_eos(barotrop(uold(i,irho)))
  end do

#if NDUST>0
  ! Dust related primitive quantities
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,idust)
  !$OMP DO
  do idust = 1,ndust
     do i  = 1,ncells
        q(i,irhod(idust)) = uold(i,irhod(idust))
        q(i,ivd(idust))   = uold(i,ivd(idust))/uold(i,irhod(idust))
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
#endif
#if GRAVITY==1  
  call mtot
#endif  
#if NDUST>0  
  call distribution_dust(.false.)
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
  implicit none
  integer :: i,idust
  do i = 1 ,ncells
     uold(i,irho) = q(i,irho)
     uold(i,iv)   = q(i,irho)*q(i,iv)
#if NDUST>0     
     do idust = 1,ndust
        uold(i,irhod(idust)) = q(i,irhod(idust))
        uold(i,ivd(idust))   = q(i,irhod(idust))*q(i,ivd(idust))
     end do
#endif     
  end do
end subroutine primtoc



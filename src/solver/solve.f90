!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This solves the equations of hydro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solve(verbose)
   use parameters
   use commons
   use units
   implicit none
   logical :: verbose
   integer::clock_rate, clock_max, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
   real(dp):: tall
   call system_clock(t1, clock_rate, clock_max)

   call apply_boundaries !Boundaries are applied here.

   call system_clock(t2, clock_rate, clock_max)
   call ctoprim

   call system_clock(t3, clock_rate, clock_max)

    if (force_kick) call update_force_setup

#if NDUST>0
   ! Re-calc distribution
   call compute_tstop  !Re-calc distribution
#endif

   call system_clock(t4, clock_rate, clock_max)

   ! We compute the stability timestep
   call courant
   if (force_kick) call kick

   call system_clock(t5, clock_rate, clock_max)
   

   ! Predictor step. Variables are estimated at cell interfaces and half dt
   call predictor
   call system_clock(t6, clock_rate, clock_max)

   ! Flux are computed and added to u_prim
   call add_delta_u
   call system_clock(t7, clock_rate, clock_max)

   ! Source terms are computed and added to u_prim
   call source_terms

   call system_clock(t8, clock_rate, clock_max)

#if NDUST>0
   ! Dust step (dynamics, growth, charging)
   if (drag) call dust_drag ! Second half kick
   if (growth) call dust_growth(verbose)
#if NDUSTPSCAL > 0
   if (growth_step) call dust_growth_stepinski! Dust growth with Stepinski /!\ dust size is in the first pscal
#endif
#endif
   call system_clock(t9, clock_rate, clock_max)

   ! Setup related modifs
   call setup_inloop
   call system_clock(t10, clock_rate, clock_max)

   t21 = t21 + real(t2 - t1)/real(clock_rate)
   t32 = t32 + real(t3 - t2)/real(clock_rate)
   t43 = t43 + real(t4 - t3)/real(clock_rate)
   t54 = t54 + real(t5 - t4)/real(clock_rate)
   t65 = t65 + real(t6 - t5)/real(clock_rate)
   t76 = t76 + real(t7 - t6)/real(clock_rate)
   t87 = t87 + real(t8 - t7)/real(clock_rate)
   t98 = t98 + real(t9 - t8)/real(clock_rate)
   t109 = t109 + real(t10 - t9)/real(clock_rate)

   if (verbose) then
      write (*, *) "Time spent in each routines: cumulative percentage & real time of current timestep & real time cumulative  "
      tall = t21 + t32 + t43 + t54 + t65 + t76 + t87 + t98 + t109
      write (*, *) 'Boundaries    ', t21/tall*100., '(%) ,', t21, ' seconds'
      write (*, *) 'Ctoprim                 ', t32/tall*100., '(%) ,', t32, ' seconds'
      write (*, *) 'Charging + stoping time  + forces ', t43/tall*100., '(%) ,', t43, ' seconds'
      write (*, *) 'Courant                   ', t54/tall*100., '(%) ,', t54, ' seconds'
      write (*, *) 'Predictor step            ', t65/tall*100., '(%) ,', t65, ' seconds'
      write (*, *) 'Corrector step            ', t76/tall*100., '(%) ,', t76, ' seconds'
      write (*, *) 'Source terms              ', t87/tall*100., '(%) ,', t87, ' seconds'
      write (*, *) 'Dust                      ', t98/tall*100., '(%) ,', t98, ' seconds'
      write (*, *) 'Force kick + setup action ', t109/tall*100., '(%) ,', t109, ' seconds'
      write (*, *) 'Running time is ', tall, ' seconds'

   end if

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
   integer :: idust, ipscal, ix,iy
   real(dp):: ekin, cs_eos, barotrop

   ! Gas related primitive quantities
   !$omp parallel do schedule(RUNTIME) default(shared) private(ekin,idust, ipscal, ix,iy)
   do iy = 1, ny_max
      do ix = 1, nx_max
         q(irho,ix,iy) = max(u_prim(irho,ix,iy), smallr)
         q(ivx,ix,iy)  = u_prim(ivx,ix,iy)/u_prim(irho,ix,iy)
         q(ivy,ix,iy)  = u_prim(ivy,ix,iy)/u_prim(irho,ix,iy)
         q(ivz,ix,iy)  = u_prim(ivz,ix,iy)/u_prim(irho,ix,iy)
#if NDUST>0
         do idust = 1, ndust
            q(irhod(idust),ix,iy) = u_prim(irhod(idust),ix,iy)
            q(ivdx(idust),ix,iy)  = u_prim(ivdx(idust),ix,iy)/u_prim(irhod(idust),ix,iy)
            q(ivdy(idust),ix,iy)  = u_prim(ivdy(idust),ix,iy)/u_prim(irhod(idust),ix,iy)
            q(ivdz(idust),ix,iy)  = u_prim(ivdz(idust),ix,iy)/u_prim(irhod(idust),ix,iy)
#if NDUSTPSCAL>0
            do ipscal = 1, ndustpscal
               q(idust_pscal(idust,ipscal),ix,iy) = u_prim(idust_pscal(idust, ipscal),ix,iy)/u_prim(irhod(idust),ix,iy)
            end do
#endif
         end do
#endif
      end do
   end do


   if (non_standard_eos == 1) then
      !$omp parallel do schedule(RUNTIME) default(shared) private(ix,iy)
      do iy = 1, ny_max
         do ix = 1, nx_max
            cs(ix,iy)   = cs_eos(barotrop(q(irho,ix,iy)))
            q(iP,ix,iy) = u_prim(irho,ix,iy)*cs(ix,iy)**2
         end do
      end do

   else if (iso_cs == 1 ) then 
      !$omp parallel do schedule(RUNTIME) default(shared) private(ix,iy)
      do iy = 1, ny_max
         do ix = 1, nx_max
            q(iP,ix,iy) = u_prim(irho,ix,iy)*cs(ix,iy)**2
         end do
      end do

   else
     !$omp parallel do schedule(RUNTIME) default(shared) private(ix,iy,ekin)
      do iy = 1, ny_max
         do ix = 1, nx_max
            ekin        = half*u_prim(irho,ix,iy)*((u_prim(ivx,ix,iy)/u_prim(irho,ix,iy))**2.0) + half*u_prim(irho,ix,iy)*((u_prim(ivy,ix,iy)/u_prim(irho,ix,iy))**2.0) + half*u_prim(irho,ix,iy)*((u_prim(ivz,ix,iy)/u_prim(irho,ix,iy))**2.0)
            q(iP,ix,iy) = max((gamma - 1.0d0)*(u_prim(iP,ix,iy) - ekin), smallp) 
         end do
      end do

   endif

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
   integer :: idust, ix,iy
   integer :: ipscal

   !$omp parallel do default(shared) schedule(RUNTIME) private(idust, ipscal, ix,iy)
   do iy = 1, ny_max
      do ix = 1, nx_max
         u_prim(irho,ix,iy) = q(irho,ix,iy)
         u_prim(ivx,ix,iy)  = q(irho,ix,iy)*q(ivx,ix,iy)
         u_prim(ivy,ix,iy)  = q(irho,ix,iy)*q(ivy,ix,iy)
         u_prim(ivz,ix,iy)  = q(irho,ix,iy)*q(ivz,ix,iy)
         u_prim(iP,ix,iy)   = q(iP,ix,iy)/(gamma-1.0d0)+half* q(irho,ix,iy)*q(ivx,ix,iy)**2 + half* q(irho,ix,iy)*q(ivy,ix,iy)**2 + half* q(irho,ix,iy)*q(ivz,ix,iy)**2
#if NDUST>0
         do idust = 1, ndust
            u_prim(irhod(idust),ix,iy) = q(irhod(idust),ix,iy)
            u_prim(ivdx(idust),ix,iy)  = q(irhod(idust),ix,iy)*q(ivdx(idust),ix,iy)
            u_prim(ivdy(idust),ix,iy)  = q(irhod(idust),ix,iy)*q(ivdy(idust),ix,iy)
            u_prim(ivdz(idust),ix,iy)  = q(irhod(idust),ix,iy)*q(ivdz(idust),ix,iy)
#if NDUSTPSCAL>0
            do ipscal = 1, ndustpscal
               u_prim(idust_pscal(idust, ipscal),ix,iy) = q(irhod(idust),ix,iy)*q(idust_pscal(idust,ipscal),ix,iy)
            end do
#endif
         end do
#endif
      end do
   end do

end subroutine primtoc


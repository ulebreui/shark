! SHARK program developed by Ugo Lebreuilly. The base hydro solver is largely inspired from DUMSPHERINO (Hennebelle 2021)
! It includes gas + dust + dust growth + self-gravity
program shark
 
  use parameters
  use commons
  implicit none
  integer::clock_rate, clock_max,t1, t2

  ! Hello message to the user
  call prompt

  ! Parameter reading
  call read_params

  ! Setup initialisation
  call setup
  

  call system_clock ( t1, clock_rate, clock_max )
  ! Actual time loop
  call time_loop

  call system_clock ( t2, clock_rate, clock_max )

  write ( *, * ) 'Elapsed real time ', real ( t2- t1 ) / real ( clock_rate ) 
  write (*, *) 'Boundaries  ', t21, ' seconds'
  write (*, *) 'Ctoprim                    ', t32, ' seconds'
  write (*, *) 'Charging + stoping time +force   ', t43, ' seconds'
  write (*, *) 'Courant                   ', t54, ' seconds'
  write (*, *) 'Predictor step            ', t65, ' seconds'
  write (*, *) 'Corrector step            ', t76, ' seconds'
  write (*, *) 'Source terms              ', t87, ' seconds'
  write (*, *) 'Dust                      ', t98, ' seconds'
  write (*, *) 'Force kick + setup action ', t109, ' seconds'
  open (20,file=trim('metrics.dat'))   
  write(20,*) real ( t2- t1 ) / real ( clock_rate )  
  write(20,*) t21
  write(20,*) t32
  write(20,*) t43
  write(20,*) t54
  write(20,*) t65
  write(20,*) t76
  write(20,*) t87
  write(20,*) t98
  write(20,*) t109
  close(20)
end program shark
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Main routine of the code, this is the actual
!time loop were the integration is performed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine time_loop
  use parameters
  use commons
  use units
  implicit none
  integer :: icount,iout,istep
  logical :: verbose,continue_sim,outputing
  

  icount = 0
  istep  = 0
  time   = 0.0d0
  ! We dump a first output corresponding to the ICs
  iout=1
  if(nrestart>0)restarting=1

  if(restarting>0) iout = nrestart
  if(nrestart>0) then 
    call restart(nrestart)
  endif
 
  continue_sim = .true.
  outputing    = .false.
  call setup_preloop ! Anything that must be done before the time loop and that is setup dependent
  
  !Actual time loop, continues until continue_sim=.false.
  call flag_continue(continue_sim)

  do while(continue_sim)

     verbose= .false.
     icount = icount+1
     istep  = istep+1
     if(restarting.eq.0) then
        call check_output(icount,iout,outputing,verbose)
     else
        restarting=0
     endif
     
#if NDUST>0
     if(kernel_type>0) verbose=.true.
#endif     
     !Verify and dump output if needed
     if(verbose) then
        print *, "************************************************************************************"
        print *, "************************************************************************************"
        print *, "************************************************************************************"

        print *, 'Outputing the Time step : ', istep,' Output ', iout-1

        print *, "************************************************************************************"
        print *, "************************************************************************************"
        print *, "************************************************************************************"
        print *, 'time is     ', time
        print *, 'timestep is ', dt
        
     endif
     
     !Actual solving of equations
  

     call solve(verbose)


     !New time
     !print*, 'dt=', dt
     time = time +dt
     call flag_continue(continue_sim)

  enddo
  !Last computation of primitive variables before dumping a last output
  call ctoprim
  !Final dump
  call output(iout)
          
  print *, "The simulation is over. Enjoy your results."


end subroutine time_loop

subroutine prompt
  USE OMP_LIB
  implicit none
    write ( *, * ) '-------------------------------------------------------------'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '-------------------------------------------------------------'
    write ( *, * ) '  ssss hh   hh     aaa      rrrrrr   kk   kk                 '
    write ( *, * ) ' sss   hh   hh    aa aa     rr  rrr  kk  kk                  '
    write ( *, * ) '  sss  hhhhhhh   aaaaaaa    rrrrr    kk kk                   '
    write ( *, * ) '   sss hh   hh  aaa    aaa  rr  rr   kk  kk                  '
    write ( *, * ) ' ssss  hh   hh aaa      aaa rr    rr kk   kk   (LITE)        '
    write ( *, * ) '-------------------------------------------------------------'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!000!!!!!!!!!!DDD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!00!!000!!!!!!!!D!!DDD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!000!!!!!!!D!!!!DDD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!!000!!!!!!!!!D!!!!!!D!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!00!!!!!!!!!!!D!!!!!!D!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!000!!!!!!!!!!!!D!!!DDD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!00000000!!!!!!!!DDDD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '-------------------------------------------------------------'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!Ugo Lebreuilly            !!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!CEA SACLAY                !!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!ugo.lebreuilly@cea.fr     !!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '-------------------------------------------------------------'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write ( *, * ) '------------------------------------------------------------'
    write ( *, * ) '                                  /'                                 
    write ( *, * ) '                               ,::|'
    write ( *, * ) '                             /::::|'
    write ( *, * ) '                            ,:::::\                                      _..'
    write ( *, * ) '         ____........-------,..:::::                                   -  /'
    write ( *, * ) ' _.--"""". . . .      .   .  .  .  ""`-._                            -  .;'
    write ( *, * ) '<. - :::::O......  ...   . . .. . .  .  .""--._                  ,- . .;'
    write ( *, * ) ' `-._  ` `":`:`:`::||||:::::::::::::::::.:. .  ""--._ ,/|     ,- .  .;'
    write ( *, * ) '     """_---       // :::.. ````:`:`::::::::::.:.:.:. .`-`._- .   .;'
    write ( *, * ) '         ""--.__      (       \               ` ``:`:``:::: .   .;'
    write ( *, * ) '                "\""--.:-.     `.                             .:'
    write ( *, * ) '                  \  /    `-._   `.""-----.,-..::(--"".\""`.  `:\'
    write ( *, * ) '                   `/         `-._ \          `-:\          `. `:\             '
    write ( *, * ) '                                   ""            "            `-._)            '


#if OPENMP==1
    write (*,*) ' You are using ',OMP_GET_MAX_THREADS(),' threads'
    !write(*,*) ' You are using ',OMP_SCHEDULE,' scheduling'
    !call omp_set_schedule(omp_sched_dynamic, 4)
    write(*,*) ' You are using static scheduling'

    !call omp_set_schedule(omp_sched_static,8)
    call omp_set_schedule(omp_sched_static,8)

#endif    
end subroutine prompt

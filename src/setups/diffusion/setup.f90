subroutine setup
  use parameters
  use commons
  use units
  implicit none


  integer :: i
  real(dp) :: diffusion_timescale
  allocate(By_old(1:ncells))
  allocate(Bz_old(1:ncells))
  allocate(D_test(1:ncells))
  D_test=0.0d0
  Bz_old=0.0d0
  By_old=0.0d0

  box_l = box_l/unit_l


  call allocate_init
  call gridinit(box_l,box_l)
  q=0.0d0
  iso_cs=1

  static = .true.

  do i = 1,ncells

      q(i,irho) = rho_0/unit_d 

      cs(i)=cs_0/unit_v

      q(i,iP)=q(i,irho)*cs(i)**2
  end do

  i=int(ncells/2)

#if MHD==1
  q(i,iBy)=1.0d0 !Set a Dirac at box center
  q(i,iBz)=1.0d0

  if(diffusion_test) then 
    D_test(:)=D_test_value !Diffusion coefficient
    print*, 'D_test=',D_test
  endif

#endif




  diffusion_timescale=box_l**2/D_test_value  !Diffusion_coefficient == length**2/time
  tend = 0.01*diffusion_timescale 
  print*,'tend=',tend


  call apply_boundaries
  call primtoc

  !deallocate(D_test)
end subroutine setup


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
  integer :: ilun
  real(dp):: info1
  
  read(ilun,*) info1
  time=info1
  close(ilun)
end subroutine read_setup_info


subroutine read_setup_params(ilun,nmlfile)
  use parameters
  use commons
  implicit none
  character(len=70):: nmlfile
  integer :: io,ilun
  logical::nml_ok
  namelist/setup_params/D_test_value,diffusion_test,box_l,rho_0,cs_0
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   print *, "Setup namelist reading  !"
   read(13,setup_params,IOSTAT=io)
   rewind(13)
   if (io/=0) then
      write(*,*) 'Invalid line in the setup namelist'
      stop
   end if
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"

   
 end subroutine read_setup_params

 subroutine flag_continue(continue_sim)
   use parameters
   use commons
   use units
   implicit none
   logical :: continue_sim
   integer :: i
   real(dp) :: cfl_diff
   real(dp) :: dxB,A_mean,Fp,Fm,S


   !Here you can add flags to kill the simulation

    dt= tend/100 !for this test, static=.true. --> the code won't go into courant routine

    !user-defined CFL
    cfl_diff=0.5d0
    dt = min(dt,cfl_diff*(box_l/ncells)**2/D_test_value)



#if MHD==1
  if(time==0.0) return

  if(diffusion_test) then 

    !if static=.true., the rest of source_terms won't be called/computed


          Bz_old=q(:,iBz) !Need for buffer variable: to ensure cell i updated with neighboring cells taken at t = n (and not n+1)
          By_old=q(:,iBy)




        do i=1,ncells  
            if(active_cell(i)==1) then


              call hyper_diffusion_induction_eq(D_test(:),Bz_old(:),dx(i,1),dx(i,1),dx(i,1),i,iBz)
              call hyper_diffusion_induction_eq(D_test(:),By_old(:),dx(i,1),dx(i,1),dx(i,1),i,iBy)

              ! dxB=(Bz_old(i+1)-Bz_old(i))/dx(i,1)

              ! A_mean=(D_test(i)+D_test(i+1))/2

              ! Fp=A_mean*dxB

              ! print*, 'dxB=', dxB 
              ! print*, 'A_mean=', A_mean            
              ! print*, 'Fp=', Fp 

              ! dxB=(Bz_old(i)-Bz_old(i-1))/dx(i,1)

              ! A_mean=(D_test(i-1)+D_test(i))/2

              ! Fm=A_mean*dxB

              ! S=(Fp-Fm)/dx(i,1)

              ! u_prim(i,iBz) = u_prim(i,iBz) + dt*S  

              ! print*, 'dxB=', dxB 
              ! print*, 'A_mean=', A_mean            
              ! print*, 'Fm=', Fm

              ! print*, 'S=', S 

              ! print*, 'Bz=', u_prim(i,iBz)           
           





              !call hyper_diffusion_induction_eq(q(i,iBz),D_test(i),D_test(i-1),D_test(i+1),q(i,iBz),q(i-1,iBz),q(i+1,iBz),dx(i,1),dx(i,1),dx(i,1),i)
              !call hyper_diffusion_induction_eq(q(i,iBy),D_test(i),D_test(i-1),D_test(i+1),q(i,iBy),q(i-1,iBy),q(i+1,iBy),dx(i,1),dx(i,1),dx(i,1),i)
            endif
        end do
  endif
#endif

   if(time>=tend) continue_sim=.false.

 end subroutine flag_continue

 subroutine check_output(icount,iout,outputing,verbose)
   use parameters
   use commons
   use units
   implicit none
   integer :: icount,iout
   logical :: outputing,verbose

        !We make an output at a certain frequency rate of for specific values of the density. This can be tuned at will
     if(icount.eq.freq_out) then
        print *, "time =",time,' tend = ', tend
        verbose=.true.
        outputing=.true.
        icount=0
     endif

     if(time.eq.0) outputing=.true.   
     if(outputing)call output(iout)
     if(outputing) iout=iout+1
     if(outputing) print *, "Outputing data "
     if(outputing) print *, "Total mass is", sum(u_prim(:,irho))
     if(outputing) print *, "Total momentum is", sum(u_prim(:,ivx)+u_prim(:,ivy)+u_prim(:,ivz))
     if(outputing) print *, "Total energy is", sum(u_prim(:,iP))

     outputing=.false.

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
   !call output(1)
   return
end subroutine setup_inloop

 subroutine update_force_setup
   use parameters
   use commons
   use units
   implicit none
  
end subroutine update_force_setup



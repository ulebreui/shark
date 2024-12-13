subroutine setup
  use parameters
  use commons
  use units
  implicit none


  integer :: i
  real(dp) :: diffusion_timescale
  real(dp),EXTERNAL :: Barenblatt_Pattle_Dirac_CI
  real(dp),EXTERNAL :: Barenblatt_Pattle_Lebreuilly_CI
  real(dp),EXTERNAL :: Initial_condition



  box_l = box_l/unit_l


  call allocate_init
  call gridinit(box_l,box_l)
  q=0.0d0
  iso_cs=1

  static = .true.

  do i = 1,ncells

      q(i,irho) = 1.0 

      cs(i)=cs_0/unit_v

      q(i,iP)=q(i,irho)*cs(i)**2
  end do


#if MHD==1


  do i=1,ncells

    q(i,iBy) = Initial_condition(position(i,1),x_wid,By_init) !Initial condition, close to a Dirac


  end do
 
#endif

print*,'By',q(:,iBy)


  !diffusion_timescale=box_l**2/(D_test_value * m * MAXVAL(q(:,iBy))**(m-1))  !Diffusion_coefficient == length**2/time
  diffusion_timescale=box_l**2/(2 * MAXVAL(q(:,iBy)))  !for Lebreuilly CI

  tend = 1*diffusion_timescale 
  print*,'tend=',tend


  call apply_boundaries
  call primtoc

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
  namelist/setup_params/D_test_value,diffusion_test,box_l,rho_0,cs_0,m,Cte,t_init,x_wid,By_init
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
   integer :: i,ivar
   real(dp) :: cfl_diff
   real(dp) :: dxB,A_mean,Fp,Fm
  real(dp), dimension(1:nvar) :: S_diff
  real(dp), dimension(1:ncells) :: D_test

  D_test=0.0d0


   do ivar=1,nvar
        S_diff(ivar)=0.0d0
    end do 


   !Here you can add flags to kill the simulation


    dt= tend/1000


   if (time==0.0d0) return


    do i=1,ncells  


        !D_test(i)=D_test_value * m * q(i,iBy)**(m-1) !Diffusion coefficient must be computed in ghost cells to solve diffusion
        D_test(i)= 2 * q(i,iBy) !Diffusion coefficient must be computed in ghost cells to solve diffusion


          !for this test, static=.true. --> the code won't go into courant routine

          !user-defined CFL
          cfl_diff=0.4d0

          dt = min(dt,cfl_diff*(box_l/ncells)**2/D_test(i))


          !print *,'dt diff', cfl_diff*(box_l/ncells)**2/D_test(i)
          !print*,'box/ncells',(box_l/ncells)
          !print*,'dx(i,1)',dx(i,1)


        !endif
    end do

      !print *, 'tend/100000',tend/100000
      !print *,'courant',D_test(:)



#if MHD==1
  !if(time==0.0) return

  if(diffusion_test) then 

    !if static=.true., the rest of source_terms won't be called/computed




        do i=1,ncells  
            if(active_cell(i)==1) then

              !print*, 'D_test',D_test(i)

              call hyper_diffusion_induction_eq_test(S_diff(:),D_test(:),q(:,iBy),dx(i,1),dx(i,1),dx(i,1),i,iBy)
              u_prim(i,iBy) = u_prim(i,iBy) + S_diff(iBy)*dt
              !print*,'B increment',S_diff(iBy)*dt  
             !print*,'dt',dt  




            endif
        end do
  endif
#endif

   if(time>=tend) continue_sim=.false.

 end subroutine flag_continue



subroutine hyper_diffusion_induction_eq_test(S_diff,A,B,Delta_x,Delta_xm,Delta_xp,i,i_variable) !In the form: dt(V)=dx(Adx(B))=dx(F)

   use parameters
   use commons
   use units
implicit none
integer, intent(in) :: i,i_variable
real(dp) :: dxB,A_mean,Fp,Fm,S
real(dp) :: Delta_x,Delta_xm,Delta_xp
!real(dp) :: Ai,Aim1,Aip1,Bi,Bim1,Bip1,Vi

real(dp), dimension(1:ncells),intent(in) :: A
real(dp), dimension(1:ncells),intent(in) :: B
real(dp), dimension(1:nvar), intent(inout) :: S_diff

!allocate(S_U(1:ncells,1:nvar))

dxB=0.0d0
A_mean=0.0d0
Fp=0d0

!print*,'B',B
!print*,'A',A

!dx(F)=F(i+1/2)-F(i-1/2)/Delta_x=Fp-Fm/Delta_x


!First: F(i+1/2)

!compute dxB (i+1/2) at cell interface
dxB=(B(i+1)-B(i))/Delta_xp
!print*,'dxB',dxB
!dxB=(Bip1-Bi)/Delta_xp

!Compute A (i+1/2) mean at cell interface
A_mean=(A(i)+A(i+1))/2
!print*,'A_mean',A_mean

!A_mean=(Ai+Aip1)/2

!Infer F(i+1/2)
Fp=A_mean*dxB
!print*,'Fp',Fp



dxB=0.0d0
A_mean=0.0d0
Fm=0d0

!Second: F(i-1/2)

!compute dxB (i-1/2) at cell interface
dxB=(B(i)-B(i-1))/Delta_xm
!dxB=(Bi-Bim1)/Delta_xm

!Compute A (i-1/2) mean at cell interface
A_mean=(A(i-1)+A(i))/2
!A_mean=(Aim1+Ai)/2

!Infer F(i-1/2)
Fm=A_mean*dxB
!print*,'Fm',Fm

!print*,'(Fp-Fm)/Delta_x',(Fp-Fm)/Delta_x

!Overall source term: 
!print*,'i',i
S_diff(i_variable)=(Fp-Fm)/Delta_x
!print*,'S_diff',S_diff



end subroutine hyper_diffusion_induction_eq_test


real (dp) function Barenblatt_Pattle_Dirac_CI(x,t,m_lin,C)

  use units

  real(dp),intent(in) :: x
  real(dp),intent(in) :: t
  real(dp),intent(in) :: C
  integer,intent(in) :: m_lin
  real(dp) :: alpha
  real(dp) :: beta
  real(dp) :: factor



  !Typical value for non-linearity ): m_lin = 1 --> dt(u) = dx(dx(u)). m_lin = 2 --> dt(u) = dx(2udx(u))
  alpha = 1.0d0 / (1.0d0 + m_lin)
  !print*,'  alpha',  alpha

 

  beta = alpha*(m_lin-1)/(2 * m_lin)
  factor = (C - beta * ((20*x-10) ** 2) / t**(2 * alpha))

  !C: constant arbitrarily fixed

  !Barenblatt-Pattle solution
  Barenblatt_Pattle_Dirac_CI = DMAX1(0.d0,t ** (-alpha) *( factor ) ** (1 / (m_lin - 1)))
  ! print*,'C',C
  ! print*,'m_lin',m_lin
  ! print*,'  alpha',  alpha
  ! print*,'beta',beta
  ! print*,'factor',factor



  !print *,'Pattle',t ** (-alpha) *( factor ) ** (1 / (m_lin - 1))

end function Barenblatt_Pattle_Dirac_CI



real (dp) function Barenblatt_Pattle_Lebreuilly_CI(x,t,x_width,By_initial)

!Equation of the form: dt(By) = 2dx(Bydx(By))

  use units

  real(dp),intent(in) :: x
  real(dp),intent(in) :: t
  real(dp),intent(in) :: x_width !Control the width of the initial condition
  real(dp),intent(in) :: By_initial
  real(dp) :: C


  C = (By_initial*x_width/dsqrt(6.0d0))**(2/3)  
 


  !Barenblatt-Pattle solution
  Barenblatt_Pattle_Lebreuilly_CI = DMAX1(0.d0,(2*t) ** (-1/3) *( C - 1/6 * (x-0.5)**2/(2*t)**(2/3) ))
  ! print*,'C',C


end function Barenblatt_Pattle_Lebreuilly_CI


real (dp) function Initial_condition(x,x_width,By_initial)
!Initial condition for which the above function is solution to.


  use units

  real(dp),intent(in) :: x
  real(dp),intent(in) :: x_width !Control the width of the initial condition
  real(dp),intent(in) :: By_initial



  Initial_condition = DMAX1(0.d0,By_initial*(1 - ((x-0.5)/x_width)**2 ))


end function Initial_condition



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



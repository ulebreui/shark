subroutine Source_terms
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust,ix,iy
  integer :: ixx, iyy,icell,ivar
  real(dp) :: ts,ekin,lap_x_u,lap_y_u,lap_x_v,lap_y_v,cs_eos,barotrop
  real(dp), dimension(:,:)  , allocatable :: S_U
  real(dp), dimension(1:nvar) :: S_diff


  if(static) then
     return 
  endif

  allocate(S_U(1:ncells,1:nvar))
  S_U=0.0d0

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,ix,iy,idust,ts,lap_x_u,lap_y_u,lap_x_v,lap_y_v,ivar,S_diff)
  !$OMP DO

! #if NDUST>0
!   call compute_tstop
! #endif

  do i=1,ncells
    if(active_cell(i)==1) then

#if GEOM==1
      S_U(i,ivx) = S_U(i,ivx)-dt*(- 2.0d0*q(i,irho)*cs(i)**2./radii_c(i))
#endif  
#if GEOM==2
     S_U(i,ivx)  = S_U(i,ivx)+dt*(q(i,irho)*q(i,ivy)**2+q(i,iP))/radii_c(i)
     S_U(i,ivy)  = S_U(i,ivy)-dt*(q(i,irho)*q(i,ivy)*q(i,ivx)/radii_c(i))

#if NDUST>0
    do idust=1,ndust
      S_U(i,ivdx(idust))=S_U(i,ivdx(idust))+dt*(q(i,irhod(idust))*q(i,ivdy(idust))**2/radii_c(i))
      S_U(i,ivdy(idust))=S_U(i,ivdy(idust))-dt*(q(i,irhod(idust))*q(i,ivdy(idust))*q(i,ivdx(idust))/radii_c(i))
    end do
#endif
#endif 

#if GEOM==4 
! vr is still vx but vphi is now vz !
    !No source term to vy (which is vz here)
     S_U(i,ivx)  = S_U(i,ivx)+dt*(q(i,irho)*q(i,ivz)**2+q(i,iP))/radii_c(i)
     S_U(i,ivz)  = S_U(i,ivz)-dt*(q(i,irho)*q(i,ivz)*q(i,ivx)/radii_c(i))

#if NDUST>0
    do idust=1,ndust
      S_U(i,ivdx(idust))=S_U(i,ivdx(idust))+dt*(q(i,irhod(idust))*q(i,ivdz(idust))**2/radii_c(i))
      S_U(i,ivdz(idust))=S_U(i,ivdz(idust))-dt*(q(i,irhod(idust))*q(i,ivdz(idust))*q(i,ivdx(idust))/radii_c(i))
    end do
#endif
#endif 


#if GRAVITY==1     
      S_U(i,ivx)= S_U(i,ivx)-dt*(q(i,irho)*Mc(i)/(radii_c(i)**2.+(l_soft/unit_l)**2.))
#if NDUST>0     
     do idust=1,ndust
         S_U(i,ivdx(idust))= S_U(i,ivdx(idust))-dt*(q(i,irhod(idust))*Mc(i)/(radii_c(i)**2.+(l_soft/unit_l)**2.)) 
      end do
#endif
#endif

#if NDUST>0
      if(dust_back_reaction .and. iso_cs<1) then
        do idust = 1,ndust
          S_U(i,iP) = S_U(i,iP) + dt*(q(i,ivx)*(q(i,ivx)-q(i,ivdx(idust))))/tstop(i,idust) &
          & + dt*(q(i,ivy)*(q(i,ivy)-q(i,ivdy(idust))))/tstop(i,idust)&
          & + dt*(q(i,ivz)*(q(i,ivz)-q(i,ivdz(idust))))/tstop(i,idust)
     
        end do
        endif
#endif
       !Here we add the viscosity source terme
        if(viscosity) then
          ix=ixx(i)
          iy=iyy(i)
#if GEOM==0
          lap_x_u = (q(icell(ix+1,iy),ivx)+q(icell(ix-1,iy),ivx)-2.0d0*q(i,ivx))/dx(i,1)**2
          S_U(i,ivx) = S_U(i,ivx) + q(i,irho)*dt*eta_visc(i)*lap_x_u
          lap_y_u = (q(icell(ix,iy+1),ivx)+q(icell(ix,iy-1),ivx)-2.0d0*q(i,ivx))/dx(i,2)**2
          S_U(i,ivx) = S_U(i,ivx) + q(i,irho)*dt*eta_visc(i)*lap_y_u

          lap_x_v = (q(icell(ix+1,iy),ivy)+q(icell(ix-1,iy),ivy)-2.0d0*q(i,ivy))/dx(i,1)**2
          lap_y_v = (q(icell(ix,iy+1),ivy)+q(icell(ix,iy-1),ivy)-2.0d0*q(i,ivy))/dx(i,2)**2
          S_U(i,ivy) = S_U(i,ivy) + q(i,irho)*dt*eta_visc(i)*(lap_x_v+lap_y_v)
#endif
#if GEOM==2
          lap_x_u = (q(icell(ix+1,iy),ivx)+q(icell(ix-1,iy),ivx)-2.0d0*q(i,ivx))/dx(i,1)**2+(q(icell(ix+1,iy),ivx)+q(icell(ix-1,iy),ivx))/radii_c(i)/dx(i,1)/2.0d0-q(i,ivx)/radii_c(i)/radii_c(i)
          S_U(i,ivx) = S_U(i,ivx) + q(i,irho)*dt*eta_visc(i)*lap_x_u
          lap_x_v = (q(icell(ix+1,iy),ivy)+q(icell(ix-1,iy),ivy)-2.0d0*q(i,ivy))/dx(i,1)**2-q(i,ivy)/radii_c(i)**2+(q(icell(ix+1,iy),ivy)+q(icell(ix-1,iy),ivy))/radii_c(i)/dx(i,1)/2.0d0
          S_U(i,ivy) = S_U(i,ivy) + q(i,irho)*dt*eta_visc(i)*lap_x_v
#endif 
       endif !viscosity if


#if NDUST>0
!Non-ideal MHD hyper diffusion source term in induction equation
#if MHD==1
#if GEOM==0
if (dusty_nonideal_MHD) then


    do ivar=1,nvar
        S_diff(ivar)=0.0d0
    end do


    call effective_diffusion_coef_induction !To compute effective diffusion coeffs

    ! print*, "eta_eff_yy=",eta_eff_yy
    ! print*, "eta_eff_yz=",eta_eff_yz


    !!!Keep variables in right hand part of the equation a time n (By et By are coupled)!!! --> ok since we work with q() and update u_prim

    !Hyper diffusion term for By
    call hyper_diffusion_induction_eq(S_diff(:),eta_eff_yy(:),q(:,iBy),dx(i,1),dx(i,1),dx(i,1),i,iBy) !!Resistivities must be in cm^2/s
    S_U(i,iBy)=S_U(i,iBy)+S_diff(iBy)

    call hyper_diffusion_induction_eq(S_diff(:),eta_eff_yz(:),q(:,iBz),dx(i,1),dx(i,1),dx(i,1),i,iBy)
    S_U(i,iBy)=S_U(i,iBy)+S_diff(iBy)


    !Hyper diffusion term for Bz
    call hyper_diffusion_induction_eq(S_diff(:),eta_eff_zy(:),q(:,iBy),dx(i,1),dx(i,1),dx(i,1),i,iBz)
    S_U(i,iBz)=S_U(i,iBz)+S_diff(iBz)

    call hyper_diffusion_induction_eq(S_diff(:),eta_eff_zz(:),q(:,iBz),dx(i,1),dx(i,1),dx(i,1),i,iBz)
    S_U(i,iBz)=S_U(i,iBz)+S_diff(iBz)


    ! !Lorentz forces
    ! call electric_field !Compute electrical current as well as ions and electrons velocity
    ! do idust=1,ndust
    !     S_U(i,ivdx(idust)) = zd(i,idust)*e_el_stat/clight*(q(i,irhod(idust))/mdust(i,idust))*(E_x(i) + q(i,ivdy(idust))*q(i,iBz) - q(i,ivdz(idust))*q(i,iBy))*dt !Check expression
    !     S_U(i,ivdy(idust)) = zd(i,idust)*e_el_stat/clight*(q(i,irhod(idust))/mdust(i,idust))*(E_y(i) + q(i,ivdz(idust))*q(i,iBx) - q(i,ivdx(idust))*q(i,iBz))*dt
    !     S_U(i,ivdz(idust)) = zd(i,idust)*e_el_stat/clight*(q(i,irhod(idust))/mdust(i,idust))*(E_z(i) + q(i,ivdx(idust))*q(i,iBy) - q(i,ivdy(idust))*q(i,iBx))*dt
    ! end do

    ! !Now the gas: friction with ions and electrons translates into Lorentz forces (because of their neglected inertia: balance between friction and Lorentz force)
    ! S_U(i,ivx) = (e_el_stat/clight*ni(i)*(E_x(i) + v_i_y(i)*q(i,iBz) - v_i_z(i)*q(i,iBy)) - e_el_stat/clight*ne(i)*(E_x(i) + v_e_y(i)*q(i,iBz) - v_e_z(i)*q(i,iBy)))*dt !ions + electrons. Check units
    ! S_U(i,ivy) = (e_el_stat/clight*ni(i)*(E_y(i) + v_i_z(i)*q(i,iBx) - v_i_x(i)*q(i,iBz)) - e_el_stat/clight*ne(i)*(E_y(i) + v_e_z(i)*q(i,iBx) - v_e_x(i)*q(i,iBz)))*dt !ions + electrons. Check units
    ! S_U(i,ivz) = (e_el_stat/clight*ni(i)*(E_z(i) + v_i_x(i)*q(i,iBy) - v_i_y(i)*q(i,iBx)) - e_el_stat/clight*ne(i)*(E_z(i) + v_e_x(i)*q(i,iBy) - v_e_y(i)*q(i,iBx)) )*dt!ions + electrons. Check units



endif
#endif
#endif
#endif

    endif !active cell if
  end do !i loop 
  !$OMP END DO
  !$OMP END PARALLEL


  !Update state vector
  u_prim=u_prim+S_U
  deallocate(S_U)

end subroutine Source_terms

subroutine hyper_diffusion_induction_eq(S_diff,A,B,Delta_x,Delta_xm,Delta_xp,i,ivar) !In the form: dt(V)=dx(Adx(B))=dx(F)

use parameters
use commons
use OMP_LIB
use units
implicit none
integer :: i,ivar
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
!dx(F)=F(i+1/2)-F(i-1/2)/Delta_x=Fp-Fm/Delta_x


!First: F(i+1/2)

!compute dxB (i+1/2) at cell interface
dxB=(B(i+1)-B(i))/Delta_xp
!dxB=(Bip1-Bi)/Delta_xp

!Compute A (i+1/2) mean at cell interface
A_mean=(A(i)+A(i+1))/2
!A_mean=(Ai+Aip1)/2

!Infer F(i+1/2)
Fp=A_mean*dxB



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

!Overall source term: 
S_diff(ivar)=(Fp-Fm)/Delta_x




!If testing diffusion (setup diffusion): Comment previous line and substitute with: u_prim=u_prim+S_U*dt plus remove S_U from input param 



end subroutine hyper_diffusion_induction_eq



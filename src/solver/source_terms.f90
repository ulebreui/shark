subroutine Source_terms
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: i,idust,ix,iy
  integer :: ixx, iyy,icell
  real(dp) :: ts,ekin,lap_x_u,lap_y_u,lap_x_v,lap_y_v,cs_eos,barotrop
  real(dp), dimension(:,:)  , allocatable :: S_U

  if(static) return
  allocate(S_U(1:ncells,1:nvar))
  S_U=0.0d0
! #if NDUST>0
!   call compute_tstop
! #endif
  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED)&
  !$OMP PRIVATE(i,ix,iy,idust,ts,lap_x_u,lap_y_u,lap_x_v,lap_y_v)
  !$OMP DO
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


    call effective_diffusion_coef_induction !To compute effective diffusion coeffs

    !Hyper diffusion term for By
    call hyper_diffusion_induction_eq(S_U(:,iBy),eta_eff_yy(:),q(:,iBy),dx(i,1),dx(i,1),dx(i,1),i) !!Resistivities must be in cm^2/s

    call hyper_diffusion_induction_eq(S_U(:,iBy),eta_eff_yz(:),q(:,iBz),dx(i,1),dx(i,1),dx(i,1),i)


    !Hyper diffusion term for Bz
    call hyper_diffusion_induction_eq(S_U(:,iBz),eta_eff_zy(:),q(:,iBy),dx(i,1),dx(i,1),dx(i,1),i)

    call hyper_diffusion_induction_eq(S_U(:,iBz),eta_eff_zz(:),q(:,iBz),dx(i,1),dx(i,1),dx(i,1),i)

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

subroutine hyper_diffusion_induction_eq(V,A,B,Delta_x,Delta_xm,Delta_xp,i) !In the form: dt(V)=dx(Adx(B))=dx(F)
!Non-ideal MHD hyper diffusion source term in induction equation
use parameters
use commons
use OMP_LIB
use units
implicit none
integer :: i
real(dp) :: dxB,A_mean,Fp,Fm,S
real(dp) :: Delta_x,Delta_xm,Delta_xp
real(dp), dimension(1:ncells) :: V
real(dp), dimension(1:ncells) :: A
real(dp), dimension(1:ncells) :: B



!dx(F)=F(i+1/2)-F(i-1/2)/Delta_x=Fp-Fm/Delta_x


!First: F(i+1/2)

!compute dxB (i+1/2) at cell interface
dxB=(B(i+1)-B(i))/Delta_xp
!Compute A (i+1/2) mean at cell interface
A_mean=(A(i)+A(i+1))/2
!Infer F(i+1/2)
Fp=A_mean*dxB

!Second: F(i-1/2)

!compute dxB (i-1/2) at cell interface
dxB=(B(i)-B(i-1))/Delta_xm
!Compute A (i-1/2) mean at cell interface
A_mean=(A(i-1)+A(i))/2
!Infer F(i-1/2)
Fm=A_mean*dxB

!Overall source term: 
S=fp-fm/Delta_x

V(i)=V(i)+S*dt



end subroutine hyper_diffusion_induction_eq



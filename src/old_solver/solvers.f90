module hydro_solvers
contains


subroutine solver_hllc(qleft,qright,flx,csl,csr,idim)
    use parameters
    use commons

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    real(dp) :: csl,csr
    integer  :: idim,idust

    real(dp) :: ustar, Estarleft,Estarright,Pstar,rhostarleft,rhostarright
    real(dp) :: S_lft,S_rgt,hllc_l,hllc_r,r_o,u_o,P_o,e_o,lambda_llf_d


    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt,P_lft,P_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt,E_lft,E_rgt,lambda_llf_g
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt

    !Primitive variablres
    !Density
    rho_rgt   = qright(irho) !rho
    rho_lft   = qleft(irho)

    !Velocity
    u_rgt   = qright(index_vn(idim)) ! u
    u_lft   = qleft(index_vn(idim))
#if NY>1
    !Transverse velocity
    v_rgt   = qright(index_vt(idim)) ! v
    v_lft   = qleft(index_vt(idim))
    w_rgt   = qright(ivz)! w
    w_lft   = qleft(ivz)
#endif

    P_rgt     = qright(iP)
    P_lft     = qleft(iP)

    !Conservative variables

    mom_u_rgt    = rho_rgt * u_rgt   ! rho u
    mom_u_lft    = rho_lft * u_lft   

    !Energy
    E_rgt     = P_rgt  /(gamma-1.d0)   + half * rho_rgt * u_rgt **2
    E_lft     = P_lft  /(gamma-1.d0)   + half * rho_lft * u_lft **2

#if NY>1
    mom_v_rgt     = rho_rgt  * v_rgt ! rho v
    mom_v_lft     = rho_lft  * v_lft  ! rho v

    E_rgt   = E_rgt   + half * rho_rgt   * v_rgt **2
    E_lft   = E_lft   + half * rho_lft   * v_lft **2
    !Second transverse momentum
    mom_w_rgt      = rho_rgt * w_rgt ! rho w
    mom_w_lft      = rho_lft  * w_lft

    E_rgt    = E_rgt   + half * rho_rgt  * w_rgt   **2 ! kinetic energy of z component
    E_lft    = E_lft   + half * rho_lft  * w_lft  **2

#endif

    !Fluxs
    flx_rho_rgt  = rho_rgt  * u_rgt ! rho u or rho u r if v_r
    flx_rho_lft  = rho_lft * u_lft

    flx_mom_u_rgt = (rho_rgt * u_rgt **2 + P_rgt)  ! rho u u + P
    flx_mom_u_lft = (rho_lft * u_lft **2 + P_lft)  ! rho u u + P  

#if NY>1   
    flx_mom_v_rgt  = rho_rgt  * u_rgt * v_rgt ! rho u v
    flx_mom_v_lft  = rho_lft  * u_lft  * v_lft 

    flx_mom_w_rgt  = rho_rgt  * u_rgt * w_rgt ! rho u w
    flx_mom_w_lft  = rho_lft  * u_lft * w_lft ! rho u w
#endif

    flx_P_rgt = (E_rgt +P_rgt)   * u_rgt! (E+P) v
    flx_P_lft = (E_lft +P_lft)   * u_lft

!HLLC
    S_lft  = min(min(u_lft,u_rgt)-max(csl,csr),0.0d0)
    S_rgt  = max(max(u_lft,u_rgt)+max(csl,csr),0.0d0) 

    ! Compute lagrangian sound speed
    hllc_l = rho_lft *(u_lft-S_lft)
    hllc_r = rho_rgt *(S_rgt-u_rgt)

    !Compute star state
    ustar  = (hllc_r*u_rgt     + hllc_l*u_lft   +  (P_lft-P_rgt))/(hllc_r+hllc_l)
    Pstar  = (hllc_r*P_lft     + hllc_l*P_rgt   +  hllc_l*hllc_r*(u_lft-u_rgt))/(hllc_r+hllc_l)

    ! Left star region variables
    rhostarleft=rho_lft*(S_lft-u_lft)/(S_lft-ustar)
    estarleft  =((S_lft-u_lft)*E_lft-P_lft*u_lft+Pstar*ustar)/(S_lft-ustar)

    ! Right star region variables
    rhostarright=rho_rgt*(S_rgt-u_rgt)/(S_rgt-ustar)
    estarright  =((S_rgt-u_rgt)*E_rgt-P_rgt*u_rgt+Pstar*ustar)/(S_rgt-ustar)

    ! Sample the solution at x/t=0
    if(S_lft>0.0d0)then
          r_o=rho_lft
          u_o=u_lft
          P_o=P_lft 
          e_o=E_lft
        else if(ustar>0.0d0)then
          r_o=rhostarleft
          u_o=ustar
          P_o=Pstar
          e_o=estarleft
        else if (S_rgt>0d0)then
          r_o=rhostarright
          u_o=ustar
          P_o=Pstar
          e_o=estarright
        else
          r_o=rho_rgt
          u_o=u_rgt
          P_o=P_rgt
          e_o=E_rgt
    end if
    flx(irho)            = r_o*u_o
    flx(index_vn(idim))  = r_o*u_o*u_o+P_o
    flx(iP)              = (e_o+P_o)*u_o
#if NY>1
    if(ustar>0.0d0) then
        flx(index_vt(idim))  = r_o*u_o*v_lft
        flx(ivz)             = r_o*u_o*w_lft
    else
        flx(index_vt(idim))  = r_o*u_o*v_rgt
        flx(ivz)             = r_o*u_o*w_rgt
    endif
#endif       


end subroutine solver_hllc

#if NDUST>0
subroutine solver_dust(qleft,qright,flx,idim)
    use parameters
    use commons

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    real(dp) :: csl,csr
    integer  :: idim,idust,i_u,i_v,i_rho,i_w

    real(dp) :: S_lft,S_rgt,lambda_llf_d


    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt

    do idust=1,ndust

        i_rho= irhod(idust)
        i_u  = index_vdn(idust,idim)
        i_v  = index_vdt(idust,idim)
        i_w  = ivdz(idust)
        !print *, idust, i_rho,i_n,i_t,i_z
        !Dust density 
        rho_rgt  = qright(i_rho)
        rho_lft  = qleft(i_rho)
        !Dust momentum
        u_rgt   = qright(i_u)
        u_lft   = qleft(i_u)
#if NY>1
          !Dust transverse momentum
        v_rgt   = qright(i_v)
        v_lft   = qleft(i_v)
          !Dust second transverse momentum
        w_rgt     = qright(i_w)
        w_lft     = qleft(i_w)
#endif
        mom_u_rgt     =  rho_rgt * u_rgt
        mom_u_lft     =  rho_lft * u_lft
#if NY>1
        mom_v_rgt    =  rho_rgt * v_rgt
        mom_v_lft    =  rho_lft * v_lft
           
        mom_w_rgt    =  rho_rgt  * w_rgt
        mom_w_lft    =  rho_lft  * w_lft
#endif

        flx_rho_rgt   = rho_rgt * u_rgt
        flx_rho_lft   = rho_lft  * u_lft
 
        flx_mom_u_rgt  = rho_rgt * u_rgt**2
        flx_mom_u_lft  = rho_lft * u_lft**2
#if NY>1

        flx_mom_v_rgt = rho_rgt * u_rgt  * v_rgt
        flx_mom_v_lft = rho_lft * u_lft * v_lft
           
        flx_mom_w_rgt = rho_rgt * u_rgt  * w_rgt
        flx_mom_w_lft = rho_lft * u_lft  * w_lft
#endif

flx(i_rho)  =  0.d0
flx(i_u)    =  0.d0
flx(i_v)    =  0.d0
flx(i_w)    =  0.d0

if(u_rgt>0.0d0 .and. u_lft>0.0d0) then
    flx(i_rho)  =  flx_rho_lft 
    flx(i_u)    =  flx_mom_u_lft 
    flx(i_v)    =  flx_mom_v_lft     
    flx(i_w)    =  flx_mom_w_lft 
else if (u_lft<0.0d0 .and. u_rgt<0.0d0) then
    flx(i_rho)  =  flx_rho_rgt
    flx(i_u)    =  flx_mom_u_rgt
    flx(i_v)    =  flx_mom_v_rgt    
    flx(i_w)    =  flx_mom_w_rgt 
else if (u_lft<0.0d0 .and. u_rgt>0.0d0) then
    flx(i_rho)  =  0.d0
    flx(i_u)    =  0.d0
    flx(i_v)    =  0.d0
    flx(i_w)    =  0.d0
else if (u_lft>0.0d0 .and. u_rgt<0.0d0) then
    flx(i_rho)  =  flx_rho_lft   + flx_rho_rgt
    flx(i_u)    =  flx_mom_u_lft + flx_mom_u_rgt
    flx(i_v)    =  flx_mom_v_lft + flx_mom_v_rgt 
    flx(i_w)    =  flx_mom_w_lft + flx_mom_w_rgt 
endif

end do

end subroutine solver_dust

#endif

end module hydro_solvers
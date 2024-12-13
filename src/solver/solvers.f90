module hydro_solvers
contains

subroutine solver_llf(qleft,qright,flx,csl,csr,idim,i)
    use parameters
    use commons

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    real(dp) :: csl,csr
    integer  :: idim,idust,i

    real(dp) :: ustar, Estarleft,Estarright,Pstar,rhostarleft,rhostarright
    real(dp) :: S_lft,S_rgt,hllc_l,hllc_r,r_o,u_o,P_o,e_o,lambda_llf_d


    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt,P_lft,P_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt,E_lft,E_rgt,lambda_llf_g
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt

!For now the only solver that works for MHD==1 with the gas coupled to B (hll to be modified)
!Without dust, B is coupled to the gas
#if MHD==1 
#if NDUST==0 
    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,P_mag_lft,P_mag_rgt,mag_tension_y_lft,mag_tension_y_rgt,mag_tension_z_lft,mag_tension_z_rgt,mag_tension_x_lft,mag_tension_x_rgt 


    real(dp) ::flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt,magnetosonic_fast_rgt,magnetosonic_fast_lft

    Bx_lft   = qleft(iBx)
    Bx_rgt   = qright(iBx)
    By_lft   = qleft(iBy)
    By_rgt   = qright(iBy)
    Bz_lft   = qleft(iBz)
    Bz_rgt   = qright(iBz)

    P_mag_lft   = (Bx_lft**2+Bz_lft**2+By_lft**2)/2.0d0 
    P_mag_rgt   = (Bx_rgt**2+Bz_rgt**2+By_rgt**2)/2.0d0
   
    mag_tension_x_lft = -Bx_lft*Bx_lft
    mag_tension_x_rgt = -Bx_rgt*Bx_rgt
    mag_tension_y_lft = -By_lft*Bx_lft
    mag_tension_y_rgt = -By_rgt*Bx_rgt
    mag_tension_z_lft = -Bz_lft*Bx_lft
    mag_tension_z_rgt = -Bz_rgt*Bx_rgt


#endif
#endif


    !Primitive variables
    !Density
    rho_rgt   = qright(irho) !rho
    rho_lft   = qleft(irho)

    !Velocity
    u_rgt   = qright(index_vn(idim)) ! u
    u_lft   = qleft(index_vn(idim))

    !Transverse velocity
    v_rgt   = qright(index_vt(idim)) ! v
    v_lft   = qleft(index_vt(idim))

    w_rgt   = qright(ivz)! w
    w_lft   = qleft(ivz)


    P_rgt       = qright(iP)
    P_lft       = qleft(iP)

    !Conservative variables

    mom_u_rgt    = rho_rgt * u_rgt   ! rho u
    mom_u_lft    = rho_lft * u_lft   


    !Energy
    E_rgt     = P_rgt  /(gamma-1.d0)   + half * rho_rgt * u_rgt**2
    E_lft     = P_lft  /(gamma-1.d0)   + half * rho_lft * u_lft**2

    mom_v_rgt     = rho_rgt  * v_rgt ! rho v
    mom_v_lft     = rho_lft  * v_lft  ! rho v

    E_rgt   = E_rgt   + half * rho_rgt   * v_rgt **2
    E_lft   = E_lft   + half * rho_lft   * v_lft **2
    !Second transverse momentum
    mom_w_rgt      = rho_rgt * w_rgt ! rho w
    mom_w_lft      = rho_lft * w_lft

    E_rgt    = E_rgt   + half * rho_rgt  * w_rgt  **2 ! kinetic energy of z component
    E_lft    = E_lft   + half * rho_lft  * w_lft  **2



    !Fluxs
    flx_rho_rgt   = rho_rgt  * u_rgt ! rho u or rho u r if v_r
    flx_rho_lft   = rho_lft  * u_lft

    flx_mom_u_rgt = (rho_rgt * u_rgt **2 + P_rgt)  ! rho u u + P
    flx_mom_u_lft = (rho_lft * u_lft **2 + P_lft)  ! rho u u + P  

    flx_mom_v_rgt  = rho_rgt  * u_rgt  * v_rgt ! rho u v
    flx_mom_v_lft  = rho_lft  * u_lft  * v_lft 

    flx_mom_w_rgt  = rho_rgt  * u_rgt * w_rgt ! rho u w
    flx_mom_w_lft  = rho_lft  * u_lft * w_lft ! rho u w


    flx_P_rgt = (E_rgt +P_rgt)   * u_rgt! (E+P) v
    flx_P_lft = (E_lft +P_lft)   * u_lft



    lambda_llf_g         = max(abs(u_lft)+csl,abs(u_rgt)+csr)


#if MHD==1 
#if NDUST==0
!If no dust, field lines frozen to the gas --> magnetic conservative terms to account for
!If dust present, the gas is considered neutral and the code won't go into this part of the routine
    flx_mom_u_rgt  = flx_mom_u_rgt + P_mag_rgt + mag_tension_x_rgt
    flx_mom_u_lft  = flx_mom_u_lft + P_mag_lft + mag_tension_x_lft

    flx_mom_v_rgt = flx_mom_v_rgt + mag_tension_y_rgt
    flx_mom_v_lft = flx_mom_v_lft + mag_tension_y_lft
       
    flx_mom_w_rgt = flx_mom_w_rgt + mag_tension_z_rgt
    flx_mom_w_lft = flx_mom_w_lft + mag_tension_z_lft

    E_rgt = E_rgt + half*(Bx_rgt**2+By_rgt**2+Bz_rgt**2) ! E = epsilon + Kinetic + magnetic
    E_lft = E_lft + half*(Bx_lft**2+By_lft**2+Bz_lft**2)

    flx_P_rgt = (E_rgt + P_rgt + P_mag_rgt)   * u_rgt + Bx_rgt*(Bx_rgt*u_rgt+By_rgt*v_rgt+Bz_rgt*w_rgt) ! (E+P+Pmag) v + B(B.v)
    flx_P_lft = (E_lft + P_lft + P_mag_lft)   * u_lft + Bx_lft*(Bx_lft*u_lft+By_lft*v_lft+Bz_lft*w_lft)

    magnetosonic_fast_rgt = dsqrt(half*(csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt + dsqrt((csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt)**2-4*csr**2*Bx_rgt**2/rho_rgt))) 
    magnetosonic_fast_lft = dsqrt(half*(csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft + dsqrt((csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft)**2-4*csl**2*Bx_lft**2/rho_lft))) 

    lambda_llf_g = max(abs(u_lft)+magnetosonic_fast_lft,abs(u_rgt)+magnetosonic_fast_rgt)
#endif
#endif


    flx(irho)            = half  * (flx_rho_lft   + flx_rho_rgt)    - half*lambda_llf_g * (rho_rgt   - rho_lft)
    flx(index_vn(idim))  = half  * (flx_mom_u_lft + flx_mom_u_rgt)  - half*lambda_llf_g * (mom_u_rgt - mom_u_lft)
    flx(index_vt(idim))  = half  * (flx_mom_v_lft + flx_mom_v_rgt)  - half*lambda_llf_g * (mom_v_rgt - mom_v_lft)
    flx(ivz)             = half  * (flx_mom_w_lft + flx_mom_w_rgt)  - half*lambda_llf_g * (mom_w_rgt - mom_w_lft)

    flx(iP)              = half  * (flx_P_lft+flx_P_rgt)-half*lambda_llf_g* (E_rgt-E_lft)         

end subroutine solver_llf

subroutine solver_hll(qleft,qright,flx,csl,csr,idim,i)
    use parameters
    use commons

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    real(dp) :: csl,csr
    integer  :: idim,idust,i

    real(dp) :: ustar, Estarleft,Estarright,Pstar,rhostarleft,rhostarright
    real(dp) :: S_lft,S_rgt,hllc_l,hllc_r,r_o,u_o,P_o,e_o,lambda_llf_d


    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt,P_lft,P_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt,E_lft,E_rgt,lambda_llf_g
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt

    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,P_mag_lft,P_mag_rgt,mag_tension_y_lft,mag_tension_y_rgt,mag_tension_z_lft,mag_tension_z_rgt,mag_tension_x_lft,mag_tension_x_rgt 
    real(dp) :: magnetosonic_fast_rgt,magnetosonic_fast_lft
    real(dp) :: flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt

#if MHD==1 
#if NDUST==0 
    !Without dust, B is coupled to the gas



        Bx_lft   = qleft(iBx)
        Bx_rgt   = qright(iBx)
        By_lft   = qleft(iBy)
        By_rgt   = qright(iBy)
        Bz_lft   = qleft(iBz)
        Bz_rgt   = qright(iBz)

        P_mag_lft   = (Bx_lft**2+Bz_lft**2+By_lft**2)/2.0d0 
        P_mag_rgt   = (Bx_rgt**2+Bz_rgt**2+By_rgt**2)/2.0d0
       
        mag_tension_x_lft = -Bx_lft*Bx_lft
        mag_tension_x_rgt = -Bx_rgt*Bx_rgt
        mag_tension_y_lft = -By_lft*Bx_lft
        mag_tension_y_rgt = -By_rgt*Bx_rgt
        mag_tension_z_lft = -Bz_lft*Bx_lft
        mag_tension_z_rgt = -Bz_rgt*Bx_rgt
#endif
#endif

    
    !Primitive variables
    !Density
    rho_rgt   = qright(irho) !rho
    rho_lft   = qleft(irho)

    !Velocity
    u_rgt   = qright(index_vn(idim)) ! u
    u_lft   = qleft(index_vn(idim))

    !Transverse velocity
    v_rgt   = qright(index_vt(idim)) ! v
    v_lft   = qleft(index_vt(idim))
    w_rgt   = qright(ivz)! w
    w_lft   = qleft(ivz)


    P_rgt     = qright(iP)
    P_lft     = qleft(iP)

    !Conservative variables

    mom_u_rgt    = rho_rgt * u_rgt   ! rho u
    mom_u_lft    = rho_lft * u_lft   

    !Energy
    E_rgt     = P_rgt  /(gamma-1.d0)   + half * rho_rgt * u_rgt **2
    E_lft     = P_lft  /(gamma-1.d0)   + half * rho_lft * u_lft **2

    mom_v_rgt     = rho_rgt  * v_rgt ! rho v
    mom_v_lft     = rho_lft  * v_lft  ! rho v

    E_rgt   = E_rgt   + half * rho_rgt   * v_rgt **2
    E_lft   = E_lft   + half * rho_lft   * v_lft **2
    !Second transverse momentum
    mom_w_rgt      = rho_rgt * w_rgt ! rho w
    mom_w_lft      = rho_lft  * w_lft

    E_rgt    = E_rgt   + half * rho_rgt  * w_rgt   **2 ! kinetic energy of z component
    E_lft    = E_lft   + half * rho_lft  * w_lft  **2


    !Fluxs
    flx_rho_rgt  = rho_rgt  * u_rgt ! rho u or rho u r if v_r
    flx_rho_lft  = rho_lft * u_lft

    flx_mom_u_rgt = (rho_rgt * u_rgt **2 + P_rgt)  ! rho u u + P
    flx_mom_u_lft = (rho_lft * u_lft **2 + P_lft)  ! rho u u + P  

    flx_mom_v_rgt  = rho_rgt  * u_rgt * v_rgt ! rho u v
    flx_mom_v_lft  = rho_lft  * u_lft  * v_lft 

    flx_mom_w_rgt  = rho_rgt  * u_rgt * w_rgt ! rho u w
    flx_mom_w_lft  = rho_lft  * u_lft * w_lft ! rho u w

    flx_P_rgt = (E_rgt +P_rgt)   * u_rgt! (E+P) v
    flx_P_lft = (E_lft +P_lft)   * u_lft

#if MHD==1 
#if NDUST==0
        flx_mom_u_rgt  = flx_mom_u_rgt + P_mag_rgt + mag_tension_x_rgt
        flx_mom_u_lft  = flx_mom_u_lft + P_mag_lft + mag_tension_x_lft

        flx_mom_v_rgt = flx_mom_v_rgt + mag_tension_y_rgt
        flx_mom_v_lft = flx_mom_v_lft + mag_tension_y_lft
           
        flx_mom_w_rgt = flx_mom_w_rgt + mag_tension_z_rgt
        flx_mom_w_lft = flx_mom_w_lft + mag_tension_z_lft

        E_rgt = E_rgt + half*(Bx_rgt**2+By_rgt**2+Bz_rgt**2) ! E = epsilon + Kinetic + magnetic
        E_lft = E_lft + half*(Bx_lft**2+By_lft**2+Bz_lft**2)

        flx_P_rgt = (E_rgt + P_rgt + P_mag_rgt)   * u_rgt + Bx_rgt*(Bx_rgt*u_rgt+By_rgt*v_rgt+Bz_rgt*w_rgt) ! (E+P+Pmag) v + B(B.v)
        flx_P_lft = (E_lft + P_lft + P_mag_lft)   * u_lft + Bx_lft*(Bx_lft*u_lft+By_lft*v_lft+Bz_lft*w_lft)

#endif
#endif


    S_lft  = min(min(u_lft,u_rgt) -max(csl,csr),0.0d0)
    S_rgt  = max(max(u_lft,u_rgt) +max(csl,csr),0.0d0) 

#if MHD==1
#if NDUST==0
        magnetosonic_fast_rgt = dsqrt(half*(csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt + dsqrt((csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt)**2-4*csr**2*Bx_rgt**2/rho_rgt))) 
        magnetosonic_fast_lft = dsqrt(half*(csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft + dsqrt((csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft)**2-4*csl**2*Bx_lft**2/rho_lft))) 

        S_lft  = min(min(u_lft,u_rgt) -max(magnetosonic_fast_lft,magnetosonic_fast_rgt),0.0d0)
        S_rgt  = max(max(u_lft,u_rgt) +max(magnetosonic_fast_lft,magnetosonic_fast_rgt),0.0d0) 
#endif
#endif

    flx(irho)            = (S_rgt*flx_rho_lft  -S_lft*flx_rho_rgt  + S_rgt*S_lft*(rho_rgt-rho_lft))      / (S_rgt-S_lft)
    flx(index_vn(idim))  = (S_rgt*flx_mom_u_lft-S_lft*flx_mom_u_rgt+ S_rgt*S_lft*(mom_u_rgt-mom_u_lft))  / (S_rgt-S_lft)
    flx(index_vt(idim))  = (S_rgt*flx_mom_v_lft-S_lft*flx_mom_v_rgt+ S_rgt*S_lft*(mom_v_rgt-mom_v_lft))  / (S_rgt-S_lft)
    flx(ivz)             = (S_rgt*flx_mom_w_lft-S_lft*flx_mom_w_rgt+ S_rgt*S_lft*(mom_w_rgt-mom_w_lft))  / (S_rgt-S_lft)

    flx(iP)  = (S_rgt*flx_P_lft-S_lft*flx_P_rgt+S_rgt*S_lft*(E_rgt-E_lft))/(S_rgt-S_lft)           


end subroutine solver_hll

subroutine solver_hllc(qleft,qright,flx,csl,csr,idim,i)
    use parameters
    use commons

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    real(dp) :: csl,csr
    integer  :: idim,idust,i

    real(dp) :: ustar, Estarleft,Estarright,Pstar,rhostarleft,rhostarright
    real(dp) :: S_lft,S_rgt,hllc_l,hllc_r,r_o,u_o,P_o,e_o,lambda_llf_d


    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt,P_lft,P_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt,E_lft,E_rgt,lambda_llf_g
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt
    real(dp) :: magnetosonic_fast_rgt,magnetosonic_fast_lft


#if MHD==1
#if NDUST==0
!Without dust, B is coupled to the gas !TODO modify correspondingly the solver for MHD==1

    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,P_mag_lft,P_mag_rgt,mag_tension_y_lft,mag_tension_y_rgt,mag_tension_z_lft,mag_tension_z_rgt,mag_tension_x_lft,mag_tension_x_rgt 


    real(dp) ::flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt

    Bx_lft   = qleft(iBx)
    Bx_rgt   = qright(iBx)
    By_lft   = qleft(iBy)
    By_rgt   = qright(iBy)
    Bz_lft   = qleft(iBz)
    Bz_rgt   = qright(iBz)

    P_mag_lft   = (Bx_lft**2+Bz_lft**2+By_lft**2)/2.0d0 
    P_mag_rgt   = (Bx_rgt**2+Bz_rgt**2+By_rgt**2)/2.0d0
   
    mag_tension_x_lft = -Bx_lft*Bx_lft
    mag_tension_x_rgt = -Bx_rgt*Bx_rgt
    mag_tension_y_lft = -By_lft*Bx_lft
    mag_tension_y_rgt = -By_rgt*Bx_rgt
    mag_tension_z_lft = -Bz_lft*Bx_lft
    mag_tension_z_rgt = -Bz_rgt*Bx_rgt
#endif
#endif


    !Primitive variables
    !Density
    rho_rgt   = qright(irho) !rho
    rho_lft   = qleft(irho)

    !Velocity
    u_rgt   = qright(index_vn(idim)) ! u
    u_lft   = qleft(index_vn(idim))
    !Transverse velocity
    v_rgt   = qright(index_vt(idim)) ! v
    v_lft   = qleft(index_vt(idim))
    w_rgt   = qright(ivz)! w
    w_lft   = qleft(ivz)


    P_rgt     = qright(iP)
    P_lft     = qleft(iP)

    !Conservative variables

    mom_u_rgt    = rho_rgt * u_rgt   ! rho u
    mom_u_lft    = rho_lft * u_lft   

    mom_v_rgt     = rho_rgt  * v_rgt  ! rho v
    mom_v_lft     = rho_lft  * v_lft  ! rho v


    !Second transverse momentum
    mom_w_rgt      = rho_rgt  * w_rgt ! rho w
    mom_w_lft      = rho_lft  * w_lft

    !Energy
    E_rgt     = P_rgt  /(gamma-1.d0)   + half * rho_rgt * u_rgt **2 + half * rho_rgt   * v_rgt **2 + half * rho_rgt  * w_rgt   **2
    E_lft     = P_lft  /(gamma-1.d0)   + half * rho_lft * u_lft **2 + half * rho_lft   * v_lft **2 + half * rho_lft  * w_lft  **2


    !Fluxs
    flx_rho_rgt  = rho_rgt  * u_rgt ! rho u or rho u r if v_r
    flx_rho_lft  = rho_lft  * u_lft

    flx_mom_u_rgt = (rho_rgt * u_rgt **2 + P_rgt)  ! rho u u + P
    flx_mom_u_lft = (rho_lft * u_lft **2 + P_lft)  ! rho u u + P  

    flx_mom_v_rgt  = rho_rgt  * u_rgt * v_rgt ! rho u v
    flx_mom_v_lft  = rho_lft  * u_lft  * v_lft 

    flx_mom_w_rgt  = rho_rgt  * u_rgt * w_rgt ! rho u w
    flx_mom_w_lft  = rho_lft  * u_lft * w_lft ! rho u w


    flx_P_rgt = (E_rgt +P_rgt)   * u_rgt! (E+P) v
    flx_P_lft = (E_lft +P_lft)   * u_lft



! #if MHD==1 
! #if NDUST==0
!         flx_mom_u_rgt  = flx_mom_u_rgt + P_mag_rgt + mag_tension_x_rgt
!         flx_mom_u_lft  = flx_mom_u_lft + P_mag_lft + mag_tension_x_lft

!         flx_mom_v_rgt = flx_mom_v_rgt + mag_tension_y_rgt
!         flx_mom_v_lft = flx_mom_v_lft + mag_tension_y_lft
           
!         flx_mom_w_rgt = flx_mom_w_rgt + mag_tension_z_rgt
!         flx_mom_w_lft = flx_mom_w_lft + mag_tension_z_lft

!         E_rgt = E_rgt + half*(Bx_rgt**2+By_rgt**2+Bz_rgt**2) ! E = epsilon + Kinetic + magnetic
!         E_lft = E_lft + half*(Bx_lft**2+By_lft**2+Bz_lft**2)

!         flx_P_rgt = (E_rgt + P_rgt + P_mag_rgt)   * u_rgt + Bx_rgt*(Bx_rgt*u_rgt+By_rgt*v_rgt+Bz_rgt*w_rgt) ! (E+P+Pmag) v + B(B.v)
!         flx_P_lft = (E_lft + P_lft + P_mag_lft)   * u_lft + Bx_lft*(Bx_lft*u_lft+By_lft*v_lft+Bz_lft*w_lft)

! #endif
! #endif

    !HLLC
    S_lft  = min(min(u_lft,u_rgt)-max(csl,csr),0.0d0)
    S_rgt  = max(max(u_lft,u_rgt)+max(csl,csr),0.0d0) 

! #if MHD==1
! #if NDUST==0 
! !Does hll_c work with a magnetic field? --> nope
!         magnetosonic_fast_rgt = dsqrt(half*(csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt + dsqrt((csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt)**2-4*csr**2*Bx_rgt**2/rho_rgt))) 
!         magnetosonic_fast_lft = dsqrt(half*(csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft + dsqrt((csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft)**2-4*csl**2*Bx_lft**2/rho_lft))) 

!         S_lft  = min(min(u_lft,u_rgt) -max(magnetosonic_fast_lft,magnetosonic_fast_rgt),0.0d0)
!         S_rgt  = max(max(u_lft,u_rgt) +max(magnetosonic_fast_lft,magnetosonic_fast_rgt),0.0d0) 
! #endif
! #endif

    ! Compute lagrangian sound speed
    hllc_l = rho_lft *(u_lft-S_lft)
    hllc_r = rho_rgt *(S_rgt-u_rgt)

    ! Compute star state
    ustar  = (hllc_r*u_rgt     + hllc_l*u_lft   +  (P_lft-P_rgt))/(hllc_r+hllc_l)
    Pstar  = (hllc_r*P_lft     + hllc_l*P_rgt   +  hllc_l*hllc_r*(u_lft-u_rgt))/(hllc_r+hllc_l)

    ! Left star region variables
    rhostarleft = rho_lft*(S_lft-u_lft)/(S_lft-ustar)
    estarleft   = ((S_lft-u_lft)*E_lft-P_lft*u_lft+Pstar*ustar)/(S_lft-ustar)

    ! Right star region variables
    rhostarright = rho_rgt*(S_rgt-u_rgt)/(S_rgt-ustar)
    estarright   = ((S_rgt-u_rgt)*E_rgt-P_rgt*u_rgt+Pstar*ustar)/(S_rgt-ustar)

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

    flx(irho)                = r_o*u_o
    flx(index_vn(idim))      = r_o*u_o*u_o+P_o
    flx(iP)                  = (e_o+P_o)*u_o
    if(ustar>0.0d0) then
        flx(index_vt(idim))  = r_o*u_o*v_lft
        flx(ivz)             = r_o*u_o*w_lft
    else
        flx(index_vt(idim))  = r_o*u_o*v_rgt
        flx(ivz)             = r_o*u_o*w_rgt
    endif


end subroutine solver_hllc


#if NDUST>0
#if SOLVERDUST==0

subroutine solver_dust_Huang_Bai(qleft,qright,flx,idim,i)
    use parameters
    use commons

    !!!Does not work with a magnetic field!!!

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    integer  :: idim,idust,i_u,i_v,i_rho,i_w,ipscal,i

    real(dp) :: S_lft,S_rgt,lambda_llf_d


    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt
#if MHD==1

    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,P_mag_lft,P_mag_rgt,mag_tension_y_lft,mag_tension_y_rgt,mag_tension_z_lft,mag_tension_z_rgt,mag_tension_x_lft,mag_tension_x_rgt 


    real(dp) ::flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt

    Bx_lft   = qleft(iBx)
    Bx_rgt   = qright(iBx)
    By_lft   = qleft(iBy)
    By_rgt   = qright(iBy)
    Bz_lft   = qleft(iBz)
    Bz_rgt   = qright(iBz)

    P_mag_lft   = (Bx_lft**2+Bz_lft**2+By_lft**2)/2.0d0 
    P_mag_rgt   = (Bx_rgt**2+Bz_rgt**2+By_rgt**2)/2.0d0
   
    mag_tension_x_lft = -Bx_lft*Bx_lft
    mag_tension_x_rgt = -Bx_rgt*Bx_rgt
    mag_tension_y_lft = -By_lft*Bx_lft
    mag_tension_y_rgt = -By_rgt*Bx_rgt
    mag_tension_z_lft = -Bz_lft*Bx_lft
    mag_tension_z_rgt = -Bz_rgt*Bx_rgt
#endif


    do idust=1,ndust

        i_rho= irhod(idust)
        i_u  = index_vdn(idust,idim)
        i_v  = index_vdt(idust,idim)
        i_w  = ivdz(idust)
        !print *, idust, i_rho,i_n,i_t,i_z
        !Dust density 
        rho_rgt   = qright(i_rho)
        rho_lft   = qleft(i_rho)
        !Dust momentum
        u_rgt     = qright(i_u)
        u_lft     = qleft(i_u)
          !Dust transverse momentum
        v_rgt     = qright(i_v)
        v_lft     = qleft(i_v)
          !Dust second transverse momentum
        w_rgt     = qright(i_w)
        w_lft     = qleft(i_w)





        mom_u_rgt    =  rho_rgt * u_rgt
        mom_u_lft    =  rho_lft * u_lft
        mom_v_rgt    =  rho_rgt * v_rgt
        mom_v_lft    =  rho_lft * v_lft
        mom_w_rgt    =  rho_rgt  * w_rgt
        mom_w_lft    =  rho_lft  * w_lft


        flx_rho_rgt   = rho_rgt  * u_rgt
        flx_rho_lft   = rho_lft  * u_lft
 
        flx_mom_u_rgt  = rho_rgt * u_rgt**2
        flx_mom_u_lft  = rho_lft * u_lft**2

        flx_mom_v_rgt = rho_rgt * u_rgt  * v_rgt
        flx_mom_v_lft = rho_lft * u_lft  * v_lft
           
        flx_mom_w_rgt = rho_rgt * u_rgt  * w_rgt
        flx_mom_w_lft = rho_lft * u_lft  * w_lft

#if MHD==1


        flx_mom_u_rgt  = flx_mom_u_rgt + P_mag_rgt + mag_tension_x_rgt
        flx_mom_u_lft  = flx_mom_u_lft + P_mag_lft + mag_tension_x_lft

        flx_mom_v_rgt = flx_mom_v_rgt + mag_tension_y_rgt
        flx_mom_v_lft = flx_mom_v_lft + mag_tension_y_lft
           
        flx_mom_w_rgt = flx_mom_w_rgt + mag_tension_z_rgt
        flx_mom_w_lft = flx_mom_w_lft + mag_tension_z_lft





#endif

flx(i_rho)  =  0.d0
flx(i_u)    =  0.d0
flx(i_v)    =  0.d0
flx(i_w)    =  0.d0
#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
        flx(idust_pscal(idust,ipscal))  =  0
    end do
#endif 
! Huang & Bai solver (for both MHD==0 and MHD==1)
if(u_rgt>0.0d0 .and. u_lft>0.0d0) then
    flx(i_rho)  =  flx_rho_lft 
    flx(i_u)    =  flx_mom_u_lft 
    flx(i_v)    =  flx_mom_v_lft     
    flx(i_w)    =  flx_mom_w_lft 
#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
        flx(idust_pscal(idust,ipscal))  = qleft(idust_pscal(idust,ipscal))*flx_rho_lft
    end do
#endif 
else if (u_lft<0.0d0 .and. u_rgt<0.0d0) then
    flx(i_rho)  =  flx_rho_rgt
    flx(i_u)    =  flx_mom_u_rgt
    flx(i_v)    =  flx_mom_v_rgt    
    flx(i_w)    =  flx_mom_w_rgt 
#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
        flx(idust_pscal(idust,ipscal))  = qright(idust_pscal(idust,ipscal))*flx_rho_rgt
    end do
#endif 
else if (u_lft<0.0d0 .and. u_rgt>0.0d0) then
    flx(i_rho)  =  0.d0
    flx(i_u)    =  0.d0
    flx(i_v)    =  0.d0
    flx(i_w)    =  0.d0
#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
        flx(idust_pscal(idust,ipscal))  = 0.0d0
    end do
#endif 
else if (u_lft>0.0d0 .and. u_rgt<0.0d0) then
    flx(i_rho)  =  flx_rho_lft   + flx_rho_rgt
    flx(i_u)    =  flx_mom_u_lft + flx_mom_u_rgt
    flx(i_v)    =  flx_mom_v_lft + flx_mom_v_rgt 
    flx(i_w)    =  flx_mom_w_lft + flx_mom_w_rgt 
#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
        flx(idust_pscal(idust,ipscal))  = qleft(idust_pscal(idust,ipscal))*flx_rho_lft + qright(idust_pscal(idust,ipscal))*flx_rho_rgt
    end do
#endif 
endif

end do

end subroutine solver_dust_Huang_Bai
#endif
#endif


#if NDUST>0
#if SOLVERDUST==1
subroutine solver_dust_llf(qleft,qright,flx,idim,i)

    use parameters
    use commons

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    integer  :: idim,idust,i_u,i_v,i_rho,i_w,i

    real(dp) :: S_lft,S_rgt,lambda_llf_d,csl,csr


    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt
#if MHD==1

    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,P_mag_lft,P_mag_rgt,mag_tension_y_lft,mag_tension_y_rgt,mag_tension_z_lft,mag_tension_z_rgt,mag_tension_x_lft,mag_tension_x_rgt 


    real(dp) ::flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt

    Bx_lft   = qleft(iBx)
    Bx_rgt   = qright(iBx)
    By_lft   = qleft(iBy)
    By_rgt   = qright(iBy)
    Bz_lft   = qleft(iBz)
    Bz_rgt   = qright(iBz)

    P_mag_lft   = (Bx_lft**2+Bz_lft**2+By_lft**2)/2.0d0 
    P_mag_rgt   = (Bx_rgt**2+Bz_rgt**2+By_rgt**2)/2.0d0
   
    mag_tension_x_lft = -Bx_lft*Bx_lft
    mag_tension_x_rgt = -Bx_rgt*Bx_rgt
    mag_tension_y_lft = -By_lft*Bx_lft
    mag_tension_y_rgt = -By_rgt*Bx_rgt
    mag_tension_z_lft = -Bz_lft*Bx_lft
    mag_tension_z_rgt = -Bz_rgt*Bx_rgt
#endif


    do idust=1,ndust

        i_rho= irhod(idust)
        i_u  = index_vdn(idust,idim)
        i_v  = index_vdt(idust,idim)
        i_w  = ivdz(idust)
        !print *, idust, i_rho,i_n,i_t,i_z
        !Dust density 
        rho_rgt   = qright(i_rho)
        rho_lft   = qleft(i_rho)
        !Dust momentum
        u_rgt     = qright(i_u)
        u_lft     = qleft(i_u)
          !Dust transverse momentum
        v_rgt     = qright(i_v)
        v_lft     = qleft(i_v)
          !Dust second transverse momentum
        w_rgt     = qright(i_w)
        w_lft     = qleft(i_w)





        mom_u_rgt    =  rho_rgt * u_rgt
        mom_u_lft    =  rho_lft * u_lft
        mom_v_rgt    =  rho_rgt * v_rgt
        mom_v_lft    =  rho_lft * v_lft
        mom_w_rgt    =  rho_rgt  * w_rgt
        mom_w_lft    =  rho_lft  * w_lft


        flx_rho_rgt   = rho_rgt  * u_rgt
        flx_rho_lft   = rho_lft  * u_lft
 
        flx_mom_u_rgt  = rho_rgt * u_rgt**2
        flx_mom_u_lft  = rho_lft * u_lft**2

        flx_mom_v_rgt = rho_rgt * u_rgt  * v_rgt
        flx_mom_v_lft = rho_lft * u_lft  * v_lft
           
        flx_mom_w_rgt = rho_rgt * u_rgt  * w_rgt
        flx_mom_w_lft = rho_lft * u_lft  * w_lft

#if MHD==1


        flx_mom_u_rgt  = flx_mom_u_rgt + 1/(4*pi)*(P_mag_rgt + mag_tension_x_rgt)
        flx_mom_u_lft  = flx_mom_u_lft + 1/(4*pi)*(P_mag_lft + mag_tension_x_lft)

        flx_mom_v_rgt = flx_mom_v_rgt + 1/(4*pi)*mag_tension_y_rgt
        flx_mom_v_lft = flx_mom_v_lft + 1/(4*pi)*mag_tension_y_lft
           
        flx_mom_w_rgt = flx_mom_w_rgt + 1/(4*pi)*mag_tension_z_rgt
        flx_mom_w_lft = flx_mom_w_lft + 1/(4*pi)*mag_tension_z_lft



#endif

#if MHD==0

     lambda_llf_d        = max(abs(u_lft),abs(u_rgt))

#endif

#if MHD==1
    
     lambda_llf_d        = max(abs(u_lft)+dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft),abs(u_rgt)+dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt))

#endif

    flx(i_rho)  =  0.d0
    flx(i_u)    =  0.d0
    flx(i_v)    =  0.d0
    flx(i_w)    =  0.d0

    flx(i_rho) = half*(flx_rho_lft  + flx_rho_rgt)   - half*lambda_llf_d*(rho_rgt - rho_lft)
    flx(i_u)  = half*(flx_mom_u_lft + flx_mom_u_rgt)  - half*lambda_llf_d*(mom_u_rgt - mom_u_lft)
    flx(i_v)  = half*(flx_mom_v_lft + flx_mom_v_rgt)  - half*lambda_llf_d*(mom_v_rgt - mom_v_lft)
    flx(i_w)  = half*(flx_mom_w_lft + flx_mom_w_rgt)  - half*lambda_llf_d*(mom_w_rgt - mom_w_lft)

end do
end subroutine solver_dust_llf

#endif
#endif


#if NDUST>0
#if SOLVERDUST==2

subroutine solver_dust_hll(qleft,qright,flx,idim,i)

    use parameters
    use commons

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    integer  :: idim,idust,i_u,i_v,i_rho,i_w,i

    real(dp) :: S_lft,S_rgt,lambda_llf_d
    real(dp) :: ca_lft,ca_rgt



    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt
#if MHD==1

    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,P_mag_lft,P_mag_rgt,mag_tension_y_lft,mag_tension_y_rgt,mag_tension_z_lft,mag_tension_z_rgt,mag_tension_x_lft,mag_tension_x_rgt 


    real(dp) ::flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt

    Bx_lft   = qleft(iBx)
    Bx_rgt   = qright(iBx)
    By_lft   = qleft(iBy)
    By_rgt   = qright(iBy)
    Bz_lft   = qleft(iBz)
    Bz_rgt   = qright(iBz)

    P_mag_lft   = (Bx_lft**2+Bz_lft**2+By_lft**2)/2.0d0 
    P_mag_rgt   = (Bx_rgt**2+Bz_rgt**2+By_rgt**2)/2.0d0
   
    mag_tension_x_lft = -Bx_lft*Bx_lft
    mag_tension_x_rgt = -Bx_rgt*Bx_rgt
    mag_tension_y_lft = -By_lft*Bx_lft
    mag_tension_y_rgt = -By_rgt*Bx_rgt
    mag_tension_z_lft = -Bz_lft*Bx_lft
    mag_tension_z_rgt = -Bz_rgt*Bx_rgt
#endif


    do idust=1,ndust

        i_rho= irhod(idust)
        i_u  = index_vdn(idust,idim)
        i_v  = index_vdt(idust,idim)
        i_w  = ivdz(idust)
        !print *, idust, i_rho,i_n,i_t,i_z
        !Dust density 
        rho_rgt   = qright(i_rho)
        rho_lft   = qleft(i_rho)
        !Dust momentum
        u_rgt     = qright(i_u)
        u_lft     = qleft(i_u)
          !Dust transverse momentum
        v_rgt     = qright(i_v)
        v_lft     = qleft(i_v)
          !Dust second transverse momentum
        w_rgt     = qright(i_w)
        w_lft     = qleft(i_w)





        mom_u_rgt    =  rho_rgt * u_rgt
        mom_u_lft    =  rho_lft * u_lft
        mom_v_rgt    =  rho_rgt * v_rgt
        mom_v_lft    =  rho_lft * v_lft
        mom_w_rgt    =  rho_rgt  * w_rgt
        mom_w_lft    =  rho_lft  * w_lft


        flx_rho_rgt   = rho_rgt  * u_rgt
        flx_rho_lft   = rho_lft  * u_lft
 
        flx_mom_u_rgt  = rho_rgt * u_rgt**2
        flx_mom_u_lft  = rho_lft * u_lft**2

        flx_mom_v_rgt = rho_rgt * u_rgt  * v_rgt
        flx_mom_v_lft = rho_lft * u_lft  * v_lft
           
        flx_mom_w_rgt = rho_rgt * u_rgt  * w_rgt
        flx_mom_w_lft = rho_lft * u_lft  * w_lft

#if MHD==1

        flx_rho_rgt   = rho_rgt  * u_rgt
        flx_rho_lft   = rho_lft  * u_lft

        flx_mom_u_rgt  = flx_mom_u_rgt + 1/(4*pi)*(P_mag_rgt + mag_tension_x_rgt)
        flx_mom_u_lft  = flx_mom_u_lft + 1/(4*pi)*(P_mag_lft + mag_tension_x_lft)

        flx_mom_v_rgt = flx_mom_v_rgt + 1/(4*pi)*mag_tension_y_rgt
        flx_mom_v_lft = flx_mom_v_lft + 1/(4*pi)*mag_tension_y_lft
           
        flx_mom_w_rgt = flx_mom_w_rgt + 1/(4*pi)*mag_tension_z_rgt
        flx_mom_w_lft = flx_mom_w_lft + 1/(4*pi)*mag_tension_z_lft


#endif



#if MHD==0
    S_rgt  = max(max(u_lft,u_rgt),0.0d0) 
    S_lft  = min(min(u_lft,u_rgt),0.0d0)

    if (u_lft/=u_rgt) then

    flx(i_rho)            = (S_rgt*flx_rho_lft  -S_lft*flx_rho_rgt  + S_rgt*S_lft*(rho_rgt-rho_lft))      / (S_rgt-S_lft)
    flx(i_u)  = (S_rgt*flx_mom_u_lft-S_lft*flx_mom_u_rgt+ S_rgt*S_lft*(mom_u_rgt-mom_u_lft))  / (S_rgt-S_lft)
    flx(i_v)  = (S_rgt*flx_mom_v_lft-S_lft*flx_mom_v_rgt+ S_rgt*S_lft*(mom_v_rgt-mom_v_lft))  / (S_rgt-S_lft)
    flx(i_w)             = (S_rgt*flx_mom_w_lft-S_lft*flx_mom_w_rgt+ S_rgt*S_lft*(mom_w_rgt-mom_w_lft))  / (S_rgt-S_lft)
    end if

    if (u_lft==u_rgt) then !Switch back to llf

    lambda_llf_d        = max(abs(u_lft),abs(u_rgt))

    flx(i_rho) = half*(flx_rho_lft   + flx_rho_rgt)    - half*lambda_llf_d*(rho_rgt   - rho_lft)
    flx(i_u)  = half*(flx_mom_u_lft + flx_mom_u_rgt)  - half*lambda_llf_d*(mom_u_rgt - mom_u_lft)
    flx(i_v)  = half*(flx_mom_v_lft + flx_mom_v_rgt)  - half*lambda_llf_d*(mom_v_rgt - mom_v_lft)
    flx(i_w)  = half*(flx_mom_w_lft + flx_mom_w_rgt)  - half*lambda_llf_d*(mom_w_rgt - mom_w_lft)



    end if
#endif

#if MHD==1

!HLL    

    ca_lft =dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft)
    ca_rgt =dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt)

    S_rgt  = max(max(u_lft,u_rgt) +max(ca_lft,ca_rgt),0.0d0) 
    S_lft  = min(min(u_lft,u_rgt) -max(ca_lft,ca_rgt),0.0d0)

    flx(i_rho)            = (S_rgt*flx_rho_lft  -S_lft*flx_rho_rgt  + S_rgt*S_lft*(rho_rgt-rho_lft))      / (S_rgt-S_lft)
    flx(i_u)  = (S_rgt*flx_mom_u_lft-S_lft*flx_mom_u_rgt+ S_rgt*S_lft*(mom_u_rgt-mom_u_lft))  / (S_rgt-S_lft)
    flx(i_v)  = (S_rgt*flx_mom_v_lft-S_lft*flx_mom_v_rgt+ S_rgt*S_lft*(mom_v_rgt-mom_v_lft))  / (S_rgt-S_lft)
    flx(i_w)             = (S_rgt*flx_mom_w_lft-S_lft*flx_mom_w_rgt+ S_rgt*S_lft*(mom_w_rgt-mom_w_lft))  / (S_rgt-S_lft)

#endif


end do

end subroutine solver_dust_hll

#endif 
#endif


#if MHD==1
#if SOLVERB==0

subroutine solver_induction_llf(qleft,qright,flx,csl,csr,idim,i)
    use parameters
    use commons
    use slope_limiter

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    integer  :: idim,idust,i_u,i_v,i_rho,i_w,i,il,ir,ix,iy,icell,ixx,iyy

    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,lambda_llf_B
    real(dp) :: u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt,rho_rgt,rho_lft 
    real(dp) :: magnetosonic_fast_rgt,magnetosonic_fast_lft,csl,csr
    real(dp) :: deta_o,deta_h,deta_a,deta_o_il,deta_h_il,deta_a_il,eta_o_left,eta_h_left,eta_a_left,eta_o_right,eta_h_right,eta_a_right,dzd,dzd_il,zd_left,zd_right
    real(dp) :: dhall_i,dhall_il,hall_i_left,hall_i_right
    real(dp) ::flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt,total_dust_current_z_lft,total_dust_current_y_lft,total_dust_current_x_lft,total_dust_current_z_rgt,total_dust_current_y_rgt,total_dust_current_x_rgt,B_norm_lft,B_norm_rgt


    Bx_lft   = qleft(iBx)
    Bx_rgt   = qright(iBx)
    By_lft   = qleft(iBy)
    By_rgt   = qright(iBy)
    Bz_lft   = qleft(iBz)
    Bz_rgt   = qright(iBz)

  

#if NDUST>0

        idust=i_coupled_species !Is the grain species considered to be ideally coupled to B

        i_rho= irhod(idust)

        i_u  = index_vdn(idust,idim)
        i_v  = index_vdt(idust,idim)
        i_w  = ivdz(idust)
        !print *, idust, i_rho,i_n,i_t,i_z

        !Dust momentum
        u_rgt     = qright(i_u)
        u_lft     = qleft(i_u)
          !Dust transverse momentum
        v_rgt     = qright(i_v)
        v_lft     = qleft(i_v)
          !Dust second transverse momentum
        w_rgt     = qright(i_w)
        w_lft     = qleft(i_w)

        rho_rgt   = qright(i_rho)
        rho_lft   = qleft(i_rho)

#endif
!B coupled to the gas
#if NDUST==0 
        rho_rgt   = qright(irho)
        rho_lft   = qleft(irho)

        !Velocity
        u_rgt   = qright(index_vn(idim)) ! u
        u_lft   = qleft(index_vn(idim))
        !Transverse velocity
        v_rgt   = qright(index_vt(idim)) ! v
        v_lft   = qleft(index_vt(idim))
        w_rgt   = qright(ivz)! w
        w_lft   = qleft(ivz)


#endif


    flx_Bx_lft = 0.0d0
    flx_Bx_rgt = 0.0d0

    flx_By_lft = By_lft*u_lft - Bx_lft*v_lft
    flx_By_rgt = By_rgt*u_rgt - Bx_rgt*v_rgt

    flx_Bz_lft = Bz_lft*u_lft - Bx_lft*w_lft
    flx_Bz_rgt = Bz_rgt*u_rgt - Bx_rgt*w_rgt    



    flx(iBx)    = 0.0d0



#if NDUST>0



      if (dusty_nonideal_MHD_no_electron) then !Additional terms in the fluxes for the induction equation

            idust = i_coupled_species
        
            i_rho= irhod(idust)

            i_u  = index_vdn(idust,idim)
            i_v  = index_vdt(idust,idim)
            i_w  = ivdz(idust)
            !print *, idust, i_rho,i_n,i_t,i_z

            !Dust momentum
            u_rgt     = qright(i_u)
            u_lft     = qleft(i_u)
              !Dust transverse momentum
            v_rgt     = qright(i_v)
            v_lft     = qleft(i_v)
              !Dust second transverse momentum
            w_rgt     = qright(i_w)
            w_lft     = qleft(i_w)

            rho_rgt   = qright(i_rho)
            rho_lft   = qleft(i_rho)


            ix=ixx(i)
            iy=iyy(i)

            if(slope_type>0) then
                il = icell(ix-1,iy)
                ir = icell(ix+1,iy)
                ! dzd = slope_limit(2.0d0*(zd(i,idust) - zd(il,idust))/(dx(i,1)+dx(il,1)),2.0d0*(zd(ir,idust) - zd(i,idust))/(dx(ir,1)+dx(i,1)))
                ! dzd_il = slope_limit(2.0d0*(zd(i-1,idust) - zd(il-1,idust))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(zd(ir-1,idust) - zd(i-1,idust))/(dx(ir-1,1)+dx(i-1,1)))
                ! zd_left = zd(il,idust) + half*dzd_il*dx(il,1)
                ! zd_right = zd(i,idust) - half*dzd*dx(i,1)

                ! dni = slope_limit(2.0d0*(ni(i) - ni(il))/(dx(i,1)+dx(il,1)),2.0d0*(ni(ir) - ni(i))/(dx(ir,1)+dx(i,1)))
                ! dni_il = slope_limit(2.0d0*(ni(i-1) - ni(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(ni(ir-1) - ni(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                ! ni_left = ni(il) + half*dni_il*dx(il,1)
                ! ni_right = ni(i) - half*dni*dx(i,1)

           
                dhall_i = slope_limit(2.0d0*(Hall_i(i) - Hall_i(il))/(dx(i,1)+dx(il,1)),2.0d0*(Hall_i(ir) - Hall_i(i))/(dx(ir,1)+dx(i,1)))
                dhall_il = slope_limit(2.0d0*(Hall_i(i-1) - Hall_i(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Hall_i(ir-1) - Hall_i(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                hall_i_left = Hall_i(il) + half*dhall_il*dx(il,1)
                hall_i_right = Hall_i(i) - half*dhall_i*dx(i,1)


            endif


        B_norm_lft = dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)
        B_norm_rgt = dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2) 

        !nd_zd_over_ni_left = (rho_lft/mdust(i,idust))*zd_left/(ni_left-ne_left)
        !nd_zd_over_ni_right = (rho_rgt/mdust(i,idust))*zd_right/(ni_right-ne_right)

        !nd_zd_over_ni = (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/(ni(i)-ne(i))
        !nd_zd_over_ni = -1.0d0 !Due to electroneutrality


        !Beware: here signs are reversed with respect to predictor step (fluxes are defined in the left-hand side of the equation) 

        flx_By_lft = -Bx_lft*v_lft + By_lft*u_lft
        flx_By_rgt = -Bx_rgt*v_rgt + By_rgt*u_rgt

        flx_Bz_lft = -Bx_lft*w_lft + Bz_lft*u_lft
        flx_Bz_rgt = -Bx_rgt*w_rgt + Bz_rgt*u_rgt


        ! !Additional term
        flx_By_lft = flx_By_lft - B_norm_lft/hall_i_left*(w_lft - qleft(ivz))
        flx_By_rgt = flx_By_rgt - B_norm_rgt/hall_i_right*(w_rgt - qright(ivz))

        flx_Bz_lft = flx_Bz_lft + B_norm_lft/hall_i_left*(v_lft - qleft(index_vt(idim)))
        flx_Bz_rgt = flx_Bz_rgt + B_norm_rgt/hall_i_right*(v_rgt - qright(index_vt(idim)))

    endif


#endif


#if NDUST==0
magnetosonic_fast_rgt = dsqrt(half*(csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt + dsqrt((csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt)**2-4*csr**2*Bx_rgt**2/rho_rgt))) 
magnetosonic_fast_lft = dsqrt(half*(csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft + dsqrt((csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft)**2-4*csl**2*Bx_lft**2/rho_lft))) 
lambda_llf_B        = max(abs(u_lft)+magnetosonic_fast_lft,abs(u_rgt)+magnetosonic_fast_rgt)

#endif

#if NDUST>0

    idust=i_coupled_species !Is the grain species considered in the magnetosonic/Alfven velocity expressions

    i_rho= irhod(idust)

    i_u  = index_vdn(idust,idim)
    i_v  = index_vdt(idust,idim)
    i_w  = ivdz(idust)
    !print *, idust, i_rho,i_n,i_t,i_z

    !Dust momentum
    u_rgt     = qright(i_u)
    u_lft     = qleft(i_u)
      !Dust transverse momentum
    v_rgt     = qright(i_v)
    v_lft     = qleft(i_v)
      !Dust second transverse momentum
    w_rgt     = qright(i_w)
    w_lft     = qleft(i_w)

    rho_rgt   = qright(i_rho)
    rho_lft   = qleft(i_rho)
    
    lambda_llf_B        = max(abs(u_lft)+dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft),abs(u_rgt)+dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt))
#endif


flx(iBy) = half*(flx_By_lft   + flx_By_rgt)    - half*lambda_llf_B*(By_rgt   - By_lft)
flx(iBz) = half*(flx_Bz_lft   + flx_Bz_rgt)    - half*lambda_llf_B*(Bz_rgt   - Bz_lft)

end subroutine solver_induction_llf

#endif
#endif


#if MHD==1
#if SOLVERB==1

subroutine solver_induction_Huang_Bai(qleft,qright,flx,idim,i)
    use parameters
    use commons
    use slope_limiter

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    integer  :: idim,idust,i_u,i_v,i_rho,i_w,i,il,ir,ix,iy,icell,ixx,iyy

    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,lambda_llf_B
    real(dp) :: u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt,rho_lft,rho_rgt
    real(dp) :: dhall_i,dhall_il,hall_i_left,hall_i_right
    real(dp) ::flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt,total_dust_current_z_lft,total_dust_current_y_lft,total_dust_current_x_lft,total_dust_current_z_rgt,total_dust_current_y_rgt,total_dust_current_x_rgt,B_norm_lft,B_norm_rgt,B_norm

    Bx_lft   = qleft(iBx)
    Bx_rgt   = qright(iBx)
    By_lft   = qleft(iBy)
    By_rgt   = qright(iBy)
    Bz_lft   = qleft(iBz)
    Bz_rgt   = qright(iBz)

#if NDUST>0

        idust=i_coupled_species !Is the grain species considered to be ideally coupled to B

        i_u  = index_vdn(idust,idim)
        i_v  = index_vdt(idust,idim)
        i_w  = ivdz(idust)
        !print *, idust, i_rho,i_n,i_t,i_z

        !Dust momentum
        u_rgt     = qright(i_u)
        u_lft     = qleft(i_u)
          !Dust transverse momentum
        v_rgt     = qright(i_v)
        v_lft     = qleft(i_v)
          !Dust second transverse momentum
        w_rgt     = qright(i_w)
        w_lft     = qleft(i_w)

#endif
!B coupled to the gas
#if NDUST==0 

        rho_rgt   = qright(irho)
        rho_lft   = qleft(irho)
        !Velocity
        u_rgt   = qright(index_vn(idim)) ! u
        u_lft   = qleft(index_vn(idim))
        !Transverse velocity
        v_rgt   = qright(index_vt(idim)) ! v
        v_lft   = qleft(index_vt(idim))
        w_rgt   = qright(ivz)! w
        w_lft   = qleft(ivz)


#endif


    flx_Bx_lft = 0.0d0
    flx_Bx_rgt = 0.0d0

    flx_By_lft = By_lft*u_lft - Bx_lft*v_lft
    flx_By_rgt = By_rgt*u_rgt - Bx_rgt*v_rgt

    flx_Bz_lft = Bz_lft*u_lft - Bx_lft*w_lft
    flx_Bz_rgt = Bz_rgt*u_rgt - Bx_rgt*w_rgt    



    flx(iBx)    = 0.d0



#if NDUST>0


      if (dusty_nonideal_MHD_no_electron) then !Additional terms in the fluxes for the induction equation

            idust = i_coupled_species
        
            i_rho= irhod(idust)

            i_u  = index_vdn(idust,idim)
            i_v  = index_vdt(idust,idim)
            i_w  = ivdz(idust)
            !print *, idust, i_rho,i_n,i_t,i_z

            !Dust momentum
            u_rgt     = qright(i_u)
            u_lft     = qleft(i_u)
              !Dust transverse momentum
            v_rgt     = qright(i_v)
            v_lft     = qleft(i_v)
              !Dust second transverse momentum
            w_rgt     = qright(i_w)
            w_lft     = qleft(i_w)

            rho_rgt   = qright(i_rho)
            rho_lft   = qleft(i_rho)


            ix=ixx(i)
            iy=iyy(i)

            if(slope_type>0) then
                il = icell(ix-1,iy)
                ir = icell(ix+1,iy)
                ! dzd = slope_limit(2.0d0*(zd(i,idust) - zd(il,idust))/(dx(i,1)+dx(il,1)),2.0d0*(zd(ir,idust) - zd(i,idust))/(dx(ir,1)+dx(i,1)))
                ! dzd_il = slope_limit(2.0d0*(zd(i-1,idust) - zd(il-1,idust))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(zd(ir-1,idust) - zd(i-1,idust))/(dx(ir-1,1)+dx(i-1,1)))
                ! zd_left = zd(il,idust) + half*dzd_il*dx(il,1)
                ! zd_right = zd(i,idust) - half*dzd*dx(i,1)

                ! dni = slope_limit(2.0d0*(ni(i) - ni(il))/(dx(i,1)+dx(il,1)),2.0d0*(ni(ir) - ni(i))/(dx(ir,1)+dx(i,1)))
                ! dni_il = slope_limit(2.0d0*(ni(i-1) - ni(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(ni(ir-1) - ni(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                ! ni_left = ni(il) + half*dni_il*dx(il,1)
                ! ni_right = ni(i) - half*dni*dx(i,1)

             
                dhall_i = slope_limit(2.0d0*(Hall_i(i) - Hall_i(il))/(dx(i,1)+dx(il,1)),2.0d0*(Hall_i(ir) - Hall_i(i))/(dx(ir,1)+dx(i,1)))
                dhall_il = slope_limit(2.0d0*(Hall_i(i-1) - Hall_i(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Hall_i(ir-1) - Hall_i(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                hall_i_left = Hall_i(il) + half*dhall_il*dx(il,1)
                hall_i_right = Hall_i(i) - half*dhall_i*dx(i,1)


            endif


        B_norm_lft = dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)
        B_norm_rgt = dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2) 

        !nd_zd_over_ni_left = (rho_lft/mdust(i,idust))*zd_left/(ni_left-ne_left)
        !nd_zd_over_ni_right = (rho_rgt/mdust(i,idust))*zd_right/(ni_right-ne_right)

        !nd_zd_over_ni = (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/(ni(i)-ne(i))
        !nd_zd_over_ni = -1.0d0 !Due to electroneutrality


        !Beware: here signs are reversed with respect to predictor step (fluxes are defined in the left-hand side of the equation) 



        flx_By_lft = -Bx_lft*v_lft + By_lft*u_lft
        flx_By_rgt = -Bx_rgt*v_rgt + By_rgt*u_rgt

        flx_Bz_lft = -Bx_lft*w_lft + Bz_lft*u_lft
        flx_Bz_rgt = -Bx_rgt*w_rgt + Bz_rgt*u_rgt


        ! !Additional term
        flx_By_lft = flx_By_lft - B_norm_lft/hall_i_left*(w_lft - qleft(ivz))
        flx_By_rgt = flx_By_rgt - B_norm_rgt/hall_i_right*(w_rgt - qright(ivz))

        flx_Bz_lft = flx_Bz_lft + B_norm_lft/hall_i_left*(v_lft - qleft(index_vt(idim)))
        flx_Bz_rgt = flx_Bz_rgt + B_norm_rgt/hall_i_right*(v_rgt - qright(index_vt(idim)))
     
    endif


#endif



! Huang & Bai solver 
if(u_rgt>0.0d0 .and. u_lft>0.0d0) then
    flx(iBy)    =  flx_By_lft 
    flx(iBz)    =  flx_Bz_lft 
else if (u_lft<0.0d0 .and. u_rgt<0.0d0) then
    flx(iBy)    =  flx_By_rgt 
    flx(iBz)    =  flx_Bz_rgt
else if (u_lft<0.0d0 .and. u_rgt>0.0d0) then
    flx(iBy)    =  0.0d0
    flx(iBz)    =  0.0d0
else if (u_lft>0.0d0 .and. u_rgt<0.0d0) then
    flx(iBy)    =  flx_By_lft+flx_By_rgt 
    flx(iBz)    =  flx_Bz_lft+flx_Bz_rgt
endif

end subroutine solver_induction_Huang_Bai
#endif
#endif


#if MHD==1
#if SOLVERB==2

subroutine solver_induction_hll(qleft,qright,flx,csl,csr,idim,i)
    use parameters
    use commons
    use slope_limiter

    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    integer  :: idim,idust,i_u,i_v,i_rho,i_w,i,il,ir,ix,iy,icell,ixx,iyy

    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,lambda_llf_B
    real(dp) :: u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt,rho_rgt,rho_lft
    real(dp) :: S_lft,S_rgt
    real(dp) :: ca_lft,ca_rgt,magnetosonic_fast_lft,magnetosonic_fast_rgt,csl,csr

    real(dp) :: deta_o,deta_h,deta_a,deta_o_il,deta_h_il,deta_a_il,eta_o_left,eta_h_left,eta_a_left,eta_o_right,eta_h_right,eta_a_right,dzd,dzd_il,zd_left,zd_right

    real(dp) :: nd_zd_over_ni,nd_zd_over_ni_left,nd_zd_over_ni_right,dne,dne_il,ne_left,ne_right,dni,dni_il,ni_left,ni_right,dhall_i,dhall_il,hall_i_left,hall_i_right
    real(dp) :: flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt,total_dust_current_z_lft,total_dust_current_y_lft,total_dust_current_x_lft,total_dust_current_z_rgt,total_dust_current_y_rgt,total_dust_current_x_rgt,B_norm_lft,B_norm_rgt


    Bx_lft   = qleft(iBx)
    Bx_rgt   = qright(iBx)
    By_lft   = qleft(iBy)
    By_rgt   = qright(iBy)
    Bz_lft   = qleft(iBz)
    Bz_rgt   = qright(iBz)




#if NDUST>0

        idust=i_coupled_species !Is the grain species considered to be ideally coupled to B

        i_rho=irhod(idust)
        i_u  = index_vdn(idust,idim)
        i_v  = index_vdt(idust,idim)
        i_w  = ivdz(idust)
        !print *, idust, i_rho,i_n,i_t,i_z

        !Dust momentum
        u_rgt     = qright(i_u)
        u_lft     = qleft(i_u)
          !Dust transverse momentum
        v_rgt     = qright(i_v)
        v_lft     = qleft(i_v)
          !Dust second transverse momentum
        w_rgt     = qright(i_w)
        w_lft     = qleft(i_w)

        rho_rgt   = qright(i_rho)
        rho_lft   = qleft(i_rho)


#endif
!B coupled to the gas
#if NDUST==0

        rho_rgt   = qright(irho)
        rho_lft   = qleft(irho)
        !Velocity
        u_rgt   = qright(index_vn(idim)) ! u
        u_lft   = qleft(index_vn(idim))
        !Transverse velocity
        v_rgt   = qright(index_vt(idim)) ! v
        v_lft   = qleft(index_vt(idim))
        w_rgt   = qright(ivz)! w
        w_lft   = qleft(ivz)


#endif


    flx_Bx_lft = 0.0d0
    flx_Bx_rgt = 0.0d0
   
    flx_By_lft = By_lft*u_lft - Bx_lft*v_lft
    flx_By_rgt = By_rgt*u_rgt - Bx_rgt*v_rgt

    flx_Bz_lft = Bz_lft*u_lft - Bx_lft*w_lft
    flx_Bz_rgt = Bz_rgt*u_rgt - Bx_rgt*w_rgt    




    flx(iBx)    = 0.d0


    ! flx_By_lft = By_lft*qleft(index_vn(idim)) - Bx_lft*qleft(index_vt(idim))
    ! flx_By_rgt = By_rgt*qright(index_vn(idim)) - Bx_rgt*qright(index_vt(idim)) !!Velocity of the neutrals (gas)

    ! flx_Bz_lft = Bz_lft*qleft(index_vn(idim)) - Bx_lft*qleft(ivz)
    ! flx_Bz_rgt = Bz_rgt*qright(index_vn(idim)) - Bx_rgt*qright(ivz)




#if NDUST>0


      if (dusty_nonideal_MHD_no_electron) then !Additional terms in the fluxes for the induction equation

            idust = i_coupled_species
        
            i_rho= irhod(idust)

            i_u  = index_vdn(idust,idim)
            i_v  = index_vdt(idust,idim)
            i_w  = ivdz(idust)
            !print *, idust, i_rho,i_n,i_t,i_z

            !Dust momentum
            u_rgt     = qright(i_u)
            u_lft     = qleft(i_u)
              !Dust transverse momentum
            v_rgt     = qright(i_v)
            v_lft     = qleft(i_v)
              !Dust second transverse momentum
            w_rgt     = qright(i_w)
            w_lft     = qleft(i_w)

            rho_rgt   = qright(i_rho)
            rho_lft   = qleft(i_rho)


            ix=ixx(i)
            iy=iyy(i)

            if(slope_type>0) then
                il = icell(ix-1,iy)
                ir = icell(ix+1,iy)
                ! dzd = slope_limit(2.0d0*(zd(i,idust) - zd(il,idust))/(dx(i,1)+dx(il,1)),2.0d0*(zd(ir,idust) - zd(i,idust))/(dx(ir,1)+dx(i,1)))
                ! dzd_il = slope_limit(2.0d0*(zd(i-1,idust) - zd(il-1,idust))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(zd(ir-1,idust) - zd(i-1,idust))/(dx(ir-1,1)+dx(i-1,1)))
                ! zd_left = zd(il,idust) + half*dzd_il*dx(il,1)
                ! zd_right = zd(i,idust) - half*dzd*dx(i,1)

                ! dni = slope_limit(2.0d0*(ni(i) - ni(il))/(dx(i,1)+dx(il,1)),2.0d0*(ni(ir) - ni(i))/(dx(ir,1)+dx(i,1)))
                ! dni_il = slope_limit(2.0d0*(ni(i-1) - ni(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(ni(ir-1) - ni(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                ! ni_left = ni(il) + half*dni_il*dx(il,1)
                ! ni_right = ni(i) - half*dni*dx(i,1)

                dhall_i = slope_limit(2.0d0*(Hall_i(i) - Hall_i(il))/(dx(i,1)+dx(il,1)),2.0d0*(Hall_i(ir) - Hall_i(i))/(dx(ir,1)+dx(i,1)))
                dhall_il = slope_limit(2.0d0*(Hall_i(i-1) - Hall_i(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Hall_i(ir-1) - Hall_i(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                hall_i_left = Hall_i(il) + half*dhall_il*dx(il,1)
                hall_i_right = Hall_i(i) - half*dhall_i*dx(i,1)


            endif


        B_norm_lft = dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)
        B_norm_rgt = dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2) 

        !nd_zd_over_ni_left = (rho_lft/mdust(i,idust))*zd_left/(ni_left-ne_left)
        !nd_zd_over_ni_right = (rho_rgt/mdust(i,idust))*zd_right/(ni_right-ne_right)

        !nd_zd_over_ni = (q(i,irhod(idust))/mdust(i,idust))*zd(i,idust)/(ni(i)-ne(i))
        !nd_zd_over_ni = -1.0d0 !Due to electroneutrality


        !Beware: here signs are reversed with respect to predictor step (fluxes are defined in the left-hand side of the equation) 


        flx_By_lft = -Bx_lft*v_lft + By_lft*u_lft
        flx_By_rgt = -Bx_rgt*v_rgt + By_rgt*u_rgt

        flx_Bz_lft = -Bx_lft*w_lft + Bz_lft*u_lft
        flx_Bz_rgt = -Bx_rgt*w_rgt + Bz_rgt*u_rgt


        !Additional term
        flx_By_lft = flx_By_lft - B_norm_lft/hall_i_left*(w_lft - qleft(ivz))
        flx_By_rgt = flx_By_rgt - B_norm_rgt/hall_i_right*(w_rgt - qright(ivz))

        flx_Bz_lft = flx_Bz_lft + B_norm_lft/hall_i_left*(v_lft - qleft(index_vt(idim)))
        flx_Bz_rgt = flx_Bz_rgt + B_norm_rgt/hall_i_right*(v_rgt - qright(index_vt(idim)))

        if (only_Hall_effect) then

            flx_By_lft = -Bx_lft*v_lft + By_lft*u_lft
            flx_By_rgt = -Bx_rgt*v_rgt + By_rgt*u_rgt

            flx_Bz_lft = -Bx_lft*w_lft + Bz_lft*u_lft
            flx_Bz_rgt = -Bx_rgt*w_rgt + Bz_rgt*u_rgt


        endif
     


    endif


#endif



!HLL

#if NDUST>0

    idust=i_coupled_species !Is the grain species considered in the magnetosonic/Alfven velocity expressions

    i_rho= irhod(idust)

    i_u  = index_vdn(idust,idim)
    i_v  = index_vdt(idust,idim)
    i_w  = ivdz(idust)
    !print *, idust, i_rho,i_n,i_t,i_z

    !Dust momentum
    u_rgt     = qright(i_u)
    u_lft     = qleft(i_u)
      !Dust transverse momentum
    v_rgt     = qright(i_v)
    v_lft     = qleft(i_v)
      !Dust second transverse momentum
    w_rgt     = qright(i_w)
    w_lft     = qleft(i_w)

    rho_rgt   = qright(i_rho)
    rho_lft   = qleft(i_rho)




    ca_lft =dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft)
    ca_rgt =dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt)

    ! magnetosonic_fast_rgt = dsqrt(half*(csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt + dsqrt((csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt)**2-4*csr**2*Bx_rgt**2/rho_rgt))) 
    ! magnetosonic_fast_lft = dsqrt(half*(csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft + dsqrt((csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft)**2-4*csl**2*Bx_lft**2/rho_lft))) 


    S_rgt  = max(max(u_lft,u_rgt) +max(ca_lft,ca_rgt),0.0d0) 
    S_lft  = min(min(u_lft,u_rgt) -max(ca_lft,ca_rgt),0.0d0)

    ! S_rgt  = max(max(u_lft,u_rgt) +max(magnetosonic_fast_lft,magnetosonic_fast_rgt),0.0d0) 
    ! S_lft  = min(min(u_lft,u_rgt) -max(magnetosonic_fast_lft,magnetosonic_fast_rgt),0.0d0)
#endif

#if NDUST==0
    magnetosonic_fast_rgt = dsqrt(half*(csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt + dsqrt((csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt)**2-4*csr**2*Bx_rgt**2/rho_rgt))) 
    magnetosonic_fast_lft = dsqrt(half*(csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft + dsqrt((csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft)**2-4*csl**2*Bx_lft**2/rho_lft))) 
    S_rgt  = max(max(u_lft,u_rgt) +max(magnetosonic_fast_lft,magnetosonic_fast_rgt),0.0d0) 
    S_lft  = min(min(u_lft,u_rgt) -max(magnetosonic_fast_lft,magnetosonic_fast_rgt),0.0d0)
#endif
    flx(iBy)            = (S_rgt*flx_By_lft  -S_lft*flx_By_rgt  + S_rgt*S_lft*(By_rgt-By_lft))      / (S_rgt-S_lft)
    flx(iBz)            = (S_rgt*flx_Bz_lft  -S_lft*flx_Bz_rgt  + S_rgt*S_lft*(Bz_rgt-Bz_lft))      / (S_rgt-S_lft)

end subroutine solver_induction_hll
#endif
#endif
end module hydro_solvers
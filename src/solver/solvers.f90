module hydro_solvers
contains

!!!MHD wave fans may or may not include 1/4pi coefficients, depending on your setup and choice of units. Be careful!!!

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
    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,P_mag_lft,P_mag_rgt,mag_tension_y_lft,mag_tension_y_rgt,mag_tension_z_lft,mag_tension_z_rgt,mag_tension_x_lft,mag_tension_x_rgt,ca_lft,ca_rgt,c_fast_rgt,c_fast_lft


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

    magnetosonic_fast_rgt = dsqrt(csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/(4*pi*rho_rgt)) 
    magnetosonic_fast_lft = dsqrt(csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/(4*pi*rho_lft)) 

    lambda_llf_g = max(abs(u_lft)+magnetosonic_fast_lft,abs(u_rgt)+magnetosonic_fast_rgt)
#endif

#if NDUST>0
    if (ideal_MHD .or. dusty_nonideal_MHD) then
    !if (ideal_MHD) then !Recouple too gas even in presence of dust
                print*,"godunov ideal MHD"


    !If dust and ideal_MHD==true, field lines frozen to the gas and dust does not backreact on B.
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


        !!!This is the magnetodonic mode for theta_B = 0 --> reduces to a soundwave!!!
        !magnetosonic_fast_rgt = dsqrt(half*(csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt + dsqrt((csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt)**2-4*csr**2*(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt))) 
        !magnetosonic_fast_lft = dsqrt(half*(csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft + dsqrt((csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft)**2-4*csl**2*(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft))) 

        !!!This is the fast mode, obtained for theta_B = pi / 2. This is safe (although maybe diffusive) because it is the maximum speed possible!!
        ca_lft =dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft)
        ca_rgt =dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt)
        c_fast_lft = dsqrt(csl**2 + ca_lft**2) 
        c_fast_rgt = dsqrt(csr**2 + ca_rgt**2)

        lambda_llf_g = max(abs(u_lft)+c_fast_lft,abs(u_rgt)+c_fast_rgt)

    endif
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
    real(dp) :: flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt,ca_rgt,ca_lft,c_fast_lft,c_fast_rgt

#if MHD==1 
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

#if NDUST>0
    if (ideal_MHD .or. dusty_nonideal_MHD) then
    !if (ideal_MHD) then !Recouple too gas even in presence of dust    !If dust and ideal_MHD==true, field lines frozen to the gas and dust does not backreact on B.
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

    endif
#endif

#endif


    S_lft  = min(min(u_lft,u_rgt) -max(csl,csr),0.0d0)
    S_rgt  = max(max(u_lft,u_rgt) +max(csl,csr),0.0d0) 

#if MHD==1
#if NDUST==0

        !!!This is the magnetosonic mode for theta_B = 0 --> reduces to a soundwave!!!
        !magnetosonic_fast_rgt = dsqrt(half*(csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt + dsqrt((csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt)**2-4*csr**2*(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt))) 
        !magnetosonic_fast_lft = dsqrt(half*(csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft + dsqrt((csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft)**2-4*csl**2*(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft))) 

        !!!This is the fast mode, obtained for theta_B = pi / 2. This is safe (although maybe diffusive) because it is the maximum speed possible!!
        ca_lft =dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft)
        ca_rgt =dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt)
        c_fast_lft = dsqrt(csl**2 + ca_lft**2) 
        c_fast_rgt = dsqrt(csr**2 + ca_rgt**2)


        S_lft  = min(min(u_lft,u_rgt) -max(c_fast_lft,c_fast_rgt),0.0d0)
        S_rgt  = max(max(u_lft,u_rgt) +max(c_fast_lft,c_fast_rgt),0.0d0) 
#endif

#if NDUST>0
    if (ideal_MHD .or. dusty_nonideal_MHD) then
    !if (ideal_MHD) then !Recouple too gas even in presence of dust

        !!!This is the magnetosonic mode for theta_B = 0 --> reduces to a soundwave!!!
        !magnetosonic_fast_rgt = dsqrt(half*(csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt + dsqrt((csr**2+(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt)**2-4*csr**2*(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/rho_rgt))) 
        !magnetosonic_fast_lft = dsqrt(half*(csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft + dsqrt((csl**2+(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft)**2-4*csl**2*(Bx_lft**2+By_lft**2+Bz_lft**2)/rho_lft))) 

        !!!This is the fast mode, obtained for theta_B = pi / 2. This is safe (although maybe diffusive) because it is the maximum speed possible!!
        ca_lft =dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft)
        ca_rgt =dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt)
        c_fast_lft = dsqrt(csl**2 + ca_lft**2) 
        c_fast_rgt = dsqrt(csr**2 + ca_rgt**2)

        S_lft  = min(min(u_lft,u_rgt) -max(c_fast_lft,c_fast_rgt),0.0d0)
        S_rgt  = max(max(u_lft,u_rgt) +max(c_fast_lft,c_fast_rgt),0.0d0) 

    endif
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
    integer  :: idim,idust,i_u,i_v,i_rho,i_w,i,ipscal

    real(dp) :: S_lft,S_rgt,lambda_llf_d,csl,csr


    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt,flx_pscal_lft,flx_pscal_rgt
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
    
    if (ideal_MHD .eqv. .false.) then !if ideal_MHD.eqv. .true., dust is considered neutral
        flx_mom_u_rgt  = flx_mom_u_rgt + 1/(4*pi)*(P_mag_rgt + mag_tension_x_rgt)
        flx_mom_u_lft  = flx_mom_u_lft + 1/(4*pi)*(P_mag_lft + mag_tension_x_lft)

        flx_mom_v_rgt = flx_mom_v_rgt + 1/(4*pi)*mag_tension_y_rgt
        flx_mom_v_lft = flx_mom_v_lft + 1/(4*pi)*mag_tension_y_lft
           
        flx_mom_w_rgt = flx_mom_w_rgt + 1/(4*pi)*mag_tension_z_rgt
        flx_mom_w_lft = flx_mom_w_lft + 1/(4*pi)*mag_tension_z_lft
    endif


#endif

#if MHD==0

     lambda_llf_d        = max(abs(u_lft),abs(u_rgt))

#endif

#if MHD==1
    if (ideal_MHD .eqv. .false.) then    
     lambda_llf_d        = max(abs(u_lft)+dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft),abs(u_rgt)+dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt))
    endif

    if (ideal_MHD .eqv. .true.) then    
     lambda_llf_d        = max(abs(u_lft),abs(u_rgt))
    endif
#endif

    flx(i_rho)  =  0.d0
    flx(i_u)    =  0.d0
    flx(i_v)    =  0.d0
    flx(i_w)    =  0.d0


    flx(i_rho) = half*(flx_rho_lft  + flx_rho_rgt)   - half*lambda_llf_d*(rho_rgt - rho_lft)
    flx(i_u)  = half*(flx_mom_u_lft + flx_mom_u_rgt)  - half*lambda_llf_d*(mom_u_rgt - mom_u_lft)
    flx(i_v)  = half*(flx_mom_v_lft + flx_mom_v_rgt)  - half*lambda_llf_d*(mom_v_rgt - mom_v_lft)
    flx(i_w)  = half*(flx_mom_w_lft + flx_mom_w_rgt)  - half*lambda_llf_d*(mom_w_rgt - mom_w_lft)

#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
        flx_pscal_lft = qleft(idust_pscal(idust,ipscal)) * rho_lft  * u_lft
        flx_pscal_rgt = qright(idust_pscal(idust,ipscal)) * rho_rgt  * u_rgt
        flx(idust_pscal(idust,ipscal))  = half*(flx_pscal_lft  + flx_pscal_rgt)   - half*lambda_llf_d*(rho_rgt*qright(idust_pscal(idust,ipscal)) - rho_lft*qleft(idust_pscal(idust,ipscal)))
    end do
#endif 

end do
end subroutine solver_dust_llf

#endif
#endif


#if NDUST>0
#if SOLVERDUST==2

subroutine solver_dust_hll(qleft,qright,csl,csr,flx,idim,i)

    use parameters
    use commons
    use slope_limiter


    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    integer  :: idim,idust,i_u,i_v,i_rho,i_w,i,il,ir,ix,iy,icell,ixx,iyy,ipscal
    real(dp) :: csl,csr,P_lft,P_rgt
    real(dp) :: S_lft,S_rgt,lambda_llf_d
    real(dp) :: ca_lft,ca_rgt,cw_rgt,cw_lft,magnetosonic_fast_rgt,magnetosonic_fast_lft



    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt,flx_pscal_lft,flx_pscal_rgt
#if MHD==1

    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,P_mag_lft,P_mag_rgt,mag_tension_y_lft,mag_tension_y_rgt,mag_tension_z_lft,mag_tension_z_rgt,mag_tension_x_lft,mag_tension_x_rgt,c_fast_lft,c_fast_rgt 


    real(dp) ::flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt,deta_Hall,deta_Hall_l,eta_Hall_y_left,eta_Hall_y_right
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

#if DUST_PRESSURE==1

        P_rgt       = qright(iPd(idust))
        P_lft       = qleft(iPd(idust))

        flx_mom_u_rgt  = rho_rgt * u_rgt**2 + P_rgt
        flx_mom_u_lft  = rho_lft * u_lft**2 + P_lft

#endif
        flx_mom_v_rgt = rho_rgt * u_rgt  * v_rgt
        flx_mom_v_lft = rho_lft * u_lft  * v_lft
           
        flx_mom_w_rgt = rho_rgt * u_rgt  * w_rgt
        flx_mom_w_lft = rho_lft * u_lft  * w_lft

#if MHD==1

    if (dusty_nonideal_MHD_no_electron) then
        if (ideal_MHD .eqv. .false.) then
        if (idust==i_coupled_species) then

            flx_rho_rgt   = rho_rgt  * u_rgt
            flx_rho_lft   = rho_lft  * u_lft


            flx_mom_u_rgt  = flx_mom_u_rgt + 1/(4*pi)*(P_mag_rgt + mag_tension_x_rgt)
            flx_mom_u_lft  = flx_mom_u_lft + 1/(4*pi)*(P_mag_lft + mag_tension_x_lft)


            flx_mom_v_rgt = flx_mom_v_rgt + 1/(4*pi)*mag_tension_y_rgt
            flx_mom_v_lft = flx_mom_v_lft + 1/(4*pi)*mag_tension_y_lft
               
            flx_mom_w_rgt = flx_mom_w_rgt + 1/(4*pi)*mag_tension_z_rgt
            flx_mom_w_lft = flx_mom_w_lft + 1/(4*pi)*mag_tension_z_lft

        endif
        endif
    endif
#endif





#if MHD==0
    S_rgt  = max(max(u_lft,u_rgt),0.0d0) 
    S_lft  = min(min(u_lft,u_rgt),0.0d0)

#if DUST_PRESSURE==1
        S_lft  = min(min(u_lft,u_rgt)-max(delta_dust_cs*csl,delta_dust_cs*csr),0.0d0)
        S_rgt  = max(max(u_lft,u_rgt)+max(delta_dust_cs*csl,delta_dust_cs*csr),0.0d0)
#endif


    if (u_lft/=u_rgt) then

    flx(i_rho)            = (S_rgt*flx_rho_lft  -S_lft*flx_rho_rgt  + S_rgt*S_lft*(rho_rgt-rho_lft))      / (S_rgt-S_lft)
    flx(i_u)  = (S_rgt*flx_mom_u_lft-S_lft*flx_mom_u_rgt+ S_rgt*S_lft*(mom_u_rgt-mom_u_lft))  / (S_rgt-S_lft)
    flx(i_v)  = (S_rgt*flx_mom_v_lft-S_lft*flx_mom_v_rgt+ S_rgt*S_lft*(mom_v_rgt-mom_v_lft))  / (S_rgt-S_lft)
    flx(i_w)             = (S_rgt*flx_mom_w_lft-S_lft*flx_mom_w_rgt+ S_rgt*S_lft*(mom_w_rgt-mom_w_lft))  / (S_rgt-S_lft)



#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal

        flx_pscal_lft = qleft(idust_pscal(idust,ipscal)) * rho_lft  * u_lft
        flx_pscal_rgt = qright(idust_pscal(idust,ipscal)) * rho_rgt  * u_rgt

        flx(idust_pscal(idust,ipscal))  = (S_rgt*flx_pscal_lft  - S_lft*flx_pscal_rgt + S_rgt*S_lft*(rho_rgt*qright(idust_pscal(idust,ipscal)) - rho_lft*qleft(idust_pscal(idust,ipscal)))) / (S_rgt-S_lft)
    end do
#endif 


    end if

    if (u_lft==u_rgt) then !Switch back to llf

    lambda_llf_d        = max(abs(u_lft),abs(u_rgt))

    flx(i_rho) = half*(flx_rho_lft   + flx_rho_rgt)    - half*lambda_llf_d*(rho_rgt   - rho_lft)
    flx(i_u)  = half*(flx_mom_u_lft + flx_mom_u_rgt)  - half*lambda_llf_d*(mom_u_rgt - mom_u_lft)
    flx(i_v)  = half*(flx_mom_v_lft + flx_mom_v_rgt)  - half*lambda_llf_d*(mom_v_rgt - mom_v_lft)
    flx(i_w)  = half*(flx_mom_w_lft + flx_mom_w_rgt)  - half*lambda_llf_d*(mom_w_rgt - mom_w_lft)

#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
        flx_pscal_lft = qleft(idust_pscal(idust,ipscal)) * rho_lft  * u_lft
        flx_pscal_rgt = qright(idust_pscal(idust,ipscal)) * rho_rgt  * u_rgt
        flx(idust_pscal(idust,ipscal))  = half*(flx_pscal_lft  + flx_pscal_rgt)   - half*lambda_llf_d*(rho_rgt*qright(idust_pscal(idust,ipscal)) - rho_lft*qleft(idust_pscal(idust,ipscal)))
    end do
#endif 



    end if
#endif

#if MHD==1

    if (ideal_MHD .eqv. .false.) then

    !----------------------------------------------------------------------------------------------------------
    !For both dusty_nonideal_MHD and dusty_nonideal_MHD_no_electron, we expect dusty Alfvén waves to propagate. 
    !----------------------------------------------------------------------------------------------------------


    ! ca_lft =dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*qleft(irhod(1))) !!Here 4pi for the setup used in Vallucci-Goy+25
    ! ca_rgt =dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*qright(irhod(1)))

    ca_lft =dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft) !!Here 4pi for the setup used in Vallucci-Goy+25
    ca_rgt =dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt)


    S_rgt  = max(max(u_lft,u_rgt) +max(ca_lft,ca_rgt),0.0d0) 
    S_lft  = min(min(u_lft,u_rgt) -max(ca_lft,ca_rgt),0.0d0)
    ! S_rgt  = max(max(qleft(index_vdn(1,idim)),qright(index_vdn(1,idim))) +max(ca_lft,ca_rgt),0.0d0) 
    ! S_lft  = min(min(qleft(index_vdn(1,idim)),qright(index_vdn(1,idim))) -max(ca_lft,ca_rgt),0.0d0)



#if DUST_PRESSURE==1 
!!!TO DO!!! Add a NDUST==1. If NDUST>1, we don't consider any dust pressure.
        c_fast_lft = dsqrt((delta_dust_cs*csl)**2 + ca_lft**2) !Safer to use this
        c_fast_rgt = dsqrt((delta_dust_cs*csr)**2 + ca_rgt**2)

        S_lft  = min(min(u_lft,u_rgt) -max(c_fast_lft,c_fast_rgt),0.0d0)
        S_rgt  = max(max(u_lft,u_rgt) +max(c_fast_lft,c_fast_rgt),0.0d0) 
#endif

    endif

  if (ideal_MHD .eqv. .true.) then !Dust neutral, gas is in ideal MHD. Here, we use separate wavefans

    S_rgt  = max(max(u_lft,u_rgt),0.0d0) 
    S_lft  = min(min(u_lft,u_rgt),0.0d0)


#if DUST_PRESSURE==1 
!!!TO DO!!! Add a NDUST==1. If NDUST>1, we don't consider any dust pressure.
        c_lft = delta_dust_cs*csl 
        c_rgt = delta_dust_cs*csr

        !magnetosonic_fast_rgt = dsqrt(half*(c_fast_rgt**2 + dsqrt(c_fast_rgt**4-4*(delta_dust_cs*csr)**2*ca_rgt**2))) !In 1D along B: reduces to a simple soundwave
        !magnetosonic_fast_lft = dsqrt(half*(c_fast_lft**2 + dsqrt(c_fast_lft**4-4*(delta_dust_cs*csl)**2*ca_lft**2))) 

        S_lft  = min(min(u_lft,u_rgt) -max(c_lft,c_rgt),0.0d0)
        S_rgt  = max(max(u_lft,u_rgt) +max(c_lft,c_rgt),0.0d0) 
#endif

    endif
    



    if (S_rgt/=S_lft) then

        flx(i_rho)            = (S_rgt*flx_rho_lft  -S_lft*flx_rho_rgt  + S_rgt*S_lft*(rho_rgt-rho_lft))      / (S_rgt-S_lft)
        flx(i_u)  = (S_rgt*flx_mom_u_lft-S_lft*flx_mom_u_rgt+ S_rgt*S_lft*(mom_u_rgt-mom_u_lft))  / (S_rgt-S_lft)
        flx(i_v)  = (S_rgt*flx_mom_v_lft-S_lft*flx_mom_v_rgt+ S_rgt*S_lft*(mom_v_rgt-mom_v_lft))  / (S_rgt-S_lft)
        flx(i_w)             = (S_rgt*flx_mom_w_lft-S_lft*flx_mom_w_rgt+ S_rgt*S_lft*(mom_w_rgt-mom_w_lft))  / (S_rgt-S_lft)


#if NDUSTPSCAL>0
        do ipscal=1,ndustpscal
            flx_pscal_lft = qleft(idust_pscal(idust,ipscal)) * rho_lft  * u_lft
            flx_pscal_rgt = qright(idust_pscal(idust,ipscal)) * rho_rgt  * u_rgt
            flx(idust_pscal(idust,ipscal))  = (S_rgt*flx_pscal_lft  - S_lft*flx_pscal_rgt + S_rgt*S_lft*(rho_rgt*qright(idust_pscal(idust,ipscal)) - rho_lft*qleft(idust_pscal(idust,ipscal)))) / (S_rgt-S_lft)
        end do
#endif 

    endif

    if (S_rgt==S_lft) then

        if (ideal_MHD .eqv. .false.) then

        lambda_llf_d        = max(abs(ca_lft),abs(ca_rgt))
#if DUST_PRESSURE==1
         lambda_llf_d        = max(abs(c_fast_lft),abs(c_fast_rgt))
#endif

        endif



        if (ideal_MHD .eqv. .true.) then

            lambda_llf_d        = max(abs(u_lft),abs(u_rgt))
#if DUST_PRESSURE==1
         lambda_llf_d        = max(abs(c_lft),abs(c_rgt))
#endif

        endif



    flx(i_rho) = half*(flx_rho_lft   + flx_rho_rgt)    - half*lambda_llf_d*(rho_rgt   - rho_lft)
    flx(i_u)  = half*(flx_mom_u_lft + flx_mom_u_rgt)  - half*lambda_llf_d*(mom_u_rgt - mom_u_lft)
    flx(i_v)  = half*(flx_mom_v_lft + flx_mom_v_rgt)  - half*lambda_llf_d*(mom_v_rgt - mom_v_lft)
    flx(i_w)  = half*(flx_mom_w_lft + flx_mom_w_rgt)  - half*lambda_llf_d*(mom_w_rgt - mom_w_lft)

#if NDUSTPSCAL>0
    do ipscal=1,ndustpscal
        flx_pscal_lft = qleft(idust_pscal(idust,ipscal)) * rho_lft  * u_lft
        flx_pscal_rgt = qright(idust_pscal(idust,ipscal)) * rho_rgt  * u_rgt
        flx(idust_pscal(idust,ipscal))  = half*(flx_pscal_lft  + flx_pscal_rgt)   - half*lambda_llf_d*(rho_rgt*qright(idust_pscal(idust,ipscal)) - rho_lft*qleft(idust_pscal(idust,ipscal)))
    end do
#endif 

    endif

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
    real(dp) :: ca_lft,ca_rgt,magnetosonic_fast_lft,magnetosonic_fast_rgt,csl,csr,cw_lft,cw_rgt,dJy,dJz,dJy_l,dJz_l,Jy_left,Jy_right,Jz_left,Jz_right,c_fast_lft,c_fast_rgt

    real(dp) :: deta_o_il,deta_h_il,deta_a_il,dzd,dzd_il,zd_left,zd_right,deta_Hall,deta_Hall_l,eta_Hall_y_left,eta_Hall_y_right,eta_Hall_z_left,eta_Hall_z_right

    real(dp) :: nd_zd_over_ni,nd_zd_over_ni_left,nd_zd_over_ni_right,dne,dne_il,ne_left,ne_right,dni,dni_il,ni_left,ni_right,dhall_i,dhall_il,hall_i_left,hall_i_right
    real(dp) :: flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt,total_dust_current_z_lft,total_dust_current_y_lft,total_dust_current_x_lft,total_dust_current_z_rgt,total_dust_current_y_rgt,total_dust_current_x_rgt,B_norm_lft,B_norm_rgt
    real(dp) :: deta_o,deta_o_l,eta_o_left,eta_o_right,dJdx_tot,dJdx_tot_l,Jdx_tot_left,Jdx_tot_right,dJdy_tot,dJdy_tot_l,Jdy_tot_left,Jdy_tot_right,dJdz_tot,dJdz_tot_l,Jdz_tot_left,Jdz_tot_right
    real(dp) :: deta_H,deta_H_l,eta_H_left,eta_H_right,deta_a,deta_a_l,eta_a_left,eta_a_right,db_unit_x,db_unit_x_l,db_unit_y,db_unit_y_l,db_unit_z,db_unit_z_l,b_unit_x_left,b_unit_x_right,b_unit_y_left,b_unit_y_right,b_unit_z_left,b_unit_z_right

    Bx_lft   = qleft(iBx)
    Bx_rgt   = qright(iBx)
    By_lft   = qleft(iBy)
    By_rgt   = qright(iBy)
    Bz_lft   = qleft(iBz)
    Bz_rgt   = qright(iBz)




#if NDUST>0

    if (ideal_MHD .eqv. .false.) then
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
    endif

    if (ideal_MHD .eqv. .true.) then !There is dust, but considered neutral. Ideal MHD on the gas.

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

    endif



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


      if (dusty_nonideal_MHD_no_electron .eqv. .true.) then !Additional terms in the fluxes for the induction equation
        if (ideal_MHD .eqv. .false.) then
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

                deta_Hall = slope_limit(2.0d0*(abs(eta_eff_Hall_y(i)) - abs(eta_eff_Hall_y(il)))/(dx(i,1)+dx(il,1)),2.0d0*(abs(eta_eff_Hall_y(ir)) - abs(eta_eff_Hall_y(i)))/(dx(ir,1)+dx(i,1)))
                deta_Hall_l = slope_limit(2.0d0*(abs(eta_eff_Hall_y(i-1)) - abs(eta_eff_Hall_y(il-1)))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(abs(eta_eff_Hall_y(ir-1)) - abs(eta_eff_Hall_y(i-1)))/(dx(ir-1,1)+dx(i-1,1)))

                eta_Hall_y_left = abs(eta_eff_Hall_y(il)) + half*deta_Hall_l*dx(il,1)
                eta_Hall_y_right = abs(eta_eff_Hall_y(i))- half*deta_Hall*dx(i,1)

                eta_Hall_z_left = - eta_Hall_y_left
                eta_Hall_z_right = - eta_Hall_y_right

                dJy = slope_limit(2.0d0*(Jy(i) - Jy(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jy(ir) - Jy(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJy_l = slope_limit(2.0d0*(Jy(i-1) - Jy(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jy(ir-1) - Jy(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total current

                dJz = slope_limit(2.0d0*(Jz(i) - Jz(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jz(ir) - Jz(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJz_l = slope_limit(2.0d0*(Jz(i-1) - Jz(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jz(ir-1) - Jz(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total current

                Jy_left = Jy(il) + half*dJy_l*dx(il,1)
                Jy_right = Jy(i) - half*dJy*dx(i,1)

                Jz_left = Jz(il) + half*dJz_l*dx(il,1)
                Jz_right = Jz(i) - half*dJz*dx(i,1)




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

        if (Hall_effect) then
            !Hall term
            flx_By_lft = flx_By_lft + eta_Hall_y_left*Jy_left
            flx_By_rgt = flx_By_rgt + eta_Hall_y_right*Jy_right

            flx_Bz_lft = flx_Bz_lft - eta_Hall_z_left*Jz_left
            flx_Bz_rgt = flx_Bz_rgt - eta_Hall_z_right*Jz_right

        endif

        endif
    endif



      if (dusty_nonideal_MHD .eqv. .true.) then !Additional terms in the fluxes for the induction equation
        if (ideal_MHD .eqv. .false.) then

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

        ix=ixx(i)
        iy=iyy(i)

            if(slope_type>0) then
                il = icell(ix-1,iy)
                ir = icell(ix+1,iy)

                deta_o = slope_limit(2.0d0*(eta_o(i) - eta_o(il))/(dx(i,1)+dx(il,1)),2.0d0*(eta_o(ir) - eta_o(i))/(dx(ir,1)+dx(i,1)))
                deta_o_l = slope_limit(2.0d0*(eta_o(i-1) - eta_o(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(eta_o(ir-1) - eta_o(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                eta_o_left = eta_o(il) + half*deta_o_l*dx(il,1)
                eta_o_right = eta_o(i) - half*deta_o*dx(i,1)

                deta_H = slope_limit(2.0d0*(eta_H(i) - eta_H(il))/(dx(i,1)+dx(il,1)),2.0d0*(eta_H(ir) - eta_H(i))/(dx(ir,1)+dx(i,1)))
                deta_H_l = slope_limit(2.0d0*(eta_H(i-1) - eta_H(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(eta_H(ir-1) - eta_H(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                eta_H_left = eta_H(il) + half*deta_H_l*dx(il,1)
                eta_H_right = eta_H(i) - half*deta_H*dx(i,1)

                deta_a = slope_limit(2.0d0*(eta_a(i) - eta_a(il))/(dx(i,1)+dx(il,1)),2.0d0*(eta_a(ir) - eta_a(i))/(dx(ir,1)+dx(i,1)))
                deta_a_l = slope_limit(2.0d0*(eta_a(i-1) - eta_a(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(eta_a(ir-1) - eta_a(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                eta_a_left = eta_a(il) + half*deta_a_l*dx(il,1)
                eta_a_right = eta_a(i) - half*deta_a*dx(i,1)



                dJy = slope_limit(2.0d0*(Jy(i) - Jy(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jy(ir) - Jy(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJy_l = slope_limit(2.0d0*(Jy(i-1) - Jy(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jy(ir-1) - Jy(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total current

                dJz = slope_limit(2.0d0*(Jz(i) - Jz(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jz(ir) - Jz(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJz_l = slope_limit(2.0d0*(Jz(i-1) - Jz(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jz(ir-1) - Jz(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total current

                Jy_left = Jy(il) + half*dJy_l*dx(il,1)
                Jy_right = Jy(i) - half*dJy*dx(i,1)

                Jz_left = Jz(il) + half*dJz_l*dx(il,1)
                Jz_right = Jz(i) - half*dJz*dx(i,1)



                dJdx_tot = slope_limit(2.0d0*(Jdx_tot(i) - Jdx_tot(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jdx_tot(ir) - Jdx_tot(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJdx_tot_l = slope_limit(2.0d0*(Jdx_tot(i-1) - Jdx_tot(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jdx_tot(ir-1) - Jdx_tot(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total dust current

                dJdy_tot = slope_limit(2.0d0*(Jdy_tot(i) - Jdy_tot(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jdy_tot(ir) - Jdy_tot(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJdy_tot_l = slope_limit(2.0d0*(Jdy_tot(i-1) - Jdy_tot(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jdy_tot(ir-1) - Jdy_tot(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total dust current

                dJdz_tot = slope_limit(2.0d0*(Jdz_tot(i) - Jdz_tot(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jdz_tot(ir) - Jdz_tot(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJdz_tot_l = slope_limit(2.0d0*(Jdz_tot(i-1) - Jdz_tot(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jdz_tot(ir-1) - Jdz_tot(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total dust current


                Jdx_tot_left = Jdx_tot(il) + half*dJdx_tot_l*dx(il,1)
                Jdx_tot_right = Jdx_tot(i) - half*dJdx_tot*dx(i,1)

                Jdy_tot_left = Jdy_tot(il) + half*dJdy_tot_l*dx(il,1)
                Jdy_tot_right = Jdy_tot(i) - half*dJdy_tot*dx(i,1)

                Jdz_tot_left = Jdz_tot(il) + half*dJdz_tot_l*dx(il,1)
                Jdz_tot_right = Jdz_tot(i) - half*dJdz_tot*dx(i,1)


                db_unit_x = slope_limit(2.0d0*(b_unit_x(i) - b_unit_x(il))/(dx(i,1)+dx(il,1)),2.0d0*(b_unit_x(ir) - b_unit_x(i))/(dx(ir,1)+dx(i,1))) 
                db_unit_y = slope_limit(2.0d0*(b_unit_y(i) - b_unit_y(il))/(dx(i,1)+dx(il,1)),2.0d0*(b_unit_y(ir) - b_unit_y(i))/(dx(ir,1)+dx(i,1))) 
                db_unit_z = slope_limit(2.0d0*(b_unit_z(i) - b_unit_z(il))/(dx(i,1)+dx(il,1)),2.0d0*(b_unit_z(ir) - b_unit_z(i))/(dx(ir,1)+dx(i,1))) 

                db_unit_x_l = slope_limit(2.0d0*(b_unit_x(i-1) - b_unit_x(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(b_unit_x(ir-1) - b_unit_x(i-1))/(dx(ir-1,1)+dx(i-1,1))) 
                db_unit_y_l = slope_limit(2.0d0*(b_unit_y(i-1) - b_unit_y(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(b_unit_y(ir-1) - b_unit_y(i-1))/(dx(ir-1,1)+dx(i-1,1))) 
                db_unit_z_l = slope_limit(2.0d0*(b_unit_z(i-1) - b_unit_z(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(b_unit_z(ir-1) - b_unit_z(i-1))/(dx(ir-1,1)+dx(i-1,1))) 

                b_unit_x_left = b_unit_x(il) + half*db_unit_x_l*dx(il,1)
                b_unit_x_right = b_unit_x(i) - half*db_unit_x*dx(i,1)


                b_unit_y_left = b_unit_y(il) + half*db_unit_y_l*dx(il,1)
                b_unit_y_right = b_unit_y(i) - half*db_unit_y*dx(i,1)

                b_unit_z_left = b_unit_z(il) + half*db_unit_z_l*dx(il,1)
                b_unit_z_right = b_unit_z(i) - half*db_unit_z*dx(i,1)



            endif


        B_norm_lft = dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)
        B_norm_rgt = dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2) 


        !Beware: here signs are reversed with respect to predictor step (fluxes are defined in the left-hand side of the equation) 


        flx_By_lft = -Bx_lft*v_lft + By_lft*u_lft !In this model, this term involves the gas velocity
        flx_By_rgt = -Bx_rgt*v_rgt + By_rgt*u_rgt

        flx_Bz_lft = -Bx_lft*w_lft + Bz_lft*u_lft
        flx_Bz_rgt = -Bx_rgt*w_rgt + Bz_rgt*u_rgt


        !Additional term
        !Ohm
        flx_By_lft = flx_By_lft + clight*eta_o_left*Jdz_tot_left
        flx_By_rgt = flx_By_rgt + clight*eta_o_right*Jdz_tot_right

        flx_Bz_lft = flx_Bz_lft - clight*eta_o_left*Jdy_tot_left
        flx_Bz_rgt = flx_Bz_rgt - clight*eta_o_right*Jdy_tot_right

        !AD
        flx_By_lft = flx_By_lft + clight*eta_a_left*( (b_unit_x_left**2+b_unit_y_left**2)*Jdz_tot_left - b_unit_x_left*b_unit_z_left*Jdx_tot_left - b_unit_y_left*b_unit_z_left*Jdy_tot_left) - clight**2/(4*pi) * eta_a_left*b_unit_y_left*b_unit_z_left*(-Jy_left)
        flx_By_rgt = flx_By_rgt + clight*eta_a_right*( (b_unit_x_right**2+b_unit_y_right**2)*Jdz_tot_right - b_unit_x_right*b_unit_z_right*Jdx_tot_right - b_unit_y_right*b_unit_z_right*Jdy_tot_right) - clight**2/(4*pi) * eta_a_right*b_unit_y_right*b_unit_z_right*(-Jy_right)

        flx_Bz_lft = flx_Bz_lft - clight*eta_a_left*((b_unit_x_left**2+b_unit_z_left**2)*Jdy_tot_left - b_unit_y_left*b_unit_z_left*Jdz_tot_left - b_unit_x_left*b_unit_y_left*Jdx_tot_left) - clight**2/(4*pi) * eta_a_left*b_unit_y_left*b_unit_z_left*(Jz_left)
        flx_Bz_rgt = flx_Bz_rgt - clight*eta_a_right*((b_unit_x_right**2+b_unit_z_right**2)*Jdy_tot_right - b_unit_y_right*b_unit_z_right*Jdz_tot_right - b_unit_x_right*b_unit_y_right*Jdx_tot_right) - clight**2/(4*pi) * eta_a_right*b_unit_y_right*b_unit_z_right*(Jz_right)

       

        if (Hall_effect) then
            !Hall term
            flx_By_lft = flx_By_lft + clight*eta_H_left*b_unit_y_left*Jdx_tot_left + clight**2/(4*pi) * eta_H_left * b_unit_x_left * Jy_left - clight*eta_H_left * b_unit_x_left * Jdy_tot_left
            flx_By_rgt = flx_By_rgt + clight*eta_H_right*b_unit_y_right*Jdx_tot_right + clight**2/(4*pi) * eta_H_right * b_unit_x_right * Jy_right - clight*eta_H_right * b_unit_x_right * Jdy_tot_right

            flx_Bz_lft = flx_Bz_lft + clight*eta_H_left*b_unit_z_left*Jdx_tot_left + clight**2/(4*pi) * eta_H_left * b_unit_x_left * Jz_left - clight*eta_H_left*b_unit_x_left*Jdz_tot_left
            flx_Bz_rgt = flx_Bz_rgt + clight*eta_H_right*b_unit_z_right*Jdx_tot_right + clight**2/(4*pi) * eta_H_right * b_unit_x_right * Jz_right - clight*eta_H_right*b_unit_x_right*Jdz_tot_right

        endif

        endif
    endif

#endif



!HLL

#if NDUST==0

    ca_lft =dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft)
    ca_rgt =dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt)

    c_fast_lft = dsqrt((csl)**2 + ca_lft**2)
    c_fast_rgt = dsqrt((csr)**2 + ca_rgt**2) 


    S_rgt  = max(max(u_lft,u_rgt) +max(c_fast_lft,c_fast_rgt),0.0d0) 
    S_lft  = min(min(u_lft,u_rgt) -max(c_fast_lft,c_fast_rgt),0.0d0)

    if (S_rgt/=S_lft) then

        flx(iBy)            = (S_rgt*flx_By_lft  -S_lft*flx_By_rgt  + S_rgt*S_lft*(By_rgt-By_lft))      / (S_rgt-S_lft)
        flx(iBz)            = (S_rgt*flx_Bz_lft  -S_lft*flx_Bz_rgt  + S_rgt*S_lft*(Bz_rgt-Bz_lft))      / (S_rgt-S_lft)

    endif

    if (S_lft==S_rgt) then

        lambda_llf_B        = max(abs(u_lft)+c_fast_lft,abs(u_rgt)+c_fast_rgt)


        flx(iBy) = half*(flx_By_lft   + flx_By_rgt)    - half*lambda_llf_B*(By_rgt   - By_lft)
        flx(iBz) = half*(flx_Bz_lft   + flx_Bz_rgt)    - half*lambda_llf_B*(Bz_rgt   - Bz_lft)

    endif

#endif

#if NDUST>0

    if (ideal_MHD .eqv. .false.) then

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

        S_rgt  = max(max(u_lft,u_rgt) +max(ca_lft,ca_rgt),0.0d0) 
        S_lft  = min(min(u_lft,u_rgt) -max(ca_lft,ca_rgt),0.0d0)



#if DUST_PRESSURE==1

        c_fast_lft = dsqrt((delta_dust_cs*csl)**2 + ca_lft**2)
        c_fast_rgt = dsqrt((delta_dust_cs*csr)**2 + ca_rgt**2)

        ! magnetosonic_fast_rgt = dsqrt(half*(c_fast_rgt**2 + dsqrt(c_fast_rgt**4-4*(delta_dust_cs*csr)**2*ca_rgt**2))) !In 1D along B: reduces to a simple soundwave
        ! magnetosonic_fast_lft = dsqrt(half*(c_fast_lft**2 + dsqrt(c_fast_lft**4-4*(delta_dust_cs*csl)**2*ca_lft**2))) 
        S_lft  = min(min(u_lft,u_rgt) -max(c_fast_lft,c_fast_rgt),0.0d0)
        S_rgt  = max(max(u_lft,u_rgt) +max(c_fast_lft,c_fast_rgt),0.0d0) 
#endif


        if (Hall_effect .and. dusty_nonideal_MHD) then
        !Hall effect introduces new waves

            cw_lft = eta_H_left*pi/(2*dx(i,1)) + dsqrt((eta_H_left*pi/(2*dx(i,1)))**2 + ca_lft**2) !Whistler wave
            cw_rgt = eta_H_right*pi/(2*dx(i+1,1)) + dsqrt((eta_H_right*pi/(2*dx(i+1,1)))**2 + ca_rgt**2) !Whistler wave 

            S_rgt  = max(max(u_lft,u_rgt) +max(cw_lft,cw_rgt),0.0d0) 
            S_lft  = min(min(u_lft,u_rgt) -max(cw_lft,cw_rgt),0.0d0)

        endif

        if (Hall_effect .and. dusty_nonideal_MHD_no_electron) then
        !Hall effect introduces new waves

            cw_lft = eta_Hall_y_left*pi/(2*dx(i,1)) + dsqrt((eta_Hall_y_left*pi/(2*dx(i,1)))**2 + ca_lft**2) !Whistler wave
            cw_rgt = eta_Hall_y_right*pi/(2*dx(i+1,1)) + dsqrt((eta_Hall_y_right*pi/(2*dx(i+1,1)))**2 + ca_rgt**2) !Whistler wave 

            S_rgt  = max(max(u_lft,u_rgt) +max(cw_lft,cw_rgt),0.0d0) 
            S_lft  = min(min(u_lft,u_rgt) -max(cw_lft,cw_rgt),0.0d0)

        endif


        if (S_lft/=S_rgt) then

            flx(iBy)            = (S_rgt*flx_By_lft  -S_lft*flx_By_rgt  + S_rgt*S_lft*(By_rgt-By_lft))      / (S_rgt-S_lft)
            flx(iBz)            = (S_rgt*flx_Bz_lft  -S_lft*flx_Bz_rgt  + S_rgt*S_lft*(Bz_rgt-Bz_lft))      / (S_rgt-S_lft)

        endif


        if (S_lft==S_rgt) then


            lambda_llf_B        = max(abs(u_lft) + abs(ca_lft),abs(u_rgt) + abs(ca_rgt))

#if DUST_PRESSURE==1

            lambda_llf_B        = max(abs(u_lft) + abs(magnetosonic_fast_lft),abs(u_rgt) + abs(magnetosonic_fast_rgt)) 
#endif



            if ((Hall_effect .eqv. .true.) .and. (dusty_nonideal_MHD_no_electron .eqv. .true.) .or. (dusty_nonideal_MHD .eqv. .true.)) lambda_llf_B = max(abs(u_lft) + abs(cw_lft),abs(u_rgt) + abs(cw_rgt)) !In this particular case, only variable B includes whistler speed in the wafe fan. The B solver no longer coincides with the dust solver. No need to add the dust velocity in lambda_llf_B.


            flx(iBy) = half*(flx_By_lft   + flx_By_rgt)    - half*lambda_llf_B*(By_rgt   - By_lft)
            flx(iBz) = half*(flx_Bz_lft   + flx_Bz_rgt)    - half*lambda_llf_B*(Bz_rgt   - Bz_lft)

        endif


    endif


    if (ideal_MHD .eqv. .true.) then !Couple B to the gas


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

        flx_Bx_lft = 0.0d0
        flx_Bx_rgt = 0.0d0
       
        flx_By_lft = By_lft*u_lft - Bx_lft*v_lft
        flx_By_rgt = By_rgt*u_rgt - Bx_rgt*v_rgt

        flx_Bz_lft = Bz_lft*u_lft - Bx_lft*w_lft
        flx_Bz_rgt = Bz_rgt*u_rgt - Bx_rgt*w_rgt    




        flx(iBx)    = 0.d0

        ca_lft =dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft)
        ca_rgt =dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt)

        c_fast_lft = dsqrt((csl)**2 + ca_lft**2)
        c_fast_rgt = dsqrt((csr)**2 + ca_rgt**2)

        ! magnetosonic_fast_rgt = dsqrt(half*(c_fast_rgt**2 + dsqrt(c_fast_rgt**4-4*(csr)**2*ca_rgt**2))) !In 1D along B: reduces to a simple soundwave
        ! magnetosonic_fast_lft = dsqrt(half*(c_fast_lft**2 + dsqrt(c_fast_lft**4-4*(csl)**2*ca_lft**2))) 

        S_lft  = min(min(u_lft,u_rgt) -max(c_fast_lft,c_fast_rgt),0.0d0)
        S_rgt  = max(max(u_lft,u_rgt) +max(c_fast_lft,c_fast_rgt),0.0d0) 


        if (S_lft/=S_rgt) then

            flx(iBy)            = (S_rgt*flx_By_lft  -S_lft*flx_By_rgt  + S_rgt*S_lft*(By_rgt-By_lft))      / (S_rgt-S_lft)
            flx(iBz)            = (S_rgt*flx_Bz_lft  -S_lft*flx_Bz_rgt  + S_rgt*S_lft*(Bz_rgt-Bz_lft))      / (S_rgt-S_lft)

        endif


        if (S_lft==S_rgt) then


            lambda_llf_B      = max(abs(u_lft)+abs(c_fast_lft),abs(u_rgt)+abs(c_fast_rgt))

            flx(iBy) = half*(flx_By_lft   + flx_By_rgt)    - half*lambda_llf_B*(By_rgt   - By_lft)
            flx(iBz) = half*(flx_Bz_lft   + flx_Bz_rgt)    - half*lambda_llf_B*(Bz_rgt   - Bz_lft)


        endif


    endif 


#endif 


end subroutine solver_induction_hll
#endif
#endif


! #if MHD==1
! #if NDUST==1
! #if SOLVERB==3

! subroutine solver_Hall_hll(qleft,qright,flx,csl,csr,idim,i)
!     use parameters
!     use commons
!     use slope_limiter

!     implicit none

!     real(dp),dimension(1:nvar),intent(in) :: qright,qleft
!     real(dp),dimension(1:nvar),intent(inout) :: flx
!     integer  :: idim,idust,i_u,i_v,i_rho,i_w,i,il,ir,ix,iy,icell,ixx,iyy

!     real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,lambda_llf_B
!     real(dp) :: u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt,rho_rgt,rho_lft
!     real(dp) :: S_lft,S_rgt
!     real(dp) :: ca_lft,ca_rgt,magnetosonic_fast_lft,magnetosonic_fast_rgt,csl,csr,cw_lft,cw_rgt,dJy,dJz,dJy_l,dJz_l,Jy_left,Jy_right,Jz_left,Jz_right

!     real(dp) :: deta_o,deta_h,deta_a,deta_o_il,deta_h_il,deta_a_il,eta_o_left,eta_h_left,eta_a_left,eta_o_right,eta_h_right,eta_a_right,dzd,dzd_il,zd_left,zd_right,deta_Hall,deta_Hall_l,eta_Hall_y_left,eta_Hall_y_right,eta_Hall_z_left,eta_Hall_z_right

!     real(dp) :: nd_zd_over_ni,nd_zd_over_ni_left,nd_zd_over_ni_right,dne,dne_il,ne_left,ne_right,dni,dni_il,ni_left,ni_right,dhall_i,dhall_il,hall_i_left,hall_i_right
!     real(dp) :: flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt,total_dust_current_z_lft,total_dust_current_y_lft,total_dust_current_x_lft,total_dust_current_z_rgt,total_dust_current_y_rgt,total_dust_current_x_rgt,B_norm_lft,B_norm_rgt


!     Bx_lft   = qleft(iBx)
!     Bx_rgt   = qright(iBx)
!     By_lft   = qleft(iBy)
!     By_rgt   = qright(iBy)
!     Bz_lft   = qleft(iBz)
!     Bz_rgt   = qright(iBz)





!     flx_Bx_lft = 0.0d0
!     flx_Bx_rgt = 0.0d0
   
!     flx_By_lft = 0.0d0
!     flx_By_rgt = 0.0d0

!     flx_Bz_lft = 0.0d0
!     flx_Bz_rgt = 0.0d0   







!       if (dusty_nonideal_MHD_no_electron) then !Additional terms in the fluxes for the induction equation

!             idust = i_coupled_species
        
!             i_rho= irhod(idust)

!             i_u  = index_vdn(idust,idim)
!             i_v  = index_vdt(idust,idim)
!             i_w  = ivdz(idust)
!             !print *, idust, i_rho,i_n,i_t,i_z

!             !Dust momentum
!             u_rgt     = qright(i_u)
!             u_lft     = qleft(i_u)
!               !Dust transverse momentum
!             v_rgt     = qright(i_v)
!             v_lft     = qleft(i_v)
!               !Dust second transverse momentum
!             w_rgt     = qright(i_w)
!             w_lft     = qleft(i_w)

!             rho_rgt   = qright(i_rho)
!             rho_lft   = qleft(i_rho)


!             ix=ixx(i)
!             iy=iyy(i)

!             if(slope_type>0) then
!                 il = icell(ix-1,iy)
!                 ir = icell(ix+1,iy)
    
!                 deta_Hall = slope_limit(2.0d0*(abs(eta_eff_Hall_y(i)) - abs(eta_eff_Hall_y(il)))/(dx(i,1)+dx(il,1)),2.0d0*(abs(eta_eff_Hall_y(ir)) - abs(eta_eff_Hall_y(i)))/(dx(ir,1)+dx(i,1)))
!                 deta_Hall_l = slope_limit(2.0d0*(abs(eta_eff_Hall_y(i-1)) - abs(eta_eff_Hall_y(il-1)))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(abs(eta_eff_Hall_y(ir-1)) - abs(eta_eff_Hall_y(i-1)))/(dx(ir-1,1)+dx(i-1,1)))

!                 eta_Hall_y_left = abs(eta_eff_Hall_y(il)) + half*deta_Hall_l*dx(il,1)
!                 eta_Hall_y_right = abs(eta_eff_Hall_y(i))- half*deta_Hall*dx(i,1)

!                 eta_Hall_z_left = - eta_Hall_y_left
!                 eta_Hall_z_right = - eta_Hall_y_right

!                 dJy = slope_limit(2.0d0*(Jy(i) - Jy(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jy(ir) - Jy(i))/(dx(ir,1)+dx(i,1))) !Total current
!                 dJy_l = slope_limit(2.0d0*(Jy(i-1) - Jy(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jy(ir-1) - Jy(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total current

!                 dJz = slope_limit(2.0d0*(Jz(i) - Jz(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jz(ir) - Jz(i))/(dx(ir,1)+dx(i,1))) !Total current
!                 dJz_l = slope_limit(2.0d0*(Jz(i-1) - Jz(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jz(ir-1) - Jz(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total current

!                 Jy_left = Jy(il) + half*dJy_l*dx(il,1)
!                 Jy_right = Jy(i) - half*dJy_l*dx(il,1)

!                 Jz_left = Jz(il) + half*dJz_l*dx(il,1)
!                 Jz_right = Jz(i) - half*dJz_l*dx(il,1)




!             endif



!         !Beware: here signs are reversed with respect to predictor step (fluxes are defined in the left-hand side of the equation) 

!         !Hall term
!         flx_By_lft = flx_By_lft + eta_Hall_y_left*Jy_left
!         flx_By_rgt = flx_By_rgt + eta_Hall_y_right*Jy_right

!         flx_Bz_lft = flx_Bz_lft - eta_Hall_z_left*Jz_left
!         flx_Bz_rgt = flx_Bz_rgt - eta_Hall_z_right*Jz_right

!     endif





! !HLL


!     idust=i_coupled_species !Is the grain species considered in the magnetosonic/Alfven velocity expressions

!     i_rho= irhod(idust)

!     i_u  = index_vdn(idust,idim)
!     i_v  = index_vdt(idust,idim)
!     i_w  = ivdz(idust)
!     !print *, idust, i_rho,i_n,i_t,i_z

!     !Dust momentum
!     u_rgt     = qright(i_u)
!     u_lft     = qleft(i_u)
!       !Dust transverse momentum
!     v_rgt     = qright(i_v)
!     v_lft     = qleft(i_v)
!       !Dust second transverse momentum
!     w_rgt     = qright(i_w)
!     w_lft     = qleft(i_w)

!     rho_rgt   = qright(i_rho)
!     rho_lft   = qleft(i_rho)




!     ca_lft =dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)/dsqrt(4*pi*rho_lft)
!     ca_rgt =dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2)/dsqrt(4*pi*rho_rgt)



!     !Hall effect introduces new waves

!     cw_lft = eta_Hall_y_left*pi/(2*dx(i,1)) + dsqrt((eta_Hall_y_left*pi/(2*dx(i,1)))**2 + ca_lft**2) !Whistler wave
!     cw_rgt = eta_Hall_y_right*pi/(2*dx(i+1,1)) + dsqrt((eta_Hall_y_right*pi/(2*dx(i+1,1)))**2 + ca_rgt**2) !Whistler wave 

!     S_rgt  = max(max(u_lft,u_rgt) +max(cw_lft,cw_rgt),0.0d0) 
!     S_lft  = min(min(u_lft,u_rgt) -max(cw_lft,cw_rgt),0.0d0)



!     flx(iBy)            = flx(iBy) + (S_rgt*flx_By_lft  -S_lft*flx_By_rgt  + S_rgt*S_lft*(By_rgt-By_lft))      / (S_rgt-S_lft)
!     flx(iBz)            = flx(iBz) + (S_rgt*flx_Bz_lft  -S_lft*flx_Bz_rgt  + S_rgt*S_lft*(Bz_rgt-Bz_lft))      / (S_rgt-S_lft)

! end subroutine solver_Hall_hll
! #endif
! #endif
! #endif




#if MHD==1
#if NDUST>0
#if SOLVERB==3

subroutine solver_common_induction_multifluid_hll(qleft,qright,flx,csl,csr,idim,i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!To make it clearer, we write here the HLL solver for the MHD multifluid setup. The solver is common to B and all the dust fludis relying on the effective magnetocompressive speed defined in Verrier+26!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    real(dp) :: ca_lft,ca_rgt,magnetosonic_fast_lft,magnetosonic_fast_rgt,csl,csr,cw_lft,cw_rgt,dJy,dJz,dJy_l,dJz_l,Jy_left,Jy_right,Jz_left,Jz_right,c_fast_lft,c_fast_rgt

    real(dp) :: deta_o_il,deta_h_il,deta_a_il,dzd,dzd_il,zd_left,zd_right,deta_Hall,deta_Hall_l,eta_Hall_y_left,eta_Hall_y_right,eta_Hall_z_left,eta_Hall_z_right

    real(dp) :: nd_zd_over_ni,nd_zd_over_ni_left,nd_zd_over_ni_right,dne,dne_il,ne_left,ne_right,dni,dni_il,ni_left,ni_right,dhall_i,dhall_il,hall_i_left,hall_i_right
    real(dp) :: flx_Bx_lft,flx_Bx_rgt,flx_By_lft,flx_By_rgt,flx_Bz_lft,flx_Bz_rgt,total_dust_current_z_lft,total_dust_current_y_lft,total_dust_current_x_lft,total_dust_current_z_rgt,total_dust_current_y_rgt,total_dust_current_x_rgt,B_norm_lft,B_norm_rgt
    real(dp) :: deta_o,deta_o_l,eta_o_left,eta_o_right,dJdx_tot,dJdx_tot_l,Jdx_tot_left,Jdx_tot_right,dJdy_tot,dJdy_tot_l,Jdy_tot_left,Jdy_tot_right,dJdz_tot,dJdz_tot_l,Jdz_tot_left,Jdz_tot_right
    real(dp) :: deta_H,deta_H_l,eta_H_left,eta_H_right,deta_a,deta_a_l,eta_a_left,eta_a_right,db_unit_x,db_unit_x_l,db_unit_y,db_unit_y_l,db_unit_z,db_unit_z_l,b_unit_x_left,b_unit_x_right,b_unit_y_left,b_unit_y_right,b_unit_z_left,b_unit_z_right
    real(dp) :: c_ms_d_left, c_ms_d_right, dc_ms_d, dc_ms_d_l

    Bx_lft   = qleft(iBx)
    Bx_rgt   = qright(iBx)
    By_lft   = qleft(iBy)
    By_rgt   = qright(iBy)
    Bz_lft   = qleft(iBz)
    Bz_rgt   = qright(iBz)

  

      if (dusty_nonideal_MHD .eqv. .true.) then !Additional terms in the fluxes for the induction equation


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

        ix=ixx(i)
        iy=iyy(i)

            if(slope_type>0) then
                il = icell(ix-1,iy)
                ir = icell(ix+1,iy)

                deta_o = slope_limit(2.0d0*(eta_o(i) - eta_o(il))/(dx(i,1)+dx(il,1)),2.0d0*(eta_o(ir) - eta_o(i))/(dx(ir,1)+dx(i,1)))
                deta_o_l = slope_limit(2.0d0*(eta_o(i-1) - eta_o(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(eta_o(ir-1) - eta_o(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                eta_o_left = eta_o(il) + half*deta_o_l*dx(il,1)
                eta_o_right = eta_o(i) - half*deta_o*dx(i,1)

                deta_H = slope_limit(2.0d0*(eta_H(i) - eta_H(il))/(dx(i,1)+dx(il,1)),2.0d0*(eta_H(ir) - eta_H(i))/(dx(ir,1)+dx(i,1)))
                deta_H_l = slope_limit(2.0d0*(eta_H(i-1) - eta_H(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(eta_H(ir-1) - eta_H(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                eta_H_left = eta_H(il) + half*deta_H_l*dx(il,1)
                eta_H_right = eta_H(i) - half*deta_H*dx(i,1)

                deta_a = slope_limit(2.0d0*(eta_a(i) - eta_a(il))/(dx(i,1)+dx(il,1)),2.0d0*(eta_a(ir) - eta_a(i))/(dx(ir,1)+dx(i,1)))
                deta_a_l = slope_limit(2.0d0*(eta_a(i-1) - eta_a(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(eta_a(ir-1) - eta_a(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                eta_a_left = eta_a(il) + half*deta_a_l*dx(il,1)
                eta_a_right = eta_a(i) - half*deta_a*dx(i,1)



                dJy = slope_limit(2.0d0*(Jy(i) - Jy(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jy(ir) - Jy(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJy_l = slope_limit(2.0d0*(Jy(i-1) - Jy(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jy(ir-1) - Jy(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total current

                dJz = slope_limit(2.0d0*(Jz(i) - Jz(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jz(ir) - Jz(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJz_l = slope_limit(2.0d0*(Jz(i-1) - Jz(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jz(ir-1) - Jz(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total current

                Jy_left = Jy(il) + half*dJy_l*dx(il,1)
                Jy_right = Jy(i) - half*dJy*dx(i,1)

                Jz_left = Jz(il) + half*dJz_l*dx(il,1)
                Jz_right = Jz(i) - half*dJz*dx(i,1)



                dJdx_tot = slope_limit(2.0d0*(Jdx_tot(i) - Jdx_tot(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jdx_tot(ir) - Jdx_tot(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJdx_tot_l = slope_limit(2.0d0*(Jdx_tot(i-1) - Jdx_tot(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jdx_tot(ir-1) - Jdx_tot(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total dust current

                dJdy_tot = slope_limit(2.0d0*(Jdy_tot(i) - Jdy_tot(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jdy_tot(ir) - Jdy_tot(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJdy_tot_l = slope_limit(2.0d0*(Jdy_tot(i-1) - Jdy_tot(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jdy_tot(ir-1) - Jdy_tot(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total dust current

                dJdz_tot = slope_limit(2.0d0*(Jdz_tot(i) - Jdz_tot(il))/(dx(i,1)+dx(il,1)),2.0d0*(Jdz_tot(ir) - Jdz_tot(i))/(dx(ir,1)+dx(i,1))) !Total current
                dJdz_tot_l = slope_limit(2.0d0*(Jdz_tot(i-1) - Jdz_tot(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(Jdz_tot(ir-1) - Jdz_tot(i-1))/(dx(ir-1,1)+dx(i-1,1))) !Total dust current


                Jdx_tot_left = Jdx_tot(il) + half*dJdx_tot_l*dx(il,1)
                Jdx_tot_right = Jdx_tot(i) - half*dJdx_tot*dx(i,1)

                Jdy_tot_left = Jdy_tot(il) + half*dJdy_tot_l*dx(il,1)
                Jdy_tot_right = Jdy_tot(i) - half*dJdy_tot*dx(i,1)

                Jdz_tot_left = Jdz_tot(il) + half*dJdz_tot_l*dx(il,1)
                Jdz_tot_right = Jdz_tot(i) - half*dJdz_tot*dx(i,1)


                db_unit_x = slope_limit(2.0d0*(b_unit_x(i) - b_unit_x(il))/(dx(i,1)+dx(il,1)),2.0d0*(b_unit_x(ir) - b_unit_x(i))/(dx(ir,1)+dx(i,1))) 
                db_unit_y = slope_limit(2.0d0*(b_unit_y(i) - b_unit_y(il))/(dx(i,1)+dx(il,1)),2.0d0*(b_unit_y(ir) - b_unit_y(i))/(dx(ir,1)+dx(i,1))) 
                db_unit_z = slope_limit(2.0d0*(b_unit_z(i) - b_unit_z(il))/(dx(i,1)+dx(il,1)),2.0d0*(b_unit_z(ir) - b_unit_z(i))/(dx(ir,1)+dx(i,1))) 

                db_unit_x_l = slope_limit(2.0d0*(b_unit_x(i-1) - b_unit_x(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(b_unit_x(ir-1) - b_unit_x(i-1))/(dx(ir-1,1)+dx(i-1,1))) 
                db_unit_y_l = slope_limit(2.0d0*(b_unit_y(i-1) - b_unit_y(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(b_unit_y(ir-1) - b_unit_y(i-1))/(dx(ir-1,1)+dx(i-1,1))) 
                db_unit_z_l = slope_limit(2.0d0*(b_unit_z(i-1) - b_unit_z(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(b_unit_z(ir-1) - b_unit_z(i-1))/(dx(ir-1,1)+dx(i-1,1))) 

                b_unit_x_left = b_unit_x(il) + half*db_unit_x_l*dx(il,1)
                b_unit_x_right = b_unit_x(i) - half*db_unit_x*dx(i,1)


                b_unit_y_left = b_unit_y(il) + half*db_unit_y_l*dx(il,1)
                b_unit_y_right = b_unit_y(i) - half*db_unit_y*dx(i,1)

                b_unit_z_left = b_unit_z(il) + half*db_unit_z_l*dx(il,1)
                b_unit_z_right = b_unit_z(i) - half*db_unit_z*dx(i,1)



                dc_ms_d = slope_limit(2.0d0*(c_ms_d(i) - c_ms_d(il))/(dx(i,1)+dx(il,1)),2.0d0*(c_ms_d(ir) - c_ms_d(i))/(dx(ir,1)+dx(i,1)))
                dc_ms_d_l = slope_limit(2.0d0*(c_ms_d(i-1) - c_ms_d(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(c_ms_d(ir-1) - c_ms_d(i-1))/(dx(ir-1,1)+dx(i-1,1)))
                c_ms_d_left = c_ms_d(il) + half*dc_ms_d_l*dx(il,1)
                c_ms_d_right = c_ms_d(i) - half*dc_ms_d*dx(i,1)


            endif


        B_norm_lft = dsqrt(Bx_lft**2+By_lft**2+Bz_lft**2)
        B_norm_rgt = dsqrt(Bx_rgt**2+By_rgt**2+Bz_rgt**2) 


        !Beware: here signs are reversed with respect to predictor step (fluxes are defined in the left-hand side of the equation) 


        flx_By_lft = -Bx_lft*v_lft + By_lft*u_lft !In this model, this term involves the gas velocity
        flx_By_rgt = -Bx_rgt*v_rgt + By_rgt*u_rgt

        flx_Bz_lft = -Bx_lft*w_lft + Bz_lft*u_lft
        flx_Bz_rgt = -Bx_rgt*w_rgt + Bz_rgt*u_rgt


        !Additional term
        !Ohm
        flx_By_lft = flx_By_lft + clight*eta_o_left*Jdz_tot_left
        flx_By_rgt = flx_By_rgt + clight*eta_o_right*Jdz_tot_right

        flx_Bz_lft = flx_Bz_lft - clight*eta_o_left*Jdy_tot_left
        flx_Bz_rgt = flx_Bz_rgt - clight*eta_o_right*Jdy_tot_right

        !AD
        flx_By_lft = flx_By_lft + clight*eta_a_left*( (b_unit_x_left**2+b_unit_y_left**2)*Jdz_tot_left - b_unit_x_left*b_unit_z_left*Jdx_tot_left - b_unit_y_left*b_unit_z_left*Jdy_tot_left) - clight**2/(4*pi) * eta_a_left*b_unit_y_left*b_unit_z_left*(-Jy_left)
        flx_By_rgt = flx_By_rgt + clight*eta_a_right*( (b_unit_x_right**2+b_unit_y_right**2)*Jdz_tot_right - b_unit_x_right*b_unit_z_right*Jdx_tot_right - b_unit_y_right*b_unit_z_right*Jdy_tot_right) - clight**2/(4*pi) * eta_a_right*b_unit_y_right*b_unit_z_right*(-Jy_right)

        flx_Bz_lft = flx_Bz_lft - clight*eta_a_left*((b_unit_x_left**2+b_unit_z_left**2)*Jdy_tot_left - b_unit_y_left*b_unit_z_left*Jdz_tot_left - b_unit_x_left*b_unit_y_left*Jdx_tot_left) - clight**2/(4*pi) * eta_a_left*b_unit_y_left*b_unit_z_left*(Jz_left)
        flx_Bz_rgt = flx_Bz_rgt - clight*eta_a_right*((b_unit_x_right**2+b_unit_z_right**2)*Jdy_tot_right - b_unit_y_right*b_unit_z_right*Jdz_tot_right - b_unit_x_right*b_unit_y_right*Jdx_tot_right) - clight**2/(4*pi) * eta_a_right*b_unit_y_right*b_unit_z_right*(Jz_right)

       

        if (Hall_effect) then
            !Hall term
            flx_By_lft = flx_By_lft + clight*eta_H_left*b_unit_y_left*Jdx_tot_left + clight**2/(4*pi) * eta_H_left * b_unit_x_left * Jy_left - clight*eta_H_left * b_unit_x_left * Jdy_tot_left
            flx_By_rgt = flx_By_rgt + clight*eta_H_right*b_unit_y_right*Jdx_tot_right + clight**2/(4*pi) * eta_H_right * b_unit_x_right * Jy_right - clight*eta_H_right * b_unit_x_right * Jdy_tot_right

            flx_Bz_lft = flx_Bz_lft + clight*eta_H_left*b_unit_z_left*Jdx_tot_left + clight**2/(4*pi) * eta_H_left * b_unit_x_left * Jz_left - clight*eta_H_left*b_unit_x_left*Jdz_tot_left
            flx_Bz_rgt = flx_Bz_rgt + clight*eta_H_right*b_unit_z_right*Jdx_tot_right + clight**2/(4*pi) * eta_H_right * b_unit_x_right * Jz_right - clight*eta_H_right*b_unit_x_right*Jdz_tot_right

        endif


!HLL



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



        S_rgt  = max(max(u_lft,u_rgt) +max(c_ms_d_left,c_ms_d_right),0.0d0) 
        S_lft  = min(min(u_lft,u_rgt) -max(c_ms_d_left,c_ms_d_right),0.0d0)




        if (Hall_effect) then
        !Hall effect introduces new waves

            cw_lft = eta_H_left*pi/(2*dx(i,1)) + dsqrt((eta_H_left*pi/(2*dx(i,1)))**2 + c_ms_d_left**2) !Whistler wave
            cw_rgt = eta_H_right*pi/(2*dx(i+1,1)) + dsqrt((eta_H_right*pi/(2*dx(i+1,1)))**2 + c_ms_d_right**2) !Whistler wave 

            S_rgt  = max(max(u_lft,u_rgt) +max(cw_lft,cw_rgt),0.0d0) 
            S_lft  = min(min(u_lft,u_rgt) -max(cw_lft,cw_rgt),0.0d0)

        endif


        if (S_lft/=S_rgt) then

            flx(iBy)            = (S_rgt*flx_By_lft  -S_lft*flx_By_rgt  + S_rgt*S_lft*(By_rgt-By_lft))      / (S_rgt-S_lft)
            flx(iBz)            = (S_rgt*flx_Bz_lft  -S_lft*flx_Bz_rgt  + S_rgt*S_lft*(Bz_rgt-Bz_lft))      / (S_rgt-S_lft)

        endif


        if (S_lft==S_rgt) then


            lambda_llf_B        = max(abs(u_lft) + abs(c_ms_d_left),abs(u_rgt) + abs(c_ms_d_right))



            if (Hall_effect) lambda_llf_B = max(abs(u_lft) + abs(cw_lft),abs(u_rgt) + abs(cw_rgt)) !In this particular case, only variable B includes whistler speed in the wafe fan. The B solver no longer coincides with the dust solver. No need to add the dust velocity in lambda_llf_B.


            flx(iBy) = half*(flx_By_lft   + flx_By_rgt)    - half*lambda_llf_B*(By_rgt   - By_lft)
            flx(iBz) = half*(flx_Bz_lft   + flx_Bz_rgt)    - half*lambda_llf_B*(Bz_rgt   - Bz_lft)

        endif


    endif


end subroutine solver_common_induction_multifluid_hll
#endif
#endif
#endif



#if NDUST>0
#if MHD==1
#if SOLVERDUST==3

subroutine solver_common_dust_multifluid_hll(qleft,qright,csl,csr,flx,idim,i)

    use parameters
    use commons
    use slope_limiter


    implicit none

    real(dp),dimension(1:nvar),intent(in) :: qright,qleft
    real(dp),dimension(1:nvar),intent(inout) :: flx
    integer  :: idim,idust,i_u,i_v,i_rho,i_w,i,il,ir,ix,iy,icell,ixx,iyy,ipscal
    real(dp) :: csl,csr,P_lft,P_rgt
    real(dp) :: S_lft,S_rgt,lambda_llf_d
    real(dp) :: ca_lft,ca_rgt,cw_rgt,cw_lft,magnetosonic_fast_rgt,magnetosonic_fast_lft



    real(dp) :: rho_lft,rho_rgt,u_lft,u_rgt,v_lft,v_rgt,w_lft,w_rgt
    real(dp) :: mom_u_lft,mom_u_rgt,mom_v_lft,mom_v_rgt,mom_w_lft,mom_w_rgt
    real(dp) :: flx_rho_lft,flx_mom_u_lft,flx_mom_v_lft,flx_mom_w_lft,flx_P_lft
    real(dp) :: flx_rho_rgt,flx_mom_u_rgt,flx_mom_v_rgt,flx_mom_w_rgt,flx_P_rgt,flx_pscal_lft,flx_pscal_rgt
    real(dp) :: Bx_lft,By_lft,Bz_lft,Bx_rgt,By_rgt,Bz_rgt,P_mag_lft,P_mag_rgt,mag_tension_y_lft,mag_tension_y_rgt,mag_tension_z_lft,mag_tension_z_rgt,mag_tension_x_lft,mag_tension_x_rgt,c_fast_lft,c_fast_rgt
    real(dp) :: dc_ms_d,dc_ms_d_l,c_ms_d_left,c_ms_d_right
    
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



    !----------------------------------------------------------------------------------------------------------
    !Within the dusty_nonideal_MHD model, we expect dusty Alfvén waves (magnetocompressive modes) to propagate. 
    !----------------------------------------------------------------------------------------------------------


    ix=ixx(i)
    iy=iyy(i)

    if(slope_type>0) then

        il = icell(ix-1,iy)
        ir = icell(ix+1,iy)
        dc_ms_d = slope_limit(2.0d0*(c_ms_d(i) - c_ms_d(il))/(dx(i,1)+dx(il,1)),2.0d0*(c_ms_d(ir) - c_ms_d(i))/(dx(ir,1)+dx(i,1)))
        dc_ms_d_l = slope_limit(2.0d0*(c_ms_d(i-1) - c_ms_d(il-1))/(dx(i-1,1)+dx(il-1,1)),2.0d0*(c_ms_d(ir-1) - c_ms_d(i-1))/(dx(ir-1,1)+dx(i-1,1)))
        c_ms_d_left = c_ms_d(il) + half*dc_ms_d_l*dx(il,1)
        c_ms_d_right = c_ms_d(i) - half*dc_ms_d*dx(i,1)

    endif


    S_rgt  = max(max(u_lft,u_rgt) +max(c_ms_d_left,c_ms_d_right),0.0d0) 
    S_lft  = min(min(u_lft,u_rgt) -max(c_ms_d_left,c_ms_d_right),0.0d0)



    if (S_rgt/=S_lft) then

        flx(i_rho)            = (S_rgt*flx_rho_lft  -S_lft*flx_rho_rgt  + S_rgt*S_lft*(rho_rgt-rho_lft))      / (S_rgt-S_lft)
        flx(i_u)  = (S_rgt*flx_mom_u_lft-S_lft*flx_mom_u_rgt+ S_rgt*S_lft*(mom_u_rgt-mom_u_lft))  / (S_rgt-S_lft)
        flx(i_v)  = (S_rgt*flx_mom_v_lft-S_lft*flx_mom_v_rgt+ S_rgt*S_lft*(mom_v_rgt-mom_v_lft))  / (S_rgt-S_lft)
        flx(i_w)             = (S_rgt*flx_mom_w_lft-S_lft*flx_mom_w_rgt+ S_rgt*S_lft*(mom_w_rgt-mom_w_lft))  / (S_rgt-S_lft)


#if NDUSTPSCAL>0
        do ipscal=1,ndustpscal
            flx_pscal_lft = qleft(idust_pscal(idust,ipscal)) * rho_lft  * u_lft
            flx_pscal_rgt = qright(idust_pscal(idust,ipscal)) * rho_rgt  * u_rgt
            flx(idust_pscal(idust,ipscal))  = (S_rgt*flx_pscal_lft  - S_lft*flx_pscal_rgt + S_rgt*S_lft*(rho_rgt*qright(idust_pscal(idust,ipscal)) - rho_lft*qleft(idust_pscal(idust,ipscal)))) / (S_rgt-S_lft)
        end do
#endif 

    endif


    if (S_rgt==S_lft) then


        lambda_llf_d        = max(abs(c_ms_d_left),abs(c_ms_d_right))



        flx(i_rho) = half*(flx_rho_lft   + flx_rho_rgt)    - half*lambda_llf_d*(rho_rgt   - rho_lft)
        flx(i_u)  = half*(flx_mom_u_lft + flx_mom_u_rgt)  - half*lambda_llf_d*(mom_u_rgt - mom_u_lft)
        flx(i_v)  = half*(flx_mom_v_lft + flx_mom_v_rgt)  - half*lambda_llf_d*(mom_v_rgt - mom_v_lft)
        flx(i_w)  = half*(flx_mom_w_lft + flx_mom_w_rgt)  - half*lambda_llf_d*(mom_w_rgt - mom_w_lft)

#if NDUSTPSCAL>0
        do ipscal=1,ndustpscal
            flx_pscal_lft = qleft(idust_pscal(idust,ipscal)) * rho_lft  * u_lft
            flx_pscal_rgt = qright(idust_pscal(idust,ipscal)) * rho_rgt  * u_rgt
            flx(idust_pscal(idust,ipscal))  = half*(flx_pscal_lft  + flx_pscal_rgt)   - half*lambda_llf_d*(rho_rgt*qright(idust_pscal(idust,ipscal)) - rho_lft*qleft(idust_pscal(idust,ipscal)))
        end do
#endif 

    endif

end do

end subroutine solver_common_dust_multifluid_hll

#endif 
#endif
#endif

end module hydro_solvers
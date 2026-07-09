! =========================================================================================================================================================
! Here we solve the drag-like terms of the dust / gas momentum equations that have to be treated implicitly (within the full multidusty fluid setup).
! Since the Krapp solver is not applicable, we opt for a brute force linear system resolution (equivalent to matrix inversion although it comes with less
! numerical errors). 
! Note that each term (gyro, Ohm, AD and Hall) are treated separately, i.e. a splitting is done to simplify the system, which comes with a splitting error.
! =========================================================================================================================================================


#if MHD==1
#if NDUST>0
subroutine gyro_drift(i)
  
  use parameters
  use phys_const
  use commons
  use units
  use OMP_LIB
  use lapack_tools

  implicit none
  integer :: i,idust,k,l,m,index_species_column,index_species_row

  !Variables for coordinate system rotation
  real(dp):: alpha,beta
  real(dp), dimension(3) :: b_unit
  real(dp), dimension(3) :: b_unit_prime
  real(dp), dimension(3) :: b_unit_second
  real(dp), dimension(3) :: delta_v
  real(dp), dimension(3) :: delta_v_prime
  real(dp), dimension(3) :: delta_v_second
  real(dp), dimension(3) :: vgas
  real(dp), dimension(3) :: vgas_prime
  real(dp), dimension(3) :: vgas_second
  real(dp), dimension(3) :: vdust
  real(dp), dimension(3) :: vdust_prime
  real(dp), dimension(3) :: vdust_second
  real(dp), dimension(3) :: vgas_intermediate
  real(dp), dimension(3,3) :: rotation_z
  real(dp), dimension(3,3) :: rotation_y_prime
  real(dp), dimension(3,3) :: rotation_z_back
  real(dp), dimension(3,3) :: rotation_y_prime_back
  real(dp), dimension(2*ndust,2*ndust) :: Identity_matrix


  real(dp), dimension(:)  , allocatable :: W_drift
  real(dp), dimension(:,:)  , allocatable :: M_gyro
  real(dp), dimension(:)  , allocatable :: tau_gyr



  !Variables for LU decomposition
  integer, parameter :: n = 2*ndust    !Matrix_gyro size
  integer, parameter :: nrhs = 1       !Number of right hand side vector
  real(dp), dimension(:)  , allocatable :: right_vector !Right hand side vector in the linear system
  integer :: ipiv(n)                   !Pivot vector
  real(dp), dimension(:)  , allocatable :: output_lin_system !Vector solution
  real(dp), dimension(:)  , allocatable :: delta_v_second_x_storage


  if(static) return

  do m=1,2*ndust
    do l=1,2*ndust
      if (m==l) Identity_matrix(m,l) = 1.d0
      if (m/=l) Identity_matrix(m,l) = 0.d0  
    end do
  end do

  allocate(W_drift(1:2*ndust))
  W_drift=0.0d0

  allocate(M_gyro(1:2*ndust,1:2*ndust))
  M_gyro=0.0d0

  allocate(tau_gyr(1:ndust))
  tau_gyr=0.0d0

  allocate(output_lin_system(1:2*ndust))
  output_lin_system=0.0d0

  allocate(right_vector(1:2*ndust))
  right_vector=0.0d0

  allocate(delta_v_second_x_storage(1:ndust))
  delta_v_second_x_storage=0.0d0

    ! =============================================================================================================================================
    ! Here we deal with the gyro-drift (for Vallucci-Goy +27 setup). This gyromotion is a conservative term. Thus, sympletic schemes are required.
    ! We work with the drift velocities (vd - v) and construct a vector W that is (2*ndust) in size.
    ! =============================================================================================================================================


  ! $OMP PARALLEL &
  ! $OMP DEFAULT(SHARED)&
  ! $OMP PRIVATE(i, idust,l,m, alpha, beta, b_unit, b_unit_prime, b_unit_second,rotation_z, rotation_y_prime)
  ! $OMP DO


    b_unit(1) = b_unit_x(i)
    b_unit(2) = b_unit_y(i)
    b_unit(3) = b_unit_z(i)


    ! ========================================================================================================================
    ! Rotate the vectors to align B-vector with x-axis. In doing so, only y ans z components of the linear system will remain.
    ! ========================================================================================================================



    ! print *, 'Bx',u_prim(i,iBx)

    ! print *, 'By',u_prim(i,iBy)

    ! print *, 'Bz',u_prim(i,iBz)


    ! print *, 'b_unit', b_unit
    ! print*, 'b_unit_norm=',b_unit(1)**2+b_unit(2)**2+b_unit(3)**2


    beta = ATAN2(b_unit_y(i), b_unit_x(i)) !Better  to work with ATAN2 than ACOS (less numerical errors)
    alpha = ATAN2(dsqrt(b_unit_x(i)**2 + b_unit_y(i)**2), b_unit_z(i))

    ! print *, 'alpha2', alpha

    ! print *, 'beta2',beta


    !Create rotation matrices
    rotation_z(1,1) = cos(beta); rotation_z(1,2) = sin(beta); rotation_z(1,3) = 0d0
    rotation_z(2,1) = -sin(beta); rotation_z(2,2) = cos(beta); rotation_z(2,3) = 0d0
    rotation_z(3,1) = 0d0; rotation_z(3,2) = 0d0; rotation_z(3,3) = 1d0

    rotation_y_prime(1,1) = sin(alpha); rotation_y_prime(1,2) = 0d0; rotation_y_prime(1,3) = cos(alpha)
    rotation_y_prime(2,1) = 0d0; rotation_y_prime(2,2) = 1d0; rotation_y_prime(2,3) = 0d0
    rotation_y_prime(3,1) = -cos(alpha); rotation_y_prime(3,2) = 0d0; rotation_y_prime(3,3) = sin(alpha)

    b_unit_prime = 0.d0
    b_unit_second = 0.d0

    !Rotate b vector
    call DGEMV('N', 3, 3, 1.d0, rotation_z, 3, b_unit, 1, 0.d0, b_unit_prime, 1) !First: Rotation of an angle beta around z
    call DGEMV('N', 3, 3, 1.d0, rotation_y_prime, 3, b_unit_prime, 1, 0.d0, b_unit_second, 1) !Second: Rotation of an angle pi/2 - alpha around y'





    k = 1 !Incremental index to build W

   do idust=1,ndust

    delta_v(1) =  q(i,ivdx(idust)) - q(i,ivx)
    delta_v(2) =  q(i,ivdy(idust)) - q(i,ivy)
    delta_v(3) =  q(i,ivdz(idust)) - q(i,ivz)

    delta_v_prime = 0.d0
    delta_v_second = 0.d0

    !Rotate physical drift vector
    call DGEMV('N', 3, 3, 1.d0, rotation_z, 3, delta_v, 1, 0.d0, delta_v_prime, 1) !First: Rotation of an angle beta around z
    call DGEMV('N', 3, 3, 1.d0, rotation_y_prime, 3, delta_v_prime, 1, 0.d0, delta_v_second, 1) !Second: Rotation of an angle pi/2 - alpha around y'

    !Construct numerical drift vector W
    W_drift(k) = delta_v_second(2) !No need for x component since after the rotation the cross-product will lead to dW_x/dt = 0
    W_drift(k+1) = delta_v_second(3)

    !Store x component for later
    delta_v_second_x_storage(idust) = delta_v_second(1)

    k=k+2

    !Test
    ! print*, 'cross-product x of w and b = ', delta_v_second(2)*b_unit_second(3) - delta_v_second(3)*b_unit_second(2)
    ! print*, 'cross-product y of w and b = ', delta_v_second(3)*b_unit_second(1) - delta_v_second(1)*b_unit_second(3)
    ! print*, 'cross-product z of w and b = ', delta_v_second(1)*b_unit_second(2) - delta_v_second(2)*b_unit_second(1)


   end do

   ! print*, 'W_drift = ', W_drift


    ! ==========================================================================================
    ! Solve linear system with Crank-Nicholson scheme (simplectic implicit midpoint integrator).
    ! ==========================================================================================



    !Define gyration time
    do idust=1,ndust
      ! print*, 'zd',ABS(zd(i,idust))
      ! print*, 'abs(zd)',ABS(zd(i,idust))
      ! print*, 'mdust',mdust(i,idust)

      tau_gyr(idust) = clight * mdust(i,idust) / (ABS(zd(i,idust)) * e_el_stat * dsqrt(q(i,iBx)**2 + q(i,iBy)**2 + q(i,iBz)**2))
    end do

    !Build M_gyro

    index_species_row = 1

    do m=1,2*ndust
      if (MOD(m, 2) == 1 .and. m/=1)  index_species_row = index_species_row + 1

      index_species_column = 1

      do l=1,2*ndust
        if (MOD(l, 2) == 1 .and. l/=1)  index_species_column = index_species_column + 1

        if (MOD(m, 2) == 1) then !If row number is odd
          if (MOD(l, 2) == 1) M_gyro(m,l) = 0.d0
          if (MOD(l, 2) == 0 .and. index_species_row /= index_species_column) M_gyro(m,l) = q(i,irhod(index_species_column))/q(i,irho) * 1.d0/tau_gyr(index_species_column) * b_unit_second(1)
          if (MOD(l, 2) == 0 .and. index_species_row == index_species_column) M_gyro(m,l) = (1.d0 + q(i,irhod(index_species_column))/q(i,irho)) * 1.d0/tau_gyr(index_species_column) * b_unit_second(1)
        endif

        if (MOD(m, 2) == 0) then !If row number is even
          if (MOD(l, 2) == 0) M_gyro(m,l) = 0.d0
          if (MOD(l, 2) == 1 .and. index_species_row /= index_species_column) M_gyro(m,l) = -q(i,irhod(index_species_column))/q(i,irho) * 1.d0/tau_gyr(index_species_column) * b_unit_second(1)
          if (MOD(l, 2) == 1 .and. index_species_row == index_species_column) M_gyro(m,l) = -(1.d0 + q(i,irhod(index_species_column))/q(i,irho)) * 1.d0/tau_gyr(index_species_column) * b_unit_second(1)
        endif


      end do


    end do

    ! ====================================
    ! This is the Crank-Nicholson scheme 
    ! ====================================

    !Build left-hand side --> (Id - M_gyro dt / 2) and right-hand side --> (Id + M_gyro dt / 2)*W_drift
    !call DGEMV('N', n, n, 1.d0, Identity_matrix + M_gyro * dt / 2, n, W_drift, 1, 0.d0, right_vector, 1)

    !Solve linear system 
    !call LU_factorization_resolution(Identity_matrix - M_gyro * dt / 2, n, n, right_vector,n,nrhs,output_lin_system) !output_lin_system is now the updated W_drift (t = n+1) 
    !==============================================================================================================
    

    ! ====================================
    ! This is the Euler scheme 
    ! ====================================

    call LU_factorization_resolution(Identity_matrix - M_gyro * dt, n, n, W_drift,n,nrhs,output_lin_system) 
    !==============================================================================================================




    ! ==============================================================================================================
    ! Retrieve individual velocities via momentum conservation and rotate back to the original system of coordinates 
    ! ==============================================================================================================

    ! Back rotation matrices

    rotation_z_back(1,1) = cos(beta); rotation_z_back(1,2) = -sin(beta); rotation_z_back(1,3) = 0d0 !R(-beta)^{-1} = R(beta)
    rotation_z_back(2,1) = sin(beta); rotation_z_back(2,2) = cos(beta); rotation_z_back(2,3) = 0d0
    rotation_z_back(3,1) = 0d0; rotation_z_back(3,2) = 0d0; rotation_z_back(3,3) = 1d0

    rotation_y_prime_back(1,1) = sin(alpha); rotation_y_prime_back(1,2) = 0d0; rotation_y_prime_back(1,3) = -cos(alpha) !R(pi/2 - alpha)^{-1} = R(alpha-pi/2)
    rotation_y_prime_back(2,1) = 0d0; rotation_y_prime_back(2,2) = 1d0; rotation_y_prime_back(2,3) = 0d0
    rotation_y_prime_back(3,1) = cos(alpha); rotation_y_prime_back(3,2) = 0d0; rotation_y_prime_back(3,3) = sin(alpha)

    k = 1

    do idust=1,ndust

      !Update delta_v
      delta_v_second(1) =  delta_v_second_x_storage(idust)
      delta_v_second(2) =  output_lin_system(k)
      delta_v_second(3) =  output_lin_system(k+1)


      !Lift the degenaracy between gas and dust velocities by invoking momentum conservation
      !First: bring former (time n) gas and dust velocities to rotated base

      vdust(1) = q(i,ivdx(idust))
      vdust(2) = q(i,ivdy(idust))
      vdust(3) = q(i,ivdz(idust))

      !Rotation
      call DGEMV('N', 3, 3, 1.d0, rotation_z, 3, vdust, 1, 0.d0, vdust_prime, 1) !First: Rotation of an angle beta around z
      call DGEMV('N', 3, 3, 1.d0, rotation_y_prime, 3, vdust_prime, 1, 0.d0, vdust_second, 1) !Second: Rotation of an angle pi/2 - alpha around y'

      !First occurence only for the gas
      if (k==1) then
        vgas(1) = q(i,ivx)
        vgas(2) = q(i,ivy)
        vgas(3) = q(i,ivz)

        !Rotation
        call DGEMV('N', 3, 3, 1.d0, rotation_z, 3, vgas, 1, 0.d0, vgas_prime, 1) !First: Rotation of an angle beta around z
        call DGEMV('N', 3, 3, 1.d0, rotation_y_prime, 3, vgas_prime, 1, 0.d0, vgas_second, 1) !Second: Rotation of an angle pi/2 - alpha around y'
      endif

      !Then take the previous intermediate gas velocity to keep updating with the next dust species 
      if (k/=1) then
        vgas_second(1) = vgas_intermediate(1)
        vgas_second(2) = vgas_intermediate(2)
        vgas_second(3) = vgas_intermediate(3)
      endif





      !Intermediate gas velocity (Note that vdust_second and vgas_second are the velocities a time n (intermediate for the gas). delta_v_second is the drift at n+1)
      vgas_intermediate(:) = 1 / (1 + q(i,irhod(idust))/q(i,irho)) * (vgas_second(:) + q(i,irhod(idust))/q(i,irho) * (vdust_second(:) - delta_v_second(:))) !Individual mom. conservation by component (and for each dust fluid individually).
      !Updated dust velocity
      vdust_second(:) = delta_v_second(:) + vgas_intermediate(:)


      !Bring dust velocity back to original base
      call DGEMV('N', 3, 3, 1.d0, rotation_y_prime_back, 3, vdust_second, 1, 0.d0, vdust_prime, 1)
      call DGEMV('N', 3, 3, 1.d0, rotation_z_back, 3, vdust_prime, 1, 0.d0, vdust, 1)

      ! print*,'vdust(1)-q(i,ivdx(idust))',vdust(1)-q(i,ivdx(idust))
      ! print*,'vdust(2)-q(i,ivdy(idust))',vdust(2)-q(i,ivdy(idust))
      ! print*,'vdust(3)-q(i,ivdz(idust))',vdust(3)-q(i,ivdz(idust))

      q(i,ivdx(idust)) = vdust(1) !Update q so that the other drag terms work with the updated value
      q(i,ivdy(idust)) = vdust(2)
      q(i,ivdz(idust)) = vdust(3)

      u_prim(i,ivdx(idust)) = u_prim(i,irhod(idust)) * vdust(1) !Update u so that the CFL is computed correctly, with the updated value (needed if this drag term is called last).
      u_prim(i,ivdy(idust)) = u_prim(i,irhod(idust)) * vdust(2)
      u_prim(i,ivdz(idust)) = u_prim(i,irhod(idust)) * vdust(3)




      k = k+2

    enddo

    !Now that all dust species have been updated, we compute the updated gas velocity with the last vgas_intermediate computed
    call DGEMV('N', 3, 3, 1.d0, rotation_y_prime_back, 3, vgas_intermediate, 1, 0.d0, vgas_prime, 1)
    call DGEMV('N', 3, 3, 1.d0, rotation_z_back, 3, vgas_prime, 1, 0.d0, vgas, 1) 

    q(i,ivx) = vgas(1)
    q(i,ivy) = vgas(2)
    q(i,ivz) = vgas(3)  

    u_prim(i,ivx) = u_prim(i,irho) * vgas(1)
    u_prim(i,ivy) = u_prim(i,irho) * vgas(2)
    u_prim(i,ivz) = u_prim(i,irho) * vgas(3)      
      
    ! print*,'vgas(1) - q(i,ivx)',vgas(1) - q(i,ivx)
    ! print*,'vgas(2) - q(i,ivy)',vgas(2) - q(i,ivy)
    ! print*,'vgas(3) - q(i,ivz)',vgas(3) - q(i,ivz)











  deallocate(W_drift)
  deallocate(M_gyro)
  deallocate(tau_gyr)
  deallocate(output_lin_system)
  deallocate(right_vector)
  deallocate(delta_v_second_x_storage)


end subroutine gyro_drift




subroutine Ohmic_drag(i)
  
  use parameters
  use phys_const
  use commons
  use units
  use OMP_LIB
  use lapack_tools

  implicit none
  integer :: i,idust,l,m
  real(dp) :: zd_tot

  real(dp), dimension(ndust+1,ndust+1) :: Identity_matrix


  real(dp), dimension(:)  , allocatable :: Vx
  real(dp), dimension(:)  , allocatable :: Vy
  real(dp), dimension(:)  , allocatable :: Vz
  real(dp), dimension(:,:)  , allocatable :: M_O



  !Variables for LU decomposition
  integer, parameter :: n = ndust+1   !Matrix_Ohm size
  integer, parameter :: nrhs = 1       !Number of right hand side vector
  integer :: ipiv(n)                   !Pivot vector
  real(dp), dimension(:)  , allocatable :: output_lin_system_x !Vector solution
  real(dp), dimension(:)  , allocatable :: output_lin_system_y !Vector solution
  real(dp), dimension(:)  , allocatable :: output_lin_system_z !Vector solution



  if(static) return


  do m=1,ndust+1
    do l=1,ndust+1
      if (m==l) Identity_matrix(m,l) = 1.d0
      if (m/=l) Identity_matrix(m,l) = 0.d0  
    end do
  end do

  allocate(Vx(1:ndust+1))
  Vx=0.0d0

  allocate(Vy(1:ndust+1))
  Vy=0.0d0

  allocate(Vz(1:ndust+1))
  Vz=0.0d0

  allocate(M_O(1:ndust+1,1:ndust+1))
  M_O=0.0d0

  allocate(output_lin_system_x(1:ndust+1))
  output_lin_system_x=0.0d0

  allocate(output_lin_system_y(1:ndust+1))
  output_lin_system_y=0.0d0

  allocate(output_lin_system_z(1:ndust+1))
  output_lin_system_z=0.0d0


    ! ==========================================================================================================================================================
    ! Here we deal with the Ohmic drag term (for Vallucci-Goy +27 setup). This one is a regular dissipative drag-like term. Thus, usual schemes are well suited.
    ! ==========================================================================================================================================================

  !!$OMP PARALLEL &
  !!$OMP DEFAULT(SHARED)&
  !!$OMP PRIVATE(pnx,pny,pnz,rhon,alphak,B_norm,i,idust)
  !!$OMP DO


    ! ==============================================================================
    ! We first build the numerical velocity vector V, one for each direction.
    ! ==============================================================================


    Vx(1) = q(i,ivx)
    Vy(1) = q(i,ivy)
    Vz(1) = q(i,ivz)


    do idust=1,ndust

      Vx(idust+1) = q(i,ivdx(idust))
      Vy(idust+1) = q(i,ivdy(idust))
      Vz(idust+1) = q(i,ivdz(idust))

   end do


    ! ==========================================================================================
    ! Solve linear system with 1st order implicit Euler discretization.
    ! ==========================================================================================


    !Define total dust charge density
    zd_tot = SUM(q(i,irhod(:))/mdust(i,:) * zd(i,:) * e_el_stat)


    !Build M_O
    do m=1,1+ndust
      do l=1,1+ndust

        if (m == 1) then
          if (l == 1) M_O(m,l) = - zd_tot**2 / q(i,irho)
          if (l /= 1) M_O(m,l) = zd_tot * (q(i,irhod(l-1))/mdust(i,l-1) * zd(i,l-1) * e_el_stat) / q(i,irho)
        endif

        if (m /= 1) then
          if (l == 1) M_O(m,l) = zd_tot * (q(i,irhod(m-1))/mdust(i,m-1) * zd(i,m-1) * e_el_stat) / q(i,irhod(m-1))
          if (l /= 1) M_O(m,l) = - (q(i,irhod(m-1))/mdust(i,m-1) * zd(i,m-1) * e_el_stat) * (q(i,irhod(l-1))/mdust(i,l-1) * zd(i,l-1) * e_el_stat) / q(i,irhod(m-1))
        endif

      end do
    end do

    M_O(:,:) = eta_o(i) * M_O(:,:) 
    !Left-hand side --> (Id - M_gyro dt) and right-hand side --> V

    ! print*, 'zd',zd(i,:)
    ! print*, 'zd_tot',zd_tot
    !print*, 'eta_o(i)',eta_o(i)
    ! print*,'q,(i,irho)',q(i,irho)
    ! print*,'q,(i,irhod(idust))',q(i,irhod(:))
    ! print*, 'mdust',mdust(i,:)
    ! print*, 'M_O(1,:)=',M_O(1,:)
    ! print*, 'M_O(2,:)=',M_O(2,:)
    ! print*, 'M_O(3,:)=',M_O(3,:)


    !Solve linear system 
    call LU_factorization_resolution(Identity_matrix - M_O * dt, n, n, Vx,n,nrhs,output_lin_system_x) !output_lin_system is now the updated Vx (t = n+1) 
    call LU_factorization_resolution(Identity_matrix - M_O * dt, n, n, Vy,n,nrhs,output_lin_system_y) !output_lin_system is now the updated Vy (t = n+1) 
    call LU_factorization_resolution(Identity_matrix - M_O * dt, n, n, Vz,n,nrhs,output_lin_system_z) !output_lin_system is now the updated Vz (t = n+1) 


    ! ============================
    ! Update conservative variable
    ! ============================

    ! print*, 'output_lin_system_x - q(i,ivdx(1))=',output_lin_system_x(2) - q(i,ivdx(1))
    ! print*, 'output_lin_system_x - q(i,ivdy(1))=',output_lin_system_y(2) - q(i,ivdy(1))
    ! print*, 'output_lin_system_x - q(i,ivdz(1))=',output_lin_system_z(2) - q(i,ivdz(1))


    do idust=1,ndust

      q(i,ivdx(idust)) = output_lin_system_x(idust+1)
      q(i,ivdy(idust)) = output_lin_system_y(idust+1)
      q(i,ivdz(idust)) = output_lin_system_z(idust+1)

      u_prim(i,ivdx(idust)) = u_prim(i,irhod(idust)) * output_lin_system_x(idust+1)
      u_prim(i,ivdy(idust)) = u_prim(i,irhod(idust)) * output_lin_system_y(idust+1)
      u_prim(i,ivdz(idust)) = u_prim(i,irhod(idust)) * output_lin_system_z(idust+1)

    enddo

    !Now the gas

    q(i,ivx) = output_lin_system_x(1)
    q(i,ivy) = output_lin_system_y(1)
    q(i,ivz) = output_lin_system_z(1) 

    u_prim(i,ivx) = u_prim(i,irho) * output_lin_system_x(1)
    u_prim(i,ivy) = u_prim(i,irho) * output_lin_system_y(1)
    u_prim(i,ivz) = u_prim(i,irho) * output_lin_system_z(1)    

    ! print*, 'output_lin_system_x - q(i,ivx)=',output_lin_system_x(1) - q(i,ivx)
    ! print*, 'output_lin_system_x - q(i,ivy)=',output_lin_system_y(1) - q(i,ivy)
    ! print*, 'output_lin_system_x - q(i,ivz)=',output_lin_system_z(1) - q(i,ivz)  






  deallocate(M_O)
  deallocate(output_lin_system_x)
  deallocate(output_lin_system_y)
  deallocate(output_lin_system_z)
  deallocate(Vx)
  deallocate(Vy)
  deallocate(Vz)



end subroutine Ohmic_drag

subroutine Hall_gyromotion(i)
  
  use parameters
  use phys_const
  use commons
  use units
  use OMP_LIB
  use lapack_tools

  implicit none
  integer :: i,idust,l,m,k,count_species_column,count_species_row,index_species_column,index_species_row
  real(dp) :: zd_tot

  real(dp), dimension(:)  , allocatable :: V
  real(dp), dimension(:,:)  , allocatable :: M_H

  !Variables for LU decomposition
  integer, parameter :: n = 3*(ndust+1)   !Matrix_Hall size
  integer, parameter :: nrhs = 1       !Number of right hand side vector
  integer :: ipiv(n)                   !Pivot vector
  real(dp), dimension(:), allocatable :: right_vector_Hall
  real(dp), dimension(:)  , allocatable :: output_lin_system_Hall !Vector solution



  real(dp), dimension(n,n) :: Identity_matrix




  if(static) return


  do m=1,n
    do l=1,n
      if (m==l) Identity_matrix(m,l) = 1.d0
      if (m/=l) Identity_matrix(m,l) = 0.d0  
    end do
  end do

  allocate(V(1:n))
  V=0.0d0

  allocate(M_H(1:n,1:n))
  M_H=0.0d0

  allocate(output_lin_system_Hall(1:n))
  output_lin_system_Hall=0.0d0

  allocate(right_vector_Hall(1:n))
  right_vector_Hall=0.0d0


    ! ======================================================================================================================================================================================================================
    ! Here we deal with the Hall gyromotion term (for Vallucci-Goy +27 setup). This one is conservative. Similarly to the gyrodrag, we use a symplectic Crank-Nicholson scheme (should work but does not) --> Euler instead
    ! ======================================================================================================================================================================================================================

  !!$OMP PARALLEL &
  !!$OMP DEFAULT(SHARED)&
  !!$OMP PRIVATE(pnx,pny,pnz,rhon,alphak,B_norm,i,idust)
  !!$OMP DO


    ! ================================================
    ! We first build the numerical velocity vector V
    ! ================================================


    V(1) = q(i,ivx)
    V(2) = q(i,ivy)
    V(3) = q(i,ivz)

    k = 4
    do idust=1,ndust

      V(k) = q(i,ivdx(idust))
      V(k+1) = q(i,ivdy(idust))
      V(k+2) = q(i,ivdz(idust))

      k = k + 3

   end do


    ! ====================================================================================================
    ! Solve linear system with Crank-Nicholson scheme (simplectic implicit midpoint integrator), or Euler.
    ! ====================================================================================================


    !Define total dust charge density
    zd_tot = SUM(q(i,irhod(:))/mdust(i,:) * zd(i,:) * e_el_stat)


    ! print*, 'zd',zd(i,:)
    ! print*, 'zd_tot',zd_tot
    ! print*, 'eta_H(i)',eta_H(i)
    ! print*, 'b_unit_x(i)',b_unit_x(i)
    ! print*, 'b_unit_y(i)',b_unit_y(i)
    ! print*, 'b_unit_z(i)',b_unit_z(i)

    !  print*, 'eta_H(i)',eta_H(i)

    ! print*,'q,(i,irho)',q(i,irho)
    ! print*,'q,(i,irhod(idust))',q(i,irhod(:))
    ! print*, 'mdust',mdust(i,:)



    !Build M_H
    count_species_row = 1
    index_species_row = 0



    do m=1,n

      if (count_species_row == 4) then
        index_species_row = index_species_row + 1 !We moved the next dust species (every 3 iterations)
        count_species_row = 1 !Reset the counter
      endif

      ! print*, 'index_species_row',index_species_row

      count_species_column = 1
      index_species_column = 0

      do l=1,n

        if (count_species_column == 4) then 
          index_species_column = index_species_column + 1 !We moved the next dust species (every 3 iterations)
          count_species_column = 1 !Reset the counter
        endif
        ! print*, 'index_species_column',index_species_column
        ! print*, 'count_species_column',count_species_column


        if (index_species_row == 0 .and. index_species_column == 0) then !Top left-hand corner
          if (count_species_row == count_species_column) M_H(m,l) = 0.0d0
          if (count_species_row == 1) then
            if (count_species_column == 2) M_H(m,l) = - zd_tot**2 * b_unit_z(i)/ q(i,irho) 
            if (count_species_column == 3) M_H(m,l) = zd_tot**2 * b_unit_y(i)/ q(i,irho)
          endif 
          if (count_species_row == 2) then
            if (count_species_column == 1) M_H(m,l) = zd_tot**2 * b_unit_z(i)/ q(i,irho) 
            if (count_species_column == 3) M_H(m,l) = - zd_tot**2 * b_unit_x(i)/ q(i,irho)
          endif
          if (count_species_row == 3) then
            if (count_species_column == 1) M_H(m,l) = - zd_tot**2 * b_unit_y(i)/ q(i,irho) 
            if (count_species_column == 2) M_H(m,l) = zd_tot**2 * b_unit_x(i)/ q(i,irho)
          endif
        endif

        if (index_species_row == 0 .and. index_species_column /= 0) then
          if (count_species_row == count_species_column) M_H(m,l) = 0.0d0
          if (count_species_row == 1) then
            if (count_species_column == 2) M_H(m,l) = zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_z(i)/ q(i,irho) 
            if (count_species_column == 3) M_H(m,l) = - zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_y(i)/ q(i,irho) 
          endif 
          if (count_species_row == 2) then
            if (count_species_column == 1) M_H(m,l) = - zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_z(i)/ q(i,irho)  
            if (count_species_column == 3) M_H(m,l) = zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_x(i)/ q(i,irho) 
          endif
          if (count_species_row == 3) then
            if (count_species_column == 1) M_H(m,l) = zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_y(i)/ q(i,irho)  
            if (count_species_column == 2) M_H(m,l) = - zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_x(i)/ q(i,irho) 
          endif
        endif

        if (index_species_row /= 0 .and. index_species_column == 0) then
          if (count_species_row == count_species_column) M_H(m,l) = 0.0d0
          if (count_species_row == 1) then
            if (count_species_column == 2) M_H(m,l) = zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_z(i)/ q(i,irhod(index_species_row)) 
            if (count_species_column == 3) M_H(m,l) = - zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_y(i)/ q(i,irhod(index_species_row)) 
          endif 
          if (count_species_row == 2) then
            if (count_species_column == 1) M_H(m,l) = - zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_z(i)/ q(i,irhod(index_species_row))  
            if (count_species_column == 3) M_H(m,l) = zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_x(i)/ q(i,irhod(index_species_row)) 
          endif
          if (count_species_row == 3) then
            if (count_species_column == 1) M_H(m,l) = zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_y(i)/ q(i,irhod(index_species_row))  
            if (count_species_column == 2) M_H(m,l) = - zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_x(i)/ q(i,irhod(index_species_row)) 
          endif
        endif

        if (index_species_row /= 0 .and. index_species_column /= 0) then
          if (count_species_row == count_species_column) M_H(m,l) = 0.0d0
          if (count_species_row == 1) then
            if (count_species_column == 2) M_H(m,l) = - (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_z(i)/ q(i,irhod(index_species_row)) 
            if (count_species_column == 3) M_H(m,l) = (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_y(i)/ q(i,irhod(index_species_row)) 
          endif 
          if (count_species_row == 2) then
            if (count_species_column == 1) M_H(m,l) = (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_z(i)/ q(i,irhod(index_species_row))  
            if (count_species_column == 3) M_H(m,l) = - (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_x(i)/ q(i,irhod(index_species_row)) 
          endif
          if (count_species_row == 3) then
            if (count_species_column == 1) M_H(m,l) = - (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_y(i)/ q(i,irhod(index_species_row))  
            if (count_species_column == 2) M_H(m,l) = (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_x(i)/ q(i,irhod(index_species_row)) 
          endif
        endif

        count_species_column = count_species_column + 1

      end do

      count_species_row = count_species_row + 1

    end do

    M_H(:,:) = eta_H(i) * M_H(:,:) 


    !Left-hand side --> (Id - M_gyro dt) and right-hand side --> V


    ! print*, 'M_H(1,:)=',M_H(1,:)
    ! print*, 'M_H(2,:)=',M_H(2,:)
    ! print*, 'M_H(3,:)=',M_H(3,:)
    ! print*, 'M_H(4,:)=',M_H(4,:)
    ! print*, 'M_H(5,:)=',M_H(5,:)
    ! print*, 'M_H(6,:)=',M_H(6,:)
    ! print*, 'M_H(7,:)=',M_H(7,:)
    ! print*, 'M_H(8,:)=',M_H(8,:)
    ! print*, 'M_H(9,:)=',M_H(9,:)



    ! print*,'q(i,irho)',q(i,irho)
    ! print*,'q(i,irhod(:))',q(i,irhod(:))
    ! print*,'q,(i,ivx)',q(i,ivx)
    ! print*,'q,(i,ivy)',q(i,ivy)
    ! print*,'q,(i,ivz)',q(i,ivz)
    ! print*,'q,(i,ivdx(1))',q(i,ivdx(1))
    ! print*,'q,(i,ivdy(1))',q(i,ivdy(1))
    ! print*,'q,(i,ivdz(1))',q(i,ivdz(1))
    ! print*,'q,(i,ivdx(2))',q(i,ivdx(2))
    ! print*,'q,(i,ivdy(2))',q(i,ivdy(2))
    ! print*,'q,(i,ivdz(2))',q(i,ivdz(2))

    ! print*, 'V(:)=',V(:)
    ! print*,'dt=',dt


    ! ====================================
    ! This is the Crank-Nicholson scheme 
    ! ====================================

    ! call DGEMV('N', n, n, 1.d0, Identity_matrix + M_H * dt / 2, n, V, 1, 0.d0, right_vector_Hall, 1)

    ! print*, 'right_vector_Hall(:)=',right_vector_Hall(:)


    ! !Solve linear system 
    ! call LU_factorization_resolution(Identity_matrix - M_H * dt / 2, n, n, right_vector_Hall,n,nrhs,output_lin_system_Hall) !output_lin_system is now the updated V (t = n+1) 


    ! ====================================
    ! This is the implicit Euler scheme 
    ! ====================================

    call LU_factorization_resolution(Identity_matrix - M_H * dt, n, n, V,n,nrhs,output_lin_system_Hall) 


    !print*, 'output_lin_system_Hall(:)=',output_lin_system_Hall(:)

    ! ============================
    ! Update conservative variable
    ! ============================
    k = 4
    do idust=1,ndust


      q(i,ivdx(idust)) = output_lin_system_Hall(k)
      q(i,ivdy(idust)) = output_lin_system_Hall(k+1)
      q(i,ivdz(idust)) = output_lin_system_Hall(k+2)

      u_prim(i,ivdx(idust)) = u_prim(i,irhod(idust)) * output_lin_system_Hall(k)
      u_prim(i,ivdy(idust)) = u_prim(i,irhod(idust)) * output_lin_system_Hall(k+1)
      u_prim(i,ivdz(idust)) = u_prim(i,irhod(idust)) * output_lin_system_Hall(k+2)

      ! print*,'output_lin_system_Hall(k)-q(i,ivdx(idust))',output_lin_system_Hall(k)-q(i,ivdx(idust)) 
      ! print*,'output_lin_system_Hall(k+1)-q(i,ivdy(idust))',output_lin_system_Hall(k+1)-q(i,ivdy(idust)) 
      ! print*,'output_lin_system_Hall(k+2)-q(i,ivdz(idust))',output_lin_system_Hall(k+2)-q(i,ivdz(idust)) 


      k = k + 3
      !print*,'k',k
    enddo

    !Now the gas

    q(i,ivx) = output_lin_system_Hall(1)
    q(i,ivy) = output_lin_system_Hall(2)
    q(i,ivz) = output_lin_system_Hall(3)

    u_prim(i,ivx) = u_prim(i,irho) * output_lin_system_Hall(1)
    u_prim(i,ivy) = u_prim(i,irho) * output_lin_system_Hall(2)
    u_prim(i,ivz) = u_prim(i,irho) * output_lin_system_Hall(3)   

    ! print*,'output_lin_system_Hall(1)-q(i,ivx)',output_lin_system_Hall(1)-q(i,ivx) 
    ! print*,'output_lin_system_Hall(2)-q(i,ivy)',output_lin_system_Hall(2)-q(i,ivy) 
    ! print*,'output_lin_system_Hall(3)-q(i,ivz)',output_lin_system_Hall(3)-q(i,ivz) 







  deallocate(M_H)
  deallocate(output_lin_system_Hall)
  deallocate(right_vector_Hall)
  deallocate(V)




end subroutine Hall_gyromotion



subroutine AD_drag(i)
  
  use parameters
  use phys_const
  use commons
  use units
  use OMP_LIB
  use lapack_tools

  implicit none
  integer :: i,idust,l,m,k,count_species_column,count_species_row,index_species_column,index_species_row
  real(dp) :: zd_tot

  real(dp), dimension(:)  , allocatable :: V
  real(dp), dimension(:,:)  , allocatable :: M_AD

  !Variables for LU decomposition
  integer, parameter :: n = 3*(ndust+1)   !Matrix_AD size
  integer, parameter :: nrhs = 1       !Number of right hand side vector
  integer :: ipiv(n)                   !Pivot vector
  real(dp), dimension(:)  , allocatable :: output_lin_system_AD !Vector solution


  real(dp), dimension(n,n) :: Identity_matrix




  if(static) return


  do m=1,n
    do l=1,n
      if (m==l) Identity_matrix(m,l) = 1.d0
      if (m/=l) Identity_matrix(m,l) = 0.d0  
    end do
  end do

  allocate(V(1:n))
  V=0.0d0

  allocate(M_AD(1:n,1:n))
  M_AD=0.0d0

  allocate(output_lin_system_AD(1:n))
  output_lin_system_AD=0.0d0




    ! ========================================================================================================================================================================
    ! Here we deal with the AD drag term (for Vallucci-Goy +27 setup). This one is dissipative. Similarly to the Ohmic drag, we use a first-order Euler implicit scheme
    ! ========================================================================================================================================================================

  !!$OMP PARALLEL &
  !!$OMP DEFAULT(SHARED)&
  !!$OMP PRIVATE(pnx,pny,pnz,rhon,alphak,B_norm,i,idust)
  !!$OMP DO


    ! ================================================
    ! We first build the numerical velocity vector V
    ! ================================================


    V(1) = q(i,ivx)
    V(2) = q(i,ivy)
    V(3) = q(i,ivz)

    k = 4
    do idust=1,ndust

      V(k) = q(i,ivdx(idust))
      V(k+1) = q(i,ivdy(idust))
      V(k+2) = q(i,ivdz(idust))

      k = k + 3

   end do


    ! ================================================
    ! Solve linear system with Euler implicit scheme.
    ! ================================================


    !Define total dust charge density
    zd_tot = SUM(q(i,irhod(:))/mdust(i,:) * zd(i,:) * e_el_stat)



    ! print*, 'zd',zd(i,:)
    ! print*, 'zd_tot',zd_tot
    ! print*, 'b_unit_x(i)',b_unit_x(i)
    ! print*, 'b_unit_y(i)',b_unit_y(i)
    ! print*, 'b_unit_z(i)',b_unit_z(i)

    !  print*, 'eta_a(i)',eta_a(i)

    ! print*,'q,(i,irho)',q(i,irho)
    ! print*,'q,(i,irhod(idust))',q(i,irhod(:))
    ! print*, 'mdust',mdust(i,:)



    !Build M_AD
    count_species_row = 1
    index_species_row = 0



    do m=1,n

      if (count_species_row == 4) then
        index_species_row = index_species_row + 1 !We moved the next dust species (every 3 iterations)
        count_species_row = 1 !Reset the counter
      endif

      ! print*, 'index_species_row',index_species_row

      count_species_column = 1
      index_species_column = 0

      do l=1,n

        if (count_species_column == 4) then 
          index_species_column = index_species_column + 1 !We moved the next dust species (every 3 iterations)
          count_species_column = 1 !Reset the counter
        endif
        ! print*, 'index_species_column',index_species_column
        ! print*, 'count_species_column',count_species_column


        if (index_species_row == 0 .and. index_species_column == 0) then !Top left-hand corner
          if (count_species_row == 1 .and. count_species_column == 1) M_AD(m,l) = - zd_tot**2 * (b_unit_y(i)**2 + b_unit_z(i)**2) / q(i,irho)
          if (count_species_row == 2 .and. count_species_column == 2) M_AD(m,l) = - zd_tot**2 * (b_unit_z(i)**2 + b_unit_x(i)**2) / q(i,irho)
          if (count_species_row == 3 .and. count_species_column == 3) M_AD(m,l) = - zd_tot**2 * (b_unit_x(i)**2 + b_unit_y(i)**2) / q(i,irho)

          if (count_species_row == 1) then
            if (count_species_column == 2) M_AD(m,l) = zd_tot**2 * b_unit_x(i) * b_unit_y(i) / q(i,irho) 
            if (count_species_column == 3) M_AD(m,l) = zd_tot**2 * b_unit_x(i) * b_unit_z(i) / q(i,irho)
          endif 
          if (count_species_row == 2) then
            if (count_species_column == 1) M_AD(m,l) = zd_tot**2 * b_unit_y(i) * b_unit_x(i) / q(i,irho) 
            if (count_species_column == 3) M_AD(m,l) = zd_tot**2 * b_unit_y(i) * b_unit_z(i) / q(i,irho)
          endif 
          if (count_species_row == 3) then
            if (count_species_column == 1) M_AD(m,l) = zd_tot**2 * b_unit_z(i) * b_unit_x(i) / q(i,irho)
            if (count_species_column == 2) M_AD(m,l) = zd_tot**2 * b_unit_z(i) * b_unit_y(i) / q(i,irho)
          endif
        endif

        if (index_species_row == 0 .and. index_species_column /= 0) then
          if (count_species_row == 1 .and. count_species_column == 1) M_AD(m,l) = zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (b_unit_y(i)**2 + b_unit_z(i)**2) / q(i,irho)
          if (count_species_row == 2 .and. count_species_column == 2) M_AD(m,l) = zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (b_unit_z(i)**2 + b_unit_x(i)**2) / q(i,irho)
          if (count_species_row == 3 .and. count_species_column == 3) M_AD(m,l) = zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (b_unit_x(i)**2 + b_unit_y(i)**2) / q(i,irho)

          if (count_species_row == 1) then
            if (count_species_column == 2) M_AD(m,l) = - zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_x(i) * b_unit_y(i) / q(i,irho) 
            if (count_species_column == 3) M_AD(m,l) = - zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_x(i) * b_unit_z(i) / q(i,irho) 
          endif 
          if (count_species_row == 2) then
            if (count_species_column == 1) M_AD(m,l) = - zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_y(i) * b_unit_x(i) / q(i,irho)  
            if (count_species_column == 3) M_AD(m,l) = - zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_y(i) * b_unit_z(i) / q(i,irho) 
          endif
          if (count_species_row == 3) then
            if (count_species_column == 1) M_AD(m,l) = - zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_z(i) * b_unit_x(i) / q(i,irho)  
            if (count_species_column == 2) M_AD(m,l) = - zd_tot * (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * b_unit_z(i) * b_unit_y(i) / q(i,irho) 
          endif
        endif

        if (index_species_row /= 0 .and. index_species_column == 0) then
          if (count_species_row == 1 .and. count_species_column == 1) M_AD(m,l) = - zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * (b_unit_y(i)**2 + b_unit_z(i)**2) / q(i,irhod(index_species_row))
          if (count_species_row == 2 .and. count_species_column == 2) M_AD(m,l) = - zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * (b_unit_z(i)**2 + b_unit_x(i)**2) / q(i,irhod(index_species_row))
          if (count_species_row == 3 .and. count_species_column == 3) M_AD(m,l) = - zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * (b_unit_x(i)**2 + b_unit_y(i)**2) / q(i,irhod(index_species_row))

          if (count_species_row == 1) then
            if (count_species_column == 2) M_AD(m,l) = zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_x(i) * b_unit_y(i) / q(i,irhod(index_species_row)) 
            if (count_species_column == 3) M_AD(m,l) = zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_x(i) * b_unit_z(i) / q(i,irhod(index_species_row)) 
          endif 
          if (count_species_row == 2) then
            if (count_species_column == 1) M_AD(m,l) = zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_y(i) * b_unit_x(i) / q(i,irhod(index_species_row))  
            if (count_species_column == 3) M_AD(m,l) = zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_y(i) * b_unit_z(i) / q(i,irhod(index_species_row)) 
          endif
          if (count_species_row == 3) then
            if (count_species_column == 1) M_AD(m,l) = zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_z(i) * b_unit_x(i) / q(i,irhod(index_species_row))  
            if (count_species_column == 2) M_AD(m,l) = zd_tot * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_z(i) * b_unit_y(i) / q(i,irhod(index_species_row)) 
          endif
        endif

        if (index_species_row /= 0 .and. index_species_column /= 0) then
          if (count_species_row == 1 .and. count_species_column == 1) M_AD(m,l) = (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * (b_unit_y(i)**2 + b_unit_z(i)**2) / q(i,irhod(index_species_row))
          if (count_species_row == 2 .and. count_species_column == 2) M_AD(m,l) = (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * (b_unit_z(i)**2 + b_unit_x(i)**2) / q(i,irhod(index_species_row))
          if (count_species_row == 3 .and. count_species_column == 3) M_AD(m,l) = (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * (b_unit_x(i)**2 + b_unit_y(i)**2) / q(i,irhod(index_species_row))

          if (count_species_row == 1) then
            if (count_species_column == 2) M_AD(m,l) = - (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_x(i) * b_unit_y(i) / q(i,irhod(index_species_row)) 
            if (count_species_column == 3) M_AD(m,l) = - (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_x(i) * b_unit_z(i) / q(i,irhod(index_species_row)) 
          endif 
          if (count_species_row == 2) then
            if (count_species_column == 1) M_AD(m,l) = - (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_y(i) * b_unit_x(i) / q(i,irhod(index_species_row))
            if (count_species_column == 3) M_AD(m,l) = - (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_y(i) * b_unit_z(i) / q(i,irhod(index_species_row))
          endif
          if (count_species_row == 3) then
            if (count_species_column == 1) M_AD(m,l) = - (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_z(i) * b_unit_x(i) / q(i,irhod(index_species_row))  
            if (count_species_column == 2) M_AD(m,l) = - (q(i,irhod(index_species_column))/mdust(i,index_species_column) * zd(i,index_species_column) * e_el_stat) * (q(i,irhod(index_species_row))/mdust(i,index_species_row) * zd(i,index_species_row) * e_el_stat) * b_unit_z(i) * b_unit_y(i) / q(i,irhod(index_species_row))
          endif
        endif

        count_species_column = count_species_column + 1

      end do

      count_species_row = count_species_row + 1

    end do

    M_AD(:,:) = eta_a(i) * M_AD(:,:) 




    ! print*, 'M_AD(1,:)=',M_AD(1,:)
    ! print*, 'M_AD(2,:)=',M_AD(2,:)
    ! print*, 'M_AD(3,:)=',M_AD(3,:)
    ! print*, 'M_AD(4,:)=',M_AD(4,:)
    ! print*, 'M_AD(5,:)=',M_AD(5,:)
    ! print*, 'M_AD(6,:)=',M_AD(6,:)
    ! print*, 'M_AD(7,:)=',M_AD(7,:)
    ! print*, 'M_AD(8,:)=',M_AD(8,:)
    ! print*, 'M_AD(9,:)=',M_AD(9,:)



    ! print*,'q(i,irho)',q(i,irho)
    ! print*,'q(i,irhod(:))',q(i,irhod(:))
    ! print*,'q,(i,ivx)',q(i,ivx)
    ! print*,'q,(i,ivy)',q(i,ivy)
    ! print*,'q,(i,ivz)',q(i,ivz)
    ! print*,'q,(i,ivdx(1))',q(i,ivdx(1))
    ! print*,'q,(i,ivdy(1))',q(i,ivdy(1))
    ! print*,'q,(i,ivdz(1))',q(i,ivdz(1))
    ! print*,'q,(i,ivdx(2))',q(i,ivdx(2))
    ! print*,'q,(i,ivdy(2))',q(i,ivdy(2))
    ! print*,'q,(i,ivdz(2))',q(i,ivdz(2))

    ! print*, 'V(:)=',V(:)
    ! print*,'dt=',dt


    !Solve linear system 
    call LU_factorization_resolution(Identity_matrix - M_AD * dt, n, n, V,n,nrhs,output_lin_system_AD) 


    !print*, 'output_lin_system_AD(:)=',output_lin_system_AD(:)

    ! ============================
    ! Update conservative variable
    ! ============================
    k = 4
    do idust=1,ndust


      q(i,ivdx(idust)) = output_lin_system_AD(k)
      q(i,ivdy(idust)) = output_lin_system_AD(k+1)
      q(i,ivdz(idust)) = output_lin_system_AD(k+2)

      u_prim(i,ivdx(idust)) = u_prim(i,irhod(idust)) * output_lin_system_AD(k)
      u_prim(i,ivdy(idust)) = u_prim(i,irhod(idust)) * output_lin_system_AD(k+1)
      u_prim(i,ivdz(idust)) = u_prim(i,irhod(idust)) * output_lin_system_AD(k+2)

      ! print*,'output_lin_system_AD(k)-q(i,ivdx(idust))',output_lin_system_AD(k)-q(i,ivdx(idust)) 
      ! print*,'output_lin_system_AD(k+1)-q(i,ivdy(idust))',output_lin_system_AD(k+1)-q(i,ivdy(idust)) 
      ! print*,'output_lin_system_AD(k+2)-q(i,ivdz(idust))',output_lin_system_AD(k+2)-q(i,ivdz(idust)) 


      k = k + 3
      !print*,'k',k
    enddo

    !Now the gas

    q(i,ivx) = output_lin_system_AD(1)
    q(i,ivy) = output_lin_system_AD(2)
    q(i,ivz) = output_lin_system_AD(3)

    u_prim(i,ivx) = u_prim(i,irho) * output_lin_system_AD(1)
    u_prim(i,ivy) = u_prim(i,irho) * output_lin_system_AD(2)
    u_prim(i,ivz) = u_prim(i,irho) * output_lin_system_AD(3)   

    ! print*,'output_lin_system_AD(1)-q(i,ivx)',output_lin_system_AD(1)-q(i,ivx) 
    ! print*,'output_lin_system_AD(2)-q(i,ivy)',output_lin_system_AD(2)-q(i,ivy) 
    ! print*,'output_lin_system_AD(3)-q(i,ivz)',output_lin_system_AD(3)-q(i,ivz) 





  deallocate(M_AD)
  deallocate(output_lin_system_AD)
  deallocate(V)




end subroutine AD_drag


subroutine magnetic_drag

!!!For each term, both q & u_prim have to be updated. The former so that the next term called can work with the updated velocity. The latter so that the CFL works with the updated velocity!!!
!!N.B.: Updating q(:,ivd) is not a problem since those are local source terms. There is no communication with neighboring cells and thus no need for buffer values!!

  use commons
  use parameters


  implicit none 
  integer :: i




            do i=1,ncells
              !if(active_cell(i)==1) then


                call gyro_drift(i) !First to be called, because is very stiff for small grains.


                !!!Then, get a crude estimate of the stiffness of each term by looking at the corresponding resistivity!!!
                !!!Solve for stiffest term first!!!


                if (eta_o(i) .ge. abs(eta_H(i)) .and. abs(eta_H(i)) .ge. eta_a(i)) then 
                  call Ohmic_drag(i) 
                  call Hall_gyromotion(i) 
                  call AD_drag(i) 
                endif

                if (eta_o(i) .ge. eta_a(i) .and. eta_a(i) .ge. abs(eta_H(i))) then 
                  call Ohmic_drag(i)  
                  call AD_drag(i)
                  call Hall_gyromotion(i)
                endif

                if (eta_a(i) .ge. eta_o(i) .and. eta_o(i) .ge. abs(eta_H(i))) then 
                  call AD_drag(i)
                  call Ohmic_drag(i)  
                  call Hall_gyromotion(i)
                endif

                if (eta_a(i) .ge. abs(eta_H(i)) .and. abs(eta_H(i)) .ge. eta_o(i)) then   
                  call AD_drag(i)
                  call Hall_gyromotion(i)
                  call Ohmic_drag(i)
                endif

                if (abs(eta_H(i)) .ge. eta_o(i) .and. eta_o(i) .ge. eta_a(i)) then 
                  call Hall_gyromotion(i)
                  call Ohmic_drag(i)  
                  call AD_drag(i)
                endif


                if (abs(eta_H(i)) .ge. eta_a(i) .and. eta_a(i) .ge. eta_o(i)) then    
                  call Hall_gyromotion(i)
                  call AD_drag(i)
                  call Ohmic_drag(i)
                endif

              !endif
            end do

end subroutine magnetic_drag
#endif
#endif
module lapack_tools

   use precision
   implicit none

contains
subroutine LU_factorization_resolution(M, M_rows, M_columns, B,lda,nrhs,output_vector)

    ! Resolution of linear system of the form MX = B
    ! 1. LU factorization of matrix M. 2. Actual resolution of the simplified (factorized) system.

    ! Input:
    !   M:  the actual matrix
    !   M_rows: nb of rows
    !   M_colums: nb of colums
    !   B: right hand side vector
    !   lda: leading dimension of M: LDA >= max(1,M_rows)
    !   nrhs: The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
    ! DGETRF Returns:
    !   pivot_vector: The pivot indices 
    !   info: !>  INFO is INTEGER
    !          = 0:  successful exit
    !          < 0:  if INFO = -i, the i-th argument had an illegal value
    !          > 0:  if INFO = i, M(i,i) is exactly zero. The factorization
    !                has been completed, but the factor M is exactly
    !                singular, and division by zero will occur if it is used
    !                to solve a system of equations.
    
    !DGETRS Returns:
    !   X is computed and overwitten on B (confusing). So B is in the end the output array that we need after resolution. 
     

    implicit none
    integer, intent(in) :: M_rows
    integer, intent(in) :: M_columns
    integer, intent(in) :: lda
    integer, intent(in) :: nrhs
    real(dp), dimension(1:M_rows,1:M_columns), intent(in) :: M
    real(dp), dimension(1:M_rows), intent(in) :: B
    real(dp), dimension(1:M_rows), intent(inout) :: output_vector


    !Local variable
    integer :: info
    integer, dimension(1:M_rows) :: pivot_vector



    ! ============================================================
    ! LU factorization with pivots (DGETRF)
    ! ============================================================

    call DGETRF(M_rows,M_columns,M,lda,pivot_vector,info)
    !M is now L=U factorized

    if (info /= 0) then
        print *, 'ERROR en DGETRF: info = ', info
        if (info < 0) then
            print *, 'Argumento ', -info, ' tiene un valor ilegal'
        else if (info > 0) then
            print *, 'La matriz es singular (U(', info, ',', info, ') = 0)'
        end if
        stop
    endif

    ! ============================================================
    ! Actual resolution of the system
    ! ============================================================
    call DGETRS('N', lda, nrhs, M, lda, pivot_vector, B, lda, info)

    !Store it in the output_vector
    output_vector(:) = B(:)

    

end subroutine LU_factorization_resolution


end module lapack_tools
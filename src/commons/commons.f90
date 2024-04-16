module commons
  use precision
  use hydro_commons
#if NDUST>0
  use dust_commons
#endif
#if TURB>0
  use turb_commons
#endif
  use gravity_commons
  use ionisation_commons
  use OMP_LIB
  implicit none

  !Grid
  real(dp), dimension(:,:), allocatable :: dx
  real(dp), dimension(:,:), allocatable :: position
  real(dp), dimension(:), allocatable   :: radii
  real(dp), dimension(:), allocatable   :: radii_c

#if GEOM==1
  real(dp), dimension(:), allocatable :: theta
  real(dp), dimension(:), allocatable :: dx_l
  real(dp), dimension(:), allocatable :: dx_r
  real(dp), dimension(:), allocatable :: dx_c
  real(dp), dimension(:), allocatable :: dx_r_cell
  real(dp), dimension(:), allocatable :: dx_l_cell
#endif

#if GEOM==2
  ! Disk face-on geometry
  real(dp), dimension(:), allocatable :: phi
#if GRIDSPACE ==1
  real(dp), dimension(:), allocatable :: dx_r_cell
  real(dp), dimension(:), allocatable :: dx_l_cell
#endif
#endif

  real(dp), dimension(:,:), allocatable :: Surf
  real(dp), dimension(:), allocatable   :: dvol
  real(dp), dimension(:), allocatable   :: vol
  real(dp), dimension(:), allocatable   :: active_cell
  real(dp), dimension(:), allocatable   :: active_cell_predictor
  
  
  real(dp):: time 
  real(dp):: tend
  real(dp):: dt

  contains
  subroutine matrix_vector(A,XX,AX,n)
    use precision
    implicit none
    integer :: n,i,j
    real(dp),dimension(1:n,1:n) :: A
    real(dp),dimension(1:n) :: AX,XX
    AX=0.0d0
    do i=1,n
      do j=1,n
        AX(i) = AX(i) + A(i,j)*xx(j)
      end do
    end do
  end subroutine matrix_vector

  subroutine matrix_matrix(A,B,AB,n)
    use precision
    implicit none
    integer :: n,i,j, k
    real(dp),dimension(1:n,1:n) :: A,B,AB
    AB=0.0d0
    do j=1,n
      do i=1,n
        do k=1,n
          AB(i,j) = AB(i,j) + A(i,k)*B(k,j)
        end do
      end do
    end do
  end subroutine matrix_matrix

 subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
use precision
implicit none 
integer  :: n
real(dp),dimension(1:n,1:n) :: a, c
real(dp),dimension(1:n,1:n) :: L, U
real(dp),dimension(1:n) :: b, d, x
real(dp) :: coeff
integer :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0d0
U=0.0d0
b=0.0d0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0d0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0d0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0d0
end do
end subroutine inverse
end module commons


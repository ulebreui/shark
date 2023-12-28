! The precision of the code is set here.
MODULE real_precisions
IMPLICIT NONE
! les differentes variantes de reels disponibles (g95/gfortran)
INTEGER, PARAMETER :: sp =SELECTED_REAL_KIND(p=6)  ! simple precision 32 bits
INTEGER, PARAMETER :: dop=kind(1.0D0) ! real*8 ! double precision 64 bits
INTEGER, PARAMETER :: edp=SELECTED_REAL_KIND(p=18) ! extra double precision 80
END MODULE real_precisions

MODULE precision ! ne modifier que ce module pour choisir la precision wp
!use real_precisions, only : dp => sp ! par exemple simple precision
USE real_precisions, ONLY : dp => dop ! ou sinon double precision
!use real_precisions, only : dp => edp ! ou sinon extra double precision
implicit none

END MODULE precision

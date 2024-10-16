module newton_raphson

   use precision
   implicit none

contains
subroutine solve_newton_raphson(f, fp, y0, y, iters, maxiter, debug)

    ! Estimate the zero of f(y) using Newton's method. 
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   y0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate y satisfying f(y)=0 (assumes Newton converged!) 
    !   the number of iterations iters
     
    implicit none
    real(dp), intent(in) :: y0
    real(dp), external :: f, fp
    logical, intent(in) :: debug
    real(dp), intent(out) :: y
    integer, intent(inout) :: iters
    integer, intent(in) :: maxiter
    real(dp)  :: tol = 1.d-14



    ! Declare any local variables:
    real(dp) :: deltay, fy, fyprime
    integer :: k


    ! initial guess
    y = y0

    if (debug) then
        print 11, y
 11     format('Initial guess: y = ', e22.15)
        endif

    ! Newton iteration to find a zero of f(y) 

    do k=1,maxiter

        ! evaluate function and its derivative:
        fy = f(y)
        fyprime = fp(y)

        if (debug) then
            print 12, fy
 12         format('fy = ', e22.15)
            endif

        if (debug) then
            print 13, fyprime
 13         format('dfy/dy = ', e22.15)
            endif


        if (abs(fy) < tol) then
            exit  ! jump out of do loop
            endif

        ! compute Newton increment y:
        deltay = fy/fyprime
        if (debug) then
            print 14, deltay
 14         format('Increment: deltay = ', e22.15)
            endif

        ! update y:
        y = y - deltay

        if (debug) then
            print 15, k,y
 15         format('After', i3, ' iterations, y = ', e22.15)
            endif

        enddo


    if (k > maxiter) then
        ! might not have converged

        fy = f(y)
        if (abs(fy) > tol) then
            print *, '*** Warning: has not yet converged'
            endif
        endif 

    ! number of iterations taken:
    iters = k-1


end subroutine solve_newton_raphson
end module newton_raphson
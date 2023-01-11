double precision function B_field(rho)
   use commons
   use units, ONLY : unit_nH
   implicit none

   real(dp) :: rho

   B_field=min(B_0_lee*sqrt(rho*unit_nh/1d4),0.1d0)

 end function B_field

! 'self-consistent' B field. This is a very crude approach, for now !
 subroutine compute_B
   use commons
   use units
   implicit none
   
   real(dp) :: B_field,tff,tadd
   integer:: ii,i
   
   !External B_field
   B_cell(last_active)=B_field(uold(last_active,irho)*unit_nh)
   do i=first_active,last_active-1
      ii=ncells-i
      B_cell(ii)=B_cell(ii+1)*sqrt(uold(ii,irho)/uold(ii+1,irho))


      if(eta_A(ii).ne.0.) then
         tff=sqrt(3.*pi/(32.*grav*uold(ii,irho)*unit_d))
         tadd=4.*pi/(clight**2*eta_A(ii))*(r_c(ii)*unit_l)**2.
         if(tff>tadd) B_cell(ii)=B_cell(ii+1)
      endif
   end do

 end subroutine compute_B

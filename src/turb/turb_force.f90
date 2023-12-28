! ! Generate turbulent force /!\ Not working yet 
! subroutine turbulent_force
!   use parameters
!   use commons
!   use OMP_LIB
!   use units
!   implicit none
!   integer :: ix,iy,icell,k_mode,idim
!   real(dp):: rand_number,xx,yy, intensity
!   f_k_old = f_k
!   do k_mode = 1, n_modes
!       call random_number(intensity)
!       f_k(k_mode)       =  2.0d0*(intensity-1.0d0) * (v_rms/unit_v)/ turb_turnover
!       call random_number(intensity)
!       phase_turb(k_mode) = intensity
!   end do
!   !f_k = f_k - dt/turb_correl * f_k_old ! Substract decay

!   do iy = first_active_y,last_active_y
!     do ix = first_active,last_active
!       xx=position(icell(ix,iy),1)-half*box_l  ! Boxlen already in pc
!       a_turb(icell(ix,iy),1) = 0.0d0
! #if NY>1
!       yy=position(icell(ix,iy),2)-half*box_l_y
! #endif    
!       a_turb(icell(ix,iy),2) = 0.0d0
!       do k_mode=1,n_modes
!         a_turb(icell(ix,iy),1)= a_turb(icell(ix,iy),1) + f_k(k_mode) * cos(2.0d0*pi*real(k_mode)*xx/box_l) 
! #if NY>1
!         a_turb(icell(ix,iy),2)= a_turb(icell(ix,iy),2) + f_k(k_mode) * cos(2.0d0*pi*real(k_mode)*yy/box_l_y)
! #endif        
!       end do
!    end do
!   end do
! end subroutine turbulent_force


! Generate turbulent force 
subroutine turbulent_force
  use parameters
  use commons
  use OMP_LIB
  use units
  implicit none
  integer :: ix,iy,icell,k_mode,idim
  real(dp):: rand_number,xx,yy, intensity
  f_k_old = f_k
  do k_mode = 1, n_modes
      call random_number(intensity)
      f_k(k_mode)       =  2.0d0*(intensity-1.0d0) * (v_rms/unit_v)/ turb_turnover
      call random_number(intensity)
      phase_turb(k_mode) = intensity
  end do
  !f_k = f_k - dt/turb_correl * f_k_old ! Substract decay

  do iy = first_active_y,last_active_y
    do ix = first_active,last_active
      xx=position(icell(ix,iy),1)-half*box_l  ! Boxlen already in pc
      a_turb(icell(ix,iy),1) = 0.0d0
#if NY>1
      yy=position(icell(ix,iy),2)-half*box_l_y
#endif    
      a_turb(icell(ix,iy),2) = 0.0d0
      do k_mode=1,n_modes
        a_turb(icell(ix,iy),1)= a_turb(icell(ix,iy),1) + f_k(k_mode) * cos(2.0d0*pi*real(k_mode)*xx/box_l) 
#if NY>1
        a_turb(icell(ix,iy),2)= a_turb(icell(ix,iy),2) + f_k(k_mode) * cos(2.0d0*pi*real(k_mode)*yy/box_l_y)
#endif        
      end do
   end do
  end do
end subroutine turbulent_force


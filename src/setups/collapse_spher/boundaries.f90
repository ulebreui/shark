!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This applies the boundaries either to uold or unew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine apply_boundaries(who_app,uu,nn,nn2)
  use parameters
  use commons
  use units
  implicit none
  integer :: who_app,idust,ighost,nn,nn2,i
  real(dp), dimension (1:nn,1:nn2) :: uu
  if(who_app .eq. 1) then
     uu=uold
  else if(who_app .eq. 2) then
     uu=unew
  endif
  do ighost=1,nghost
     uu(ighost,:)= uu(nghost+1,:)
     uu(ncells-ighost+1,:)=uu(ncells-nghost,:)

     uu(ncells-ighost+1,iv)= 0.0d0!max(uu(ncells-nghost,iv),0.d0)
     uu(ighost,iv)=0.0d0
     
#if NDUST>0
     do idust=1,ndust
        uu(ncells-ighost+1,ivd(idust))=0.0d0!max(uu(ncells-nghost,ivd(idust)),0.d0)
        uu(ighost,ivd(idust))=0.0d0
      enddo
#endif 
  end do
  uu(ncells-nghost,iv)= 0.0d0!max(uu(ncells-nghost,iv),0.d0)
  uu(nghost+1,iv)=0.0d0
#if NDUST>0
     do idust=1,ndust
        uu(ncells-nghost,ivd(idust))=0.0d0!max(uu(ncells-nghost,ivd(idust)),0.d0)
        uu(nghost+1,ivd(idust))=0.0d0
     enddo
#endif
#if MHD>0
!!$     do ighost=1,nghost
!!$        uu(ncells-ighost+1,iB)=B_ext ! Imposed outer boundary for B
!!$     end do
#endif
  if(who_app .eq. 1) then
     uold=uu
  else if(who_app .eq. 2) then
     unew=uu
  endif
end subroutine apply_boundaries



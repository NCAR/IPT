#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module mass_matrix_mod

  ! Useful modules
  !---------------
  use SE_Constants  ,only: real_kind,np
  use quadrature_mod,only: quadrature_t, gauss ,gausslobatto
  use element_mod   ,only: element_t
  use edge_mod      ,only: edgebuffer_t,edgevpack,edgevunpack,freeedgebuffer,initedgebuffer  

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private

  public:: mass_matrix


  ! PFC: Hack to set global option values
  !--------------------------------------
  public :: init_mass_matrix_mod
  private:: nelemd
  integer:: nelemd


contains
  !==================================================================
  subroutine init_mass_matrix_mod(I_SEopt)
    ! PFC: Hack to set global option values
    !--------------------------------------
    use SE_Options,only: SEoptions_t
    ! Passed Variables
    !------------------
    type(SEoptions_t),intent(in):: I_SEopt

    nelemd          = I_SEopt%nelemd

   print *,' PFC: nelemd=',nelemd,I_SEopt%nelemd

    ! End Routine
    !-------------
    return
  end subroutine init_mass_matrix_mod
  !==================================================================


  !=============================================
  subroutine mass_matrix(elem)
    ! mass_matrix:
    !
    ! Compute the mass matrix for each element...
    !===============================================
    !
    ! Passed Variables
    !-------------------
    type(element_t )           :: elem(:)
    !
    ! Local Values
    !---------------
    type(EdgeBuffer_t)  :: edge
    real(kind=real_kind):: da        ! area element
    type(quadrature_t)  :: gp

    integer ie
    integer ii,jj
    integer kptr
    integer iptr

    ! Init edge buffer for boundary echecnge
    !----------------------------------------
    call initEdgeBuffer(edge,1)

    ! mass matrix on the velocity grid
    !----------------------------------
    gp=gausslobatto(np)
 
    do ie=1,nelemd
      do jj=1,np
      do ii=1,np
        ! MNL: metric term for map to reference element is now in metdet!
        !----------------------------------------------------------------
        elem(ie)%mp (ii,jj)=gp%weights(ii)*gp%weights(jj)
        elem(ie)%rmp(ii,jj)=elem(ie)%mp(ii,jj)
      end do
      end do
      kptr=0
      call edgeVpack(edge,elem(ie)%rmp,1,kptr,elem(ie)%desc)
    end do

    ! Unpack and finish rmp calc
    !---------------------------------
    do ie=1,nelemd
      kptr=0
      call edgeVunpack(edge,elem(ie)%rmp,1,kptr,elem(ie)%desc)
      do jj=1,np
      do ii=1,np
        elem(ie)%rmp(ii,jj)=1.0D0/elem(ie)%rmp(ii,jj)
      end do
      end do
    end do
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif

    ! Clean up
    !-----------
    deallocate(gp%points)
    deallocate(gp%weights)

    ! compute spherical element mass matrix
    !--------------------------------------
    do ie=1,nelemd
      do jj=1,np
      do ii=1,np
        elem(ie)%spheremp (ii,jj)=elem(ie)%mp(ii,jj)*elem(ie)%metdet(ii,jj)
        elem(ie)%rspheremp(ii,jj)=elem(ie)%spheremp(ii,jj)
      end do
      end do
      kptr=0
      call edgeVpack(edge,elem(ie)%rspheremp,1,kptr,elem(ie)%desc)
    end do

    do ie=1,nelemd
      kptr=0
      call edgeVunpack(edge,elem(ie)%rspheremp,1,kptr,elem(ie)%desc)
      do jj=1,np
      do ii=1,np
        elem(ie)%rspheremp(ii,jj)=1.0D0/elem(ie)%rspheremp(ii,jj)
      end do
      end do
    end do
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif

    ! compute the mass matrix 
    !=============================================
    ! Jose Garcia: Not sure but I think this code is just dead code
    !do ie=1,nelemd
    !  iptr=1
    !  do j=1,np
    !  do i=1,np
    !    elem(ie)%mp(i,j)=elem(ie)%mp(i,j)
    !    iptr=iptr+1
    !  end do
    !  end do
    !end do
   
    ! All don with the edge buffer
    !-----------------------------
    call FreeEdgeBuffer(edge)

    ! End Routine
    !--------------
    return
  end subroutine mass_matrix
  !=============================================

end module mass_matrix_mod

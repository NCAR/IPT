#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dof_mod

  ! Useful modules
  !----------------
  use SE_Constants  ,only: real_kind,int_kind,long_kind,np,npsq
  use quadrature_mod,only: quadrature_t
  use element_mod   ,only: element_t,index_t
  use edge_mod      ,only: longedgebuffer_t,initlongedgebuffer,freelongedgebuffer, &
                           longedgevpack   ,longedgevunpackmin

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private

  public :: genLocalDof
  public :: global_dof

  public :: UniquePoints
  public :: PutUniquePoints
  private:: UniquePoints2D
  private:: PutUniquePoints2D
  public :: UniqueNcolsP
  public :: UniqueCoords
  private:: UniquePoints3D
  private:: UniquePoints4D
  private:: PutUniquePoints3D
  private:: PutUniquePoints4D
  public :: SetElemOffset
  public :: CreateUniqueIndex
  public :: CreateMetaData
  public :: PrintDofP

  ! PFC: Hack to set global option values
  !--------------------------------------
  public :: init_dof_mod
  private:: nelem
  private:: nelemd
  integer:: nelem
  integer:: nelemd

  ! Public interface 
  !------------------------
  interface UniquePoints
     module procedure UniquePoints2D
     module procedure UniquePoints3D
     module procedure UniquePoints4D
  end interface
  interface PutUniquePoints
     module procedure PutUniquePoints2D
     module procedure PutUniquePoints3D
     module procedure PutUniquePoints4D
  end interface


contains
  !==================================================================
  subroutine init_dof_mod(I_SEopt)
    ! PFC: Hack to set global option values
    !--------------------------------------
    use SE_Options,only: SEoptions_t
    ! Passed Variables
    !------------------
    type(SEoptions_t),intent(in):: I_SEopt

    nelem           = I_SEopt%nelem
    nelemd          = I_SEopt%nelemd

    ! End Routine
    !-------------
    return
  end subroutine init_dof_mod
  !==================================================================


  !==================================================================
  subroutine genLocalDof(ig,npts,ldof)
    ! genLocalDof:
    !
    !===============================================================
    !
    ! Passed Variables
    !-------------------
    integer(kind=int_kind),intent(in   ):: ig
    integer(kind=int_kind),intent(in   ):: npts
    integer(kind=int_kind),intent(inout):: ldof(:,:)
    !
    ! Local Values
    !---------------
    integer(kind=int_kind):: ii,jj,npts2
   
    npts2=npts*npts
    do jj=1,npts
    do ii=1,npts
      ldof(ii,jj) = (ig-1)*npts2 + (jj-1)*npts + ii
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine genLocalDOF
  !==================================================================


  !==================================================================
  subroutine global_dof(elem)
    ! global_dof: Compute the global degree of freedom for each element...
    !==================================================================
    !
    ! Passed Variables
    !-------------------
    type(element_t)            :: elem(:)
    !
    ! Local Values
    !----------------
    type(LongEdgeBuffer_t):: edge
    real(kind=real_kind)  :: da                     ! area element
    type(quadrature_t)    :: gp
    integer(kind=int_kind):: ldofP(np,np,nelemd)

    integer ig,ie
    integer kptr

    ! begin code
    !-------------
    call initLongEdgeBuffer(edge,1)

    ! mass matrix on the velocity grid
    !----------------------------------
    do ie=1,nelemd
      ig  = elem(ie)%vertex%number
      kptr= 0
      call genLocalDOF(ig,np,ldofP(:,:,ie))
      call LongEdgeVpack(edge,ldofP(:,:,ie),1,kptr,elem(ie)%desc)
    end do

    do ie=1,nelemd
      ! we should unpack directly into elem(ie)%gdofV, but we dont have
      ! a VunpackMIN that takes integer*8.  gdofV integer*8 means  
      ! more than 2G grid points.
      !---------------------------------------------------------------
      kptr=0
      call LongEdgeVunpackMIN(edge,ldofP(:,:,ie),1,kptr,elem(ie)%desc)
      elem(ie)%gdofP(:,:)=ldofP(:,:,ie)
    end do

#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif

    call FreeLongEdgeBuffer(edge)
       
    ! End Routine
    !-------------
    return
  end subroutine global_dof
  !==================================================================


  !==================================================================
  subroutine UniquePoints2D(idxUnique,src,dest)
    ! UniquePoints2D:
    !
    !===============================================
    !
    ! Passed Variables
    !-------------------
    type(index_t)       :: idxUnique
    real(kind=real_kind):: src(:,:)
    real(kind=real_kind):: dest(:)
    !
    ! Local Values
    !--------------
    integer(kind=int_kind):: ii,jj,nn
    
    do nn=1,idxUnique%NumUniquePts
      ii      = idxUnique%ia(nn)
      jj      = idxUnique%ja(nn)
      dest(nn)= src(ii,jj)
    end do

    ! End Routine
    !-------------
    return
  end subroutine UniquePoints2D
  !==================================================================


  !==================================================================
  subroutine PutUniquePoints2D(idxUnique,src,dest)
    ! PutUniquePoints:  first zeros out the destination array, then fills the 
    !                   unique points of the array with values from src.  
    !                   A boundary communication should then be called to fill 
    !                   in the redundent points of the array
    !=======================================================================
    !
    ! Passed Variables
    !---------------------
    type(index_t)                   :: idxUnique
    real(kind=real_kind),intent(in ):: src(:)
    real(kind=real_kind),intent(out):: dest(:,:)
    !
    ! Local Values
    !---------------
    integer(kind=int_kind):: ii,jj,nn
    
    dest=0.0D0
    do nn=1,idxUnique%NumUniquePts
      ii         = idxUnique%ia(nn)
      jj         = idxUnique%ja(nn)
      dest(ii,jj)= src(nn)
    end do

    ! End Routine
    !-------------
    return
  end subroutine PutUniquePoints2D
  !==================================================================


  !==================================================================
  subroutine UniqueNcolsP(elem,idxUnique,cid)    
    ! UniqueNcolsP:
    !
    !===============================================================
    use element_mod,only: GetColumnIdP, element_t
    !
    ! Passed Variables
    !------------------
    type(element_t),intent(in ):: elem
    type(index_t)  ,intent(in ):: idxUnique
    integer        ,intent(out):: cid(:)
    !
    ! Local Values
    !---------------
    integer(kind=int_kind):: ii,jj,nn

    do nn=1,idxUnique%NumUniquePts
      ii     = idxUnique%ia(nn)
      jj     = idxUnique%ja(nn)
      cid(nn)= GetColumnIdP(elem,ii,jj)
    end do
    
    ! End Routine
    !-------------
    return
  end subroutine UniqueNcolsP
  !==================================================================


  !==================================================================
  subroutine UniqueCoords(idxUnique,src,lat,lon)
    ! UniqueCoords:
    !
    !=============================================================
    use coordinate_systems_mod,only: spherical_polar_t
    !
    ! Passed Variables
    !--------------------
    type(index_t)          ,intent(in ):: idxUnique
    type(spherical_polar_t)            :: src(:,:)
    real(kind=real_kind)   ,intent(out):: lat(:)
    real(kind=real_kind)   ,intent(out):: lon(:)
    !
    ! Local Values
    !----------------
    integer(kind=int_kind):: ii,jj,nn

    do nn=1,idxUnique%NumUniquePts
      ii     = idxUnique%ia(nn)
      jj     = idxUnique%ja(nn)
      lat(nn)= src(ii,jj)%lat
      lon(nn)= src(ii,jj)%lon
    end do

    ! End Routine
    !-------------
    return
  end subroutine UniqueCoords
  !==================================================================


  !==================================================================
  subroutine UniquePoints3D(idxUnique,nlyr,src,dest)
    ! UniquePoints3D:
    !
    !=============================================================
    !
    ! Passed Variables
    !------------------
    type(index_t)         :: idxUnique
    integer(kind=int_kind):: nlyr
    real(kind=real_kind)  :: src(:,:,:)
    real(kind=real_kind)  :: dest(:,:)
    !
    ! Local Values
    !--------------
    integer(kind=int_kind) :: ii,jj,kk,nn

    do nn=1,idxUnique%NumUniquePts
      ii = idxUnique%ia(nn)
      jj = idxUnique%ja(nn)
      do kk=1,nlyr
        dest(nn,kk) = src(ii,jj,kk)
      end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine UniquePoints3D
  !==================================================================


  !==================================================================
  subroutine UniquePoints4D(idxUnique,d3,d4,src,dest)
    ! UniquePoints4D:
    !
    !=========================================================
    !
    ! Passed Variables
    !--------------------
    type(index_t)         :: idxUnique
    integer(kind=int_kind):: d3,d4
    real(kind=real_kind)  :: src(:,:,:,:)
    real(kind=real_kind)  :: dest(:,:,:)
    !
    ! Local Values
    !--------------
    integer(kind=int_kind) :: ii,jj,kk,ll,nn

    do ll=1,d4
    do kk=1,d3
      do nn=1,idxUnique%NumUniquePts
        ii            = idxUnique%ia(nn)
        jj            = idxUnique%ja(nn)
        dest(nn,kk,ll)= src(ii,jj,kk,ll)
      end do
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine UniquePoints4D
  !==================================================================


  !==================================================================
  subroutine PutUniquePoints3D(idxUnique,nlyr,src,dest)
    ! PutUniquePoints:  first zeros out the destination array, then fills the 
    !                   unique points of the array with values from src.  
    !                   A boundary communication should then be called to fill 
    !                   in the redundent points of the array
    !==========================================================================
    !
    ! Passed Variables
    !--------------------
    type(index_t)                     :: idxUnique
    integer(kind=int_kind)            :: nlyr
    real(kind=real_kind)  ,intent(in ):: src(:,:)
    real(kind=real_kind)  ,intent(out):: dest(:,:,:)
    !
    ! Local Values
    !----------------
    integer(kind=int_kind) :: ii,jj,kk,nn

    dest=0.0D0
    do kk=1,nlyr
      do nn=1,idxUnique%NumUniquePts
        ii            = idxUnique%ia(nn)
        jj            = idxUnique%ja(nn)
        dest(ii,jj,kk)= src(nn,kk)
      end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine PutUniquePoints3D
  !==================================================================


  !==================================================================
  subroutine PutUniquePoints4D(idxUnique,d3,d4,src,dest)
    ! PutUniquePoints4D:
    !
    !================================================================
    !
    ! Passed Variables
    !--------------------
    type(index_t)                     :: idxUnique
    integer(kind=int_kind)            :: d3,d4
    real(kind=real_kind)  ,intent(in ):: src(:,:,:)
    real(kind=real_kind)  ,intent(out):: dest(:,:,:,:)
    !
    ! Local Values
    !----------------
    integer(kind=int_kind) :: ii,jj,kk,ll,nn

    dest=0.0D0
    do ll=1,d4
    do kk=1,d3
      do nn=1,idxunique%NumUniquePts
        ii               = idxUnique%ia(nn)
        jj               = idxUnique%ja(nn)
        dest(ii,jj,kk,ll)= src(nn,kk,ll)
      end do
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine PutUniquePoints4D
  !==================================================================


  !==================================================================
  subroutine SetElemOffset(elem,GlobalUniqueColsP)
    ! SetElemOffset:
    !
    !================================================================
     !
     ! Passed Variables
     !------------------------
     type(element_t)    :: elem(:)
     integer,intent(out):: GlobalUniqueColsP
     !
     ! Local Values
     !--------------
     logical,parameter                 :: Debug = .FALSE.
     integer(kind=int_kind),allocatable:: numElemP(:),numElem2P(:)
     integer(kind=int_kind),allocatable:: numElemV(:),numElem2V(:)
     integer(kind=int_kind),allocatable:: gOffset(:)
     integer(kind=int_kind)            :: ie,ig,ierr

     allocate(numElemP (nelem))
     allocate(numElem2P(nelem))
     allocate(gOffset  (nelem))
     numElemP (:)= 0
     numElem2P(:)= 0
     gOffset  (:)= 0

     do ie=1,nelemd
       ig          = elem(ie)%GlobalId
       numElemP(ig)= elem(ie)%idxP%NumUniquePts
     end do

     numElem2P(:) = numElemP(:)

     gOffset(1)=1
     do ig=2,nelem
       gOffset(ig) = gOffset(ig-1)+numElem2P(ig-1)
     end do

     do ie=1,nelemd
       ig                          = elem(ie)%GlobalId
       elem(ie)%idxP%UniquePtOffset= gOffset(ig)
     enddo

     GlobalUniqueColsP = gOffset(nelem)+numElem2P(nelem)-1

     deallocate(numElemP )
     deallocate(numElem2P)
     deallocate(gOffset  )

    ! End Routine
    !-------------
    return
  end subroutine SetElemOffset
  !==================================================================


  !==================================================================
  subroutine CreateUniqueIndex(ig,gdof,idx)
    ! CreateUniqueIndex:
    !
    !================================================================
    !
    ! Passed Variables
    !---------------------
    integer(kind=int_kind) :: ig
    integer(kind=long_kind):: gdof(:,:)
    type(index_t)          :: idx 
    !
    ! Local Values
    !----------------
    integer,allocatable:: ldof(:,:)
    integer            :: ii,jj,nn,npts

    ! Form the local DOF
    !-----------------------
    npts = size(gdof,dim=1)
    allocate(ldof(npts,npts))
    call genLocalDOF(ig,npts,ldof)
    
    nn=1
    do jj=1,npts
    do ii=1,npts
      ! check for point ownership
      !---------------------------
      if(gdof(ii,jj) .eq. ldof(ii,jj)) then
        idx%ia(nn) = ii
        idx%ja(nn) = jj
        nn=nn+1
      endif
    end do
    end do
    idx%NumUniquePts=nn-1
    
    deallocate(ldof)

    ! End Routine
    !-------------
    return
  end subroutine CreateUniqueIndex
  !==================================================================


  !==================================================================
  subroutine CreateMetaData(elem,subelement_corners, fdofp)
    ! CreateMetaData:
    !
    !==========================================================
    !
    ! Passed Variables
    !------------------
    type(element_t),target         :: elem(:)
    integer   ,optional,intent(out):: subelement_corners((np-1)*(np-1)*nelemd,4)
    integer(kind=int_kind),optional:: fdofp(np,np,nelemd)
    !
    ! Local Values
    !-------------
    type(index_t)          ,pointer:: idx 
    type(LongEdgeBuffer_t)         :: edge
    integer(kind=long_kind),pointer:: gdof(:,:)
    integer                        :: fdofp_local(np,np,nelemd)

    integer  ii,jj,nn,ie,base

    call initLongEdgeBuffer(edge,1)
    fdofp_local=0
    
    do ie=1,nelemd
      idx => elem(ie)%idxP
      do nn=1,idx%NumUniquePts
        ii                    = idx%ia(nn)
        jj                    = idx%ja(nn)
        fdofp_local(ii,jj,ie) = -(idx%UniquePtoffset+nn-1)
      end do
      call LongEdgeVpack(edge,fdofp_local(:,:,ie),1,0,elem(ie)%desc)
    end do

    do ie=1,nelemd
      base = (ie-1)*(np-1)*(np-1)
      call LongEdgeVunpackMIN(edge,fdofp_local(:,:,ie),1,0,elem(ie)%desc)
      if(present(subelement_corners)) then
        nn=0       
        do jj=1,np-1
        do ii=1,np-1
          nn=nn+1
          subelement_corners(base+nn,1) = -fdofp_local(ii  ,jj  ,ie)
          subelement_corners(base+nn,2) = -fdofp_local(ii  ,jj+1,ie)
          subelement_corners(base+nn,3) = -fdofp_local(ii+1,jj+1,ie)
          subelement_corners(base+nn,4) = -fdofp_local(ii+1,jj  ,ie)
        end do
        end do
      endif
    end do

    if(present(fdofp)) then
      fdofp=-fdofp_local
    endif

    ! End Routine
    !-------------
    return
  end subroutine CreateMetaData
  !==================================================================


  !==================================================================
  subroutine PrintDofP(elem)
    !  PrintDofP:    Prints the degree of freedom 
    !=================================================
    !
    ! Passed Variables
    !--------------------
    type(element_t),intent(in):: elem(:)
    !
    ! Local Values
    !---------------
    integer :: ie,nse,ii,jj

    nse = SIZE(elem)
 
    do ie=1,nse
      print *,'Element # ',elem(ie)%vertex%number
      do jj=np,1,-1
        write(6,*) (elem(ie)%gdofP(ii,jj), ii=1,np)
      end do
    end do
 10 format('I5')

    ! End Routine
    !-------------
    return
  end subroutine PrintDofP
  !==================================================================

end module dof_mod

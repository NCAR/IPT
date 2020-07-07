#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module edge_mod

  ! Useful modules
  !----------------
  use SE_Constants  ,only: int_kind, log_kind, real_kind
!PFC  use perf_mod      ,only: t_startf, t_stopf ! _EXTERNAL
  use err_exit      ,only: endrun
  use thread_mod    ,only: omp_get_num_threads, omp_get_thread_num

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private
  save

  public :: rotation_t
  public :: EdgeDescriptor_t
  public :: EdgeBuffer_t
  public :: LongEdgeBuffer_t
  public :: GhostBufferTR_t
  public :: GhostBuffer3D_t
  private:: wrap_ptr

  private:: threadsafe
  private:: edgebuff_ptrs

  public :: initEdgeBuffer
  public :: initLongEdgeBuffer
  public :: edgeDGVpack
  public :: FreeEdgeBuffer
  public :: FreeLongEdgeBuffer
  public :: edgeVpack
  public :: LongEdgeVpack
  public :: edgeVunpack
  public :: edgeVunpackVert
  public :: edgeDGVunpack
  public :: edgeVunpackMAX
  public :: edgeVunpackMIN
  public :: LongEdgeVunpackMIN
  public :: edgerotate
  public :: buffermap


  ! NOTE ON ELEMENT ORIENTATION
  ! for the edge neighbors:  
  !    we set the "reverse" flag if two elements who share an edge use a 
  !    reverse orientation.  The data is reversed during the *pack* stage
  ! For corner neighbors:  
  !    for edge buffers, there is no orientation because two corner neighbors
  !    only share a single point.
  !    For ghost cell data, there is a again two posible orientations. For
  !    this case, we set the "reverse" flag if the corner element is using 
  !    the reverse orientation.  In this case, the data is reversed during the
  !    *unpack* stage (not sure why)
  !
  ! The edge orientation is set at startup.  The corner orientation is computed
  ! at run time, via the call to compute_ghost_corner_orientation()
  ! This routine only works for meshes with at most 1 corner element.  It's
  ! not called and the corner orientation flag is not set for unstructured meshes


  !
  !
  ! Mark Taylor
  ! pack/unpack full element of data of size (nx,nx)  
  ! user specifies the size when creating the buffer 
  ! input/output arrays are cartesian, and will only unpack 1 corner element 
  ! (even if there are more when running with an unstructured grid) 
  ! This routine is used mostly for testing and to compute the orientation of
  ! an elements corner neighbors
  !
  public :: initGhostBuffer3D      ! init/free buffers used by pack/unpack full and 3D
  public :: FreeGhostBuffer3D
  public :: ghostVpackfull       
  public :: ghostVunpackfull     
  ! same as above, except orientation of element data is preserved
  ! (so boundary data for two adjacent element may not match up)
  public :: ghostVpack_unoriented   
  public :: ghostVunpack_unoriented     


  !
  ! James Overfelt
  ! pack/unpack user specifed halo region "nhc".  
  ! Does not include element edge data (assumes element edge data is C0)
  ! (appropriate for continuous GLL data where the edge data does not need to be sent)
  ! support for unstructed meshes via extra output arrays: sw,se,ne,nw
  ! This routine is currently used by surfaces_mod.F90 to construct the GLL dual grid
  !
  public :: ghostVpack3d          ! pack/unpack specifed halo size (up to 1 element)
                                  ! should be identical to ghostVpack2d except for
                                  ! shape of input array
  public :: ghostVunpack3d        ! returns v including populating halo region of v
                                  ! "extra" corner elements are returned in arrays
                                  ! sw,se,ne,nw
  ! MT TODO: this routine works for unstructed data (where the corner orientation flag is
  ! not set).  So why dont we remove all the "reverse" checks in unpack?


  !
  ! Christoph Erath
  ! pack/unpack partial element of data of size (nx,nx) with user specifed halo size nh
  ! user specifies the sizes when creating the buffer 
  ! buffer has 1 extra dimension (as compared to subroutines above) for multiple tracers
  ! input/output arrays are cartesian, and thus assume at most 1 element at each corner
  ! hence currently only supports cube-sphere grids.
  !
  ! TODO: GhostBufferTR should be removed - we only need GhostBuffer3D, if we can fix
  ! ghostVpack2d below to pass vlyr*ntrac_d instead of two seperate arguments
  !
  public :: initGhostBufferTR     ! ghostbufferTR_t
  public :: FreeGhostBufferTR     ! ghostbufferTR_t

  ! routines which including element edge data  
  ! (used for FVM arrays where edge data is not shared by neighboring elements)
  ! these routines pack/unpack element data with user specified halo size
  !
  ! THESE ROUTINES SHOULD BE MERGED 
  !
  public :: ghostVpack            ! input/output: 
  public :: ghostVunpack          ! v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac_d,timelevels)
                                  
  public :: ghostVpackR           ! used to pack/unpack SPELT "Rp".  What's this?
  public :: ghostVunpackR         ! v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac_d)

  ! routines which do NOT include element edge data
  ! (used for SPELT arrays and GLL point arrays, where edge data is shared and does not need
  ! to be sent/received.
  ! these routines pack/unpack element data with user specifed halo size
  !
  ! THESE ROUTINES CAN ALL BE REPLACED BY ghostVpack3D (if we make extra corner data arrays
  ! an optional argument).  Or at least these should be merged to 1 routine
  public :: ghostVpack2d          ! input/output:  
  public :: ghostVunpack2d        ! v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr, ntrac_d,timelevels)
                                   
                                  ! used to pack/unpack SPELT%sga.  what's this?
  public :: ghostVpack2d_single   ! input/output
  public :: ghostVunpack2d_single !   v(1-nhc:npoints+nhc,1-nhc:npoints+nhc)
                                  
                                  ! used to pack/unpack FV vertex data (velocity/grid)
  public :: ghostVpack2d_level    ! input/output
  public :: ghostVunpack2d_level  ! v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr) 

  ! PFC: Hack to set global option values
  !--------------------------------------
  public :: init_edge_mod
  private:: max_neigh_edges
  private:: max_corner_elem
  private:: ne
  private:: nelemd
  integer:: max_neigh_edges
  integer:: max_corner_elem
  integer:: ne
  integer:: nelemd

  ! Type Definitions
  !--------------------
  type rotation_t
    integer:: nbr               ! nbr direction: north south east west
    integer:: reverse           ! 0=do not reverse order : 1=reverse order
    real(kind=real_kind),pointer:: R(:,:,:) => null()  ! rotation matrix
  end type rotation_t

  type EdgeDescriptor_t
    integer(kind=int_kind)        :: use_rotation
    integer(kind=int_kind)        :: padding
    integer(kind=int_kind),pointer:: putmapP(:)       => null()
    integer(kind=int_kind),pointer:: getmapP(:)       => null()
    integer(kind=int_kind),pointer:: putmapP_ghost(:) => null()
    integer(kind=int_kind),pointer:: getmapP_ghost(:) => null()
    integer(kind=int_kind),pointer:: globalID(:)      => null()
    integer(kind=int_kind),pointer:: loc2buf(:)       => null()
    integer                       :: actual_neigh_edges
    logical(kind=log_kind),pointer:: reverse(:)       => null()
    ! Identifies list of edges that must be rotated, and how
    type(rotation_t)      ,pointer:: rot(:)           => null()
  end type EdgeDescriptor_t

  type EdgeBuffer_t
    real(kind=real_kind),pointer:: buf    (:,:) => null()
    real(kind=real_kind),pointer:: receive(:,:) => null()
    integer                     :: nlyr     ! Number of layers
    integer                     :: nbuf     ! size of the horizontal dimension of the buffers.
  end type EdgeBuffer_t

  type LongEdgeBuffer_t
    integer                       :: nlyr
    integer                       :: nbuf
    integer(kind=int_kind),pointer:: buf    (:,:) => null()
    integer(kind=int_kind),pointer:: receive(:,:) => null()
  end type LongEdgeBuffer_t

  type GhostBufferTR_t
    real(kind=real_kind),pointer:: buf    (:,:,:,:,:) => null()
    real(kind=real_kind),pointer:: receive(:,:,:,:,:) => null()
    integer                     :: nlyr ! Number of layers
    integer                     :: nbuf ! size of the horizontal dimension of the buffers.
  end type GhostBufferTR_t
  
  type GhostBuffer3D_t
    real(kind=real_kind),pointer:: buf    (:,:,:,:) => null()
    real(kind=real_kind),pointer:: receive(:,:,:,:) => null()
    integer                     :: nlyr ! Number of layers
    integer                     :: nhc  ! Number of layers of ghost cells
    integer                     :: np   ! Number of points in a cell
    integer                     :: nbuf ! size of the horizontal dimension of the buffers.
    integer                     :: elem_size ! size of 2D array (first two dimensions of buf())
  end type GhostBuffer3D_t

  ! Wrap pointer so we can make an array of them.
  !----------------------------------------------
  type wrap_ptr
    real(kind=real_kind),pointer:: ptr(:,:) => null()
  end type wrap_ptr

  ! Global Data
  !---------------------
  logical       :: threadsafe=.true.
  type(wrap_ptr):: edgebuff_ptrs(0:1)


contains
  !==================================================================
  subroutine init_edge_mod(I_SEopt)
    ! PFC: Hack to set global option values
    !--------------------------------------
    use SE_Options,only: SEoptions_t
    ! Passed Variables
    !------------------
    type(SEoptions_t),intent(in):: I_SEopt
  
    max_neigh_edges = I_SEopt%max_neigh_edges
    max_corner_elem = I_SEopt%max_corner_elem
    ne              = I_SEopt%ne     
    nelemd          = I_SEopt%nelemd

    ! End Routine
    !-------------
    return
  end subroutine init_edge_mod
  !==================================================================


  !==================================================================
  subroutine initEdgeBuffer(edge,nlyr, buf_ptr,receive_ptr)
    ! initEdgeBuffer:
    !
    ! create an Real based communication buffer
    !
    ! Notes about the buf_ptr/receive_ptr options:
    !
    ! You can pass in 1D pointers to this function. If they are not
    ! associated, they will be allocated and used as buffer space. If they
    ! are associated, their targets will be used as buffer space.
    !
    ! The pointers must not be thread-private.
    !
    ! If an EdgeBuffer_t object is initialized from pre-existing storage
    ! (i.e. buf_ptr is provided and not null), it must *not* be freed,
    ! and must not be used if the underlying storage has been deallocated.
    !
    ! All these restrictions also applied to the old newbuf and newreceive
    ! options.

    ! Workaround for NAG bug.
    ! NAG 5.3.1 dies if you use pointer bounds remapping to set
    ! a pointer that is also a component. So remap to temporary,
    ! then use that to set component pointer.
    !=====================================================
    use SE_Constants,only: np
    !
    ! Passed Variables
    !-------------------
    type(EdgeBuffer_t),target,intent(out):: edge
    integer                  ,intent(in ):: nlyr
    real(kind=real_kind),optional,pointer:: buf_ptr(:)
    real(kind=real_kind),optional,pointer:: receive_ptr(:)
    !
    ! Local Values
    !------------
    real(kind=real_kind),pointer:: tmp_ptr(:,:)
    integer                     :: nbuf,ith

    nbuf      = 4*(np+max_corner_elem)*nelemd
    edge%nlyr = nlyr
    edge%nbuf = nbuf

    ! tracer code might call initedgebuffer() with zero tracers
    !------------------------------------------------------------
    if(nlyr==0) return

!$OMP BARRIER
!   only master thread should allocate the buffer
#if (! defined ELEMENT_OPENMP)
!$OMP MASTER
#endif

    if(present(buf_ptr)) then
      ! If buffer is passed in but not allocated, allocate it.
      !-------------------------------------------------------
      if(.not. associated(buf_ptr)) allocate(buf_ptr(nlyr*nbuf))
      ! Verify dimensions
      !------------------
      if(size(buf_ptr) < nlyr*nbuf) then
        print *,'size(buf_ptr),nlyr,nbuf=',size(buf_ptr),nlyr,nbuf
        call endrun('Error: user provided edge buffer is too small')
      endif
#ifdef HAVE_F2003_PTR_BND_REMAP
      tmp_ptr(1:nlyr,1:nbuf) => buf_ptr
      edge%buf               => tmp_ptr
#else
      ! call F77 routine which will reshape array.
      !-------------------------------------------
      call remap_2D_ptr_buf(edge,nlyr,nbuf,buf_ptr)
#endif
    else
      allocate(edge%buf    (nlyr,nbuf))
    endif

    if(present(receive_ptr)) then
      ! If buffer is passed in but not allocated, allocate it.
      !-------------------------------------------------------
      if(.not. associated(receive_ptr)) allocate(receive_ptr(nlyr*nbuf))
      ! Verify dimensions
      !----------------------
      if(size(receive_ptr) < nlyr*nbuf) then
        print *,'size(receive_ptr),nlyr,nbuf=',size(receive_ptr),nlyr,nbuf
        call endrun('Error: user provided edge buffer is too small')
      endif
#ifdef HAVE_F2003_PTR_BND_REMAP
      tmp_ptr(1:nlyr,1:nbuf) => receive_ptr
      edge%receive           => tmp_ptr
#else
      ! call F77 routine which will reshape array.
      !---------------------------------------------
      call remap_2D_ptr_receive(edge,nlyr,nbuf,receive_ptr)
#endif
    else
      allocate(edge%receive(nlyr,nbuf))
    endif

    edge%buf    (:,:)=0.0D0
    edge%receive(:,:)=0.0D0
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
#endif
!   make sure all threads wait until buffer is allocated
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif

    ! sanity check on edge.  edge must NOT be thread-prive, but must be shared by all threads
    ! the calling program needs to instantiate 'edge' outside the threaded region.
    ! if 'edge' is thread private, it creates flaky openMP problems that are difficut to debug
    ! so lets try and detect it here:
    !---------------------------------------------------------------------------------------
    if(omp_get_num_threads() > 1) then
      ith=omp_get_thread_num()
      if(ith <= 1 ) then
        edgebuff_ptrs(ith)%ptr => edge%buf
      endif
#if (! defined ELEMENT_OPENMP)
      !$OMP BARRIER
      !$OMP MASTER
#endif
      if(.not. associated(edgebuff_ptrs(0)%ptr, edgebuff_ptrs(1)%ptr)) then
        call endrun('ERROR: edge struct appears to be thread-private.')
      endif
#if (! defined ELEMENT_OPENMP)
      !$OMP END MASTER
#endif
    endif

    ! End Routine 
    !-------------
    return
  end subroutine initEdgeBuffer
  !==================================================================


  !==================================================================
  subroutine initLongEdgeBuffer(edge,nlyr)
    ! initLongEdgeBuffer:
    !
    ! create an Integer based communication buffer
    !==================================================================
    use SE_Constants,only: np
    !
    ! Passed Variables
    !------------------
    type(LongEdgeBuffer_t),intent(out):: edge
    integer               ,intent(in ):: nlyr
    !
    ! Local Values
    !--------------
    integer:: nbuf

    ! sanity check for threading
    !-----------------------------
    if(omp_get_num_threads()>1) then
      call endrun('ERROR: initLongEdgeBuffer must be called before threaded reagion')
    endif

    nbuf      = 4*(np+max_corner_elem)*nelemd
    edge%nlyr = nlyr
    edge%nbuf = nbuf

    allocate(edge%buf(nlyr,nbuf))
    edge%buf(:,:)=0

    allocate(edge%receive(nlyr,nbuf))
    edge%receive(:,:)=0

    ! End Routine 
    !-------------
    return
  end subroutine initLongEdgeBuffer
  !==================================================================


  !==================================================================
  subroutine edgeDGVpack(edge,v,vlyr,kptr,desc)
    ! edgeDGVpack:
    !
    ! Pack edges of v into buf for DG stencil
    !==============================================================
    use SE_Constants,only: np
    !
    ! Passed Variables
    !---------------------
    integer             ,intent(in):: vlyr
    type(EdgeBuffer_t)             :: edge
    real(kind=real_kind),intent(in):: v(np,np,vlyr)
    integer             ,intent(in):: kptr
    type(EdgeDescriptor_t)         :: desc

    ! This code is just a wrapper call the 
    !   normal edgeVpack
    !-----------------------------------------
    call edgeVpack(edge,v,vlyr,kptr,desc)

    ! End Routine 
    !-------------
    return
  end subroutine edgeDGVpack
  !==================================================================


  !==================================================================
  subroutine FreeEdgeBuffer(edge) 
    !  FreeEdgeBuffer:
    !
    !  Freed an edge communication buffer
    !===============================================================
    !
    ! Passed Variables
    !-----------------
    type(EdgeBuffer_t),intent(inout):: edge

#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
    edge%nbuf=0
    edge%nlyr=0
    deallocate(edge%buf    )
    deallocate(edge%receive)
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
#endif

    ! End Routine 
    !-------------
    return
  end subroutine FreeEdgeBuffer
  !==================================================================


  !==================================================================
  subroutine FreeGhostBuffer3D(buffer) 
    ! FreeGhostBuffer3D:
    !
    !===========================================================
    !
    ! Passed Variables
    !--------------------
    type(Ghostbuffer3d_t),intent(inout):: buffer

#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
    buffer%nbuf=0
    buffer%nlyr=0
    deallocate(buffer%buf    )
    deallocate(buffer%receive)
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
#endif

    ! End Routine 
    !-------------
    return
  end subroutine FreeGhostBuffer3D
  !==================================================================


  !==================================================================
  subroutine FreeLongEdgeBuffer(edge) 
    !  FreeLongEdgeBuffer:
    !
    !  Freed an edge communication buffer
    !================================================================
    !
    ! Passed Variables
    !-------------------
    type(LongEdgeBuffer_t),intent(inout):: edge

    edge%nbuf=0
    edge%nlyr=0
    deallocate(edge%buf    )
    deallocate(edge%receive)

    ! End Routine 
    !-------------
    return
  end subroutine FreeLongEdgeBuffer
  !==================================================================


  !==================================================================
  subroutine edgeVpack(edge,v,vlyr,kptr,desc)
    ! edgeVpack:
    !
    !> @brief Pack edges of v into an edge buffer for boundary exchange.
    !
    !> This subroutine packs for one or more vertical layers into an edge 
    !! buffer. If the buffer associated with edge is not large enough to 
    !! hold all vertical layers you intent to pack, the method will 
    !! halt the program with a call to parallel_mod::endrun().
    !! @param[in] edge Edge Buffer into which the data will be packed.
    !! This buffer must be previously allocated with initEdgeBuffer().
    !! @param[in] v The data to be packed.
    !! @param[in] vlyr Number of vertical level coming into the subroutine
    !! for packing for input v.
    !! @param[in] kptr Vertical pointer to the place in the edge buffer where 
    !! data will be located.
    !=======================================================================
    use SE_Constants,only: np
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !-------------------
    integer               ,intent(in):: vlyr
    type(EdgeBuffer_t)               :: edge
    real(kind=real_kind)  ,intent(in):: v(np,np,vlyr)
    integer               ,intent(in):: kptr
    type(EdgeDescriptor_t),intent(in):: desc
    !
    ! Local Values
    !-------------
    logical,parameter:: UseUnroll = .TRUE.
    integer:: is,ie,in,iw
    integer:: ii,kk,ir,ll

!PFC    call t_startf('edge_pack')
    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
      threadsafe=.true.
    endif

    is = desc%putmapP(south)
    ie = desc%putmapP(east )
    in = desc%putmapP(north)
    iw = desc%putmapP(west )
    if(edge%nlyr < (kptr+vlyr) ) then
      call endrun('edgeVpack: Buffer overflow: size of the vertical dimension must be increased!')
    endif

    if((MODULO(np,2) == 0).and.(UseUnroll)) then 
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i)
#endif
      do kk=1,vlyr
      do ii=1,np,2
        edge%buf(kptr+kk,is+ii  ) = v(ii  ,   1,kk)
        edge%buf(kptr+kk,is+ii+1) = v(ii+1,   1,kk)
        edge%buf(kptr+kk,ie+ii  ) = v(np  ,ii  ,kk)
        edge%buf(kptr+kk,ie+ii+1) = v(np  ,ii+1,kk)
        edge%buf(kptr+kk,in+ii  ) = v(ii  ,np  ,kk)
        edge%buf(kptr+kk,in+ii+1) = v(ii+1,np  ,kk)
        edge%buf(kptr+kk,iw+ii  ) = v(   1,ii  ,kk)
        edge%buf(kptr+kk,iw+ii+1) = v(   1,ii+1,kk)
      end do
      end do
    else
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(kk,ii)
#endif
      do kk=1,vlyr
      do ii=1,np
        edge%buf(kptr+kk,is+ii) = v(ii, 1,kk)
        edge%buf(kptr+kk,ie+ii) = v(np,ii,kk)
        edge%buf(kptr+kk,in+ii) = v(ii,np,kk)
        edge%buf(kptr+kk,iw+ii) = v( 1,ii,kk)
      end do
      end do
    endif

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !------------------------------------------------------------------------
    if(desc%reverse(south)) then
      is = desc%putmapP(south)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(kk,ii,ir)
#endif
      do kk=1,vlyr
      do ii=1,np
        ir = np-ii+1
        edge%buf(kptr+kk,is+ir)=v(ii,1,kk)
      end do
      end do
    endif

    if(desc%reverse(east)) then
      ie = desc%putmapP(east)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(kk,ii,ir)
#endif
      do kk=1,vlyr
      do ii=1,np
        ir = np-ii+1
        edge%buf(kptr+kk,ie+ir)=v(np,ii,kk)
      end do
      end do
    endif

    if(desc%reverse(north)) then
      in = desc%putmapP(north)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(kk,ki,ir)
#endif
      do kk=1,vlyr
      do ii=1,np
        ir = np-ii+1
        edge%buf(kptr+kk,in+ir)=v(ii,np,kk)
      end do
      end do
    endif

    if(desc%reverse(west)) then
      iw = desc%putmapP(west)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(kk,ii,ir)
#endif
      do kk=1,vlyr
      do ii=1,np
        ir = np-ii+1
        edge%buf(kptr+kk,iw+ir)=v(1,ii,kk)
      end do
      end do
    endif

    ! SWEST
    !-------
    do ll=swest,swest+max_corner_elem-1
      if(desc%putmapP(ll) /= -1) then
        do kk=1,vlyr
          edge%buf(kptr+kk,desc%putmapP(ll)+1)=v(1  ,1 ,kk)
        end do
      endif
    end do

    ! SEAST
    !-------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(desc%putmapP(ll) /= -1) then
        do kk=1,vlyr
          edge%buf(kptr+kk,desc%putmapP(ll)+1)=v(np ,1 ,kk)
        end do
      endif
    end do

    ! NEAST
    !-------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(desc%putmapP(ll) /= -1) then
        do kk=1,vlyr
          edge%buf(kptr+kk,desc%putmapP(ll)+1)=v(np ,np,kk)
        end do
      endif
    end do

    ! NWEST
    !-------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(desc%putmapP(ll) /= -1) then
        do kk=1,vlyr
          edge%buf(kptr+kk,desc%putmapP(ll)+1)=v(1  ,np,kk)
        end do
      endif
    end do

!PFC    call t_stopf('edge_pack')

    ! End Routine 
    !-------------
    return
  end subroutine edgeVpack
  !==================================================================


  !==================================================================
  subroutine LongEdgeVpack(edge,v,vlyr,kptr,desc)
    ! LongEdgeVpack:
    !
    ! Pack edges of v into buf...
    !====================================================================
    use SE_Constants,only: np
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !------------------
    integer               ,intent(in):: vlyr
    type(LongEdgeBuffer_t)           :: edge
    integer(kind=int_kind),intent(in):: v(np,np,vlyr)
    integer               ,intent(in):: kptr
    type(EdgeDescriptor_t),intent(in):: desc
    !
    ! Local Values
    !--------------
    logical,parameter:: UseUnroll = .TRUE.
    integer:: is,ie,in,iw
    integer:: ii,kk,ir,ll

    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
      threadsafe=.true.
    endif

    is = desc%putmapP(south)
    ie = desc%putmapP(east )
    in = desc%putmapP(north)
    iw = desc%putmapP(west )

    if((MODULO(np,2) == 0).and.(UseUnroll)) then 
      do kk=1,vlyr
      do ii=1,np,2
        edge%buf(kptr+kk,is+ii  ) = v(ii  ,   1,kk)
        edge%buf(kptr+kk,is+ii+1) = v(ii+1,   1,kk)
        edge%buf(kptr+kk,ie+ii  ) = v(np  ,ii  ,kk)
        edge%buf(kptr+kk,ie+ii+1) = v(np  ,ii+1,kk)
        edge%buf(kptr+kk,in+ii  ) = v(ii  ,np  ,kk)
        edge%buf(kptr+kk,in+ii+1) = v(ii+1,np  ,kk)
        edge%buf(kptr+kk,iw+ii  ) = v(   1,ii  ,kk)
        edge%buf(kptr+kk,iw+ii+1) = v(   1,ii+1,kk)
      end do
      end do
    else
      do kk=1,vlyr
      do ii=1,np
        edge%buf(kptr+kk,is+ii) = v(ii,1 ,kk)
        edge%buf(kptr+kk,ie+ii) = v(np,ii,kk)
        edge%buf(kptr+kk,in+ii) = v(ii,np,kk)
        edge%buf(kptr+kk,iw+ii) = v(1 ,ii,kk)
      end do
      end do
    endif

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !------------------------------------------------------------------------
    if(desc%reverse(south)) then
      is = desc%putmapP(south)
      do kk=1,vlyr
      do ii=1,np
        ir = np-ii+1
        edge%buf(kptr+kk,is+ir)=v(ii,1,kk)
      end do
      end do
    endif

    if(desc%reverse(east)) then
      ie = desc%putmapP(east)
      do kk=1,vlyr
      do ii=1,np
        ir = np-ii+1
        edge%buf(kptr+kk,ie+ir)=v(np,ii,kk)
      end do
      end do
    endif

    if(desc%reverse(north)) then
      in = desc%putmapP(north)
      do kk=1,vlyr
      do ii=1,np
        ir = np-ii+1
        edge%buf(kptr+kk,in+ir)=v(ii,np,kk)
      end do
      end do
    endif

    if(desc%reverse(west)) then
      iw = desc%putmapP(west)
      do kk=1,vlyr
      do ii=1,np
        ir = np-ii+1
        edge%buf(kptr+kk,iw+ir)=v(1,ii,kk)
      end do
      end do
    endif

    ! SWEST
    !-------
    do ll=swest,swest+max_corner_elem-1
      if(desc%putmapP(ll) /= -1) then
        do kk=1,vlyr
          edge%buf(kptr+kk,desc%putmapP(ll)+1)=v(1  ,1 ,kk)
        end do
      endif
    end do

    ! SEAST
    !-------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(desc%putmapP(ll) /= -1) then
        do kk=1,vlyr
          edge%buf(kptr+kk,desc%putmapP(ll)+1)=v(np ,1 ,kk)
        end do
      endif
    end do

    ! NEAST
    !-------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(desc%putmapP(ll) /= -1) then
        do kk=1,vlyr
          edge%buf(kptr+kk,desc%putmapP(ll)+1)=v(np ,np,kk)
        end do
      endif
    end do

    ! NWEST
    !-------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(desc%putmapP(ll) /= -1) then
        do kk=1,vlyr
          edge%buf(kptr+kk,desc%putmapP(ll)+1)=v(1  ,np,kk)
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine LongEdgeVpack
  !==================================================================


  !==================================================================
  subroutine edgeVunpack(edge,v,vlyr,kptr,desc)
    ! edgeVunpack:
    !
    ! Unpack edges from edge buffer into v...
    !===============================================================
    use SE_Constants,only: np
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !----------------------
    integer             ,intent(in   ):: vlyr
    type(EdgeBuffer_t)  ,intent(in   ):: edge
    real(kind=real_kind),intent(inout):: v(np,np,vlyr)
    integer             ,intent(in   ):: kptr
    type(EdgeDescriptor_t)            :: desc
    !
    ! Local Values
    !--------------
    logical,parameter:: UseUnroll = .TRUE.
    integer:: is,ie,in,iw
    integer:: ii,kk,ll

!PFC    call t_startf('edge_unpack')
    threadsafe=.false.

    is=desc%getmapP(south)
    ie=desc%getmapP(east )
    in=desc%getmapP(north)
    iw=desc%getmapP(west )

    if((MODULO(np,2) == 0).and.(UseUnroll)) then 
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i)
#endif
      do kk=1,vlyr
      do ii=1,np,2
        v(ii  ,   1,kk) = v(ii  ,   1,kk)+edge%buf(kptr+kk,is+ii  )
        v(ii+1,   1,kk) = v(ii+1,   1,kk)+edge%buf(kptr+kk,is+ii+1)
        v(np  ,ii  ,kk) = v(np  ,ii  ,kk)+edge%buf(kptr+kk,ie+ii  )
        v(np  ,ii+1,kk) = v(np  ,ii+1,kk)+edge%buf(kptr+kk,ie+ii+1)
        v(ii  ,np  ,kk) = v(ii  ,np  ,kk)+edge%buf(kptr+kk,in+ii  )
        v(ii+1,np  ,kk) = v(ii+1,np  ,kk)+edge%buf(kptr+kk,in+ii+1)
        v(   1,ii  ,kk) = v(   1,ii  ,kk)+edge%buf(kptr+kk,iw+ii  )
        v(   1,ii+1,kk) = v(   1,ii+1,kk)+edge%buf(kptr+kk,iw+ii+1)
      end do
      end do
    else
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i)
#endif
      do kk=1,vlyr
      do ii=1,np
        v(ii ,1  ,kk) = v(ii ,1  ,kk)+edge%buf(kptr+kk,is+ii  )
        v(np ,ii ,kk) = v(np ,ii ,kk)+edge%buf(kptr+kk,ie+ii  )
        v(ii ,np ,kk) = v(ii ,np ,kk)+edge%buf(kptr+kk,in+ii  )
        v(1  ,ii ,kk) = v(1  ,ii ,kk)+edge%buf(kptr+kk,iw+ii  )
      end do
      end do
    endif

    ! SWEST
    !---------
    do ll=swest,swest+max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(1  ,1 ,kk)=v(1 ,1 ,kk)+edge%buf(kptr+kk,desc%getmapP(ll)+1)
        end do
      endif
    end do

    ! SEAST
    !---------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(np ,1 ,kk)=v(np,1 ,kk)+edge%buf(kptr+kk,desc%getmapP(ll)+1)
        end do
      endif
    end do

    ! NEAST
    !---------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(np ,np,kk)=v(np,np,kk)+edge%buf(kptr+kk,desc%getmapP(ll)+1)
        end do
      endif
    end do

    ! NWEST
    !---------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(1  ,np,kk)=v(1 ,np,kk)+edge%buf(kptr+kk,desc%getmapP(ll)+1)
        end do
      endif
    end do

!PFC    call t_stopf('edge_unpack')

    ! End Routine 
    !-------------
    return
  end subroutine edgeVunpack
  !==================================================================


  !==================================================================
  subroutine edgeVunpackVert(edge,v,desc)
    ! edgeVunpackVert:
    !
    !======================================================================
    use SE_Constants,only: np
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    use coordinate_systems_mod,only: cartesian3D_t
    !
    ! Passed Variables
    !--------------------
    type(EdgeBuffer_t) ,intent(inout):: edge
    type(cartesian3D_t),intent(  out):: v(:,:,:)
    type(EdgeDescriptor_t)           :: desc
    !
    ! Local Values
    !---------------
    logical,parameter:: UseUnroll = .TRUE.
    integer :: is,ie,in,iw
    integer :: ii,kk,ll

    threadsafe=.false.

    if((max_corner_elem.ne.1).and.(ne==0)) then
      ! MNL: this is used to construct the dual grid on the cube,
      !      currently only supported for the uniform grid. If
      !      this is desired on a refined grid, a little bit of
      !      work will be required.
      call endrun("edgeVunpackVert should not be called with unstructured meshes")
    endif

    is=desc%getmapP(south)
    ie=desc%getmapP(east )
    in=desc%getmapP(north)
    iw=desc%getmapP(west )

    ! N+S
    !-----------
    do ii=1,np/2
      ! North
      !----------
      v(3,ii,np)%x = edge%buf(1,in+ii) 
      v(3,ii,np)%y = edge%buf(2,in+ii) 
      v(3,ii,np)%z = edge%buf(3,in+ii) 
      ! South
      !----------
      v(2,ii,1)%x  = edge%buf(1,is+ii) 
      v(2,ii,1)%y  = edge%buf(2,is+ii) 
      v(2,ii,1)%z  = edge%buf(3,is+ii) 
    end do

    do ii=np/2+1,np
      ! North
      !-----------
      v(4,ii,np)%x = edge%buf(1,in+ii) 
      v(4,ii,np)%y = edge%buf(2,in+ii) 
      v(4,ii,np)%z = edge%buf(3,in+ii) 
      ! South
      !-----------
      v(1,ii,1)%x  = edge%buf(1,is+ii) 
      v(1,ii,1)%y  = edge%buf(2,is+ii) 
      v(1,ii,1)%z  = edge%buf(3,is+ii)        
    end do

    do ii=1,np/2
      ! East
      !----------
      v(3,np,ii)%x = edge%buf(1,ie+ii)
      v(3,np,ii)%y = edge%buf(2,ie+ii)
      v(3,np,ii)%z = edge%buf(3,ie+ii)       
      ! West
      !----------
      v(4,1,ii)%x  = edge%buf(1,iw+ii)
      v(4,1,ii)%y  = edge%buf(2,iw+ii)
      v(4,1,ii)%z  = edge%buf(3,iw+ii)
    end do

    do ii=np/2+1,np
      ! East
      !-------
      v(2,np,ii)%x = edge%buf(1,ie+ii)
      v(2,np,ii)%y = edge%buf(2,ie+ii)
      v(2,np,ii)%z = edge%buf(3,ie+ii)       
      ! West
      !-------
      v(1,1,ii)%x  = edge%buf(1,iw+ii)
      v(1,1,ii)%y  = edge%buf(2,iw+ii)
      v(1,1,ii)%z  = edge%buf(3,iw+ii)
    end do

    ! SWEST
    !--------
    do ll=swest,swest+max_corner_elem-1
      ! find the one active corner, then exist
      !---------------------------------------
      if(desc%getmapP(ll) /= -1) then 
        v(1,1,1)%x=edge%buf(1,desc%getmapP(ll)+1)
        v(1,1,1)%y=edge%buf(2,desc%getmapP(ll)+1)
        v(1,1,1)%z=edge%buf(3,desc%getmapP(ll)+1)
        exit 
      else
        v(1,1,1)%x=0_real_kind
        v(1,1,1)%y=0_real_kind
        v(1,1,1)%z=0_real_kind
      endif
    end do

    ! SEAST
    !---------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      ! find the one active corner, then exist
      if(desc%getmapP(ll) /= -1) then 
        v(2,np,1)%x=edge%buf(1,desc%getmapP(ll)+1)
        v(2,np,1)%y=edge%buf(2,desc%getmapP(ll)+1)
        v(2,np,1)%z=edge%buf(3,desc%getmapP(ll)+1)
        exit
      else
        v(2,np,1)%x=0_real_kind
        v(2,np,1)%y=0_real_kind
        v(2,np,1)%z=0_real_kind
      endif
    end do

    ! NEAST
    !----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      ! find the one active corner, then exist
      if(desc%getmapP(ll) /= -1) then 
        v(3,np,np)%x=edge%buf(1,desc%getmapP(ll)+1)
        v(3,np,np)%y=edge%buf(2,desc%getmapP(ll)+1)
        v(3,np,np)%z=edge%buf(3,desc%getmapP(ll)+1)
        exit
      else
        v(3,np,np)%x=0_real_kind
        v(3,np,np)%y=0_real_kind
        v(3,np,np)%z=0_real_kind
      endif
    end do

    ! NWEST
    !---------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      ! find the one active corner, then exist
      if(desc%getmapP(ll) /= -1) then 
        v(4,1,np)%x=edge%buf(1,desc%getmapP(ll)+1)
        v(4,1,np)%y=edge%buf(2,desc%getmapP(ll)+1)
        v(4,1,np)%z=edge%buf(3,desc%getmapP(ll)+1)
        exit
      else
        v(4,1,np)%x=0_real_kind
        v(4,1,np)%y=0_real_kind
        v(4,1,np)%z=0_real_kind
      endif
    end do

    ! Fill the missing vertex info
    !-------------------------------
    do ii=2,np/2
      ! North
      !---------
      v(4,ii,np)%x = v(3,ii-1,np)%x 
      v(4,ii,np)%y = v(3,ii-1,np)%y
      v(4,ii,np)%z = v(3,ii-1,np)%z
      ! South
      !---------
      v(1,ii,1)%x  = v(2,ii-1,1)%x 
      v(1,ii,1)%y  = v(2,ii-1,1)%y 
      v(1,ii,1)%z  = v(2,ii-1,1)%z 
    end do

    do ii=np/2+1,np-1
      ! North
      !---------
      v(3,ii,np)%x = v(4,ii+1,np)%x 
      v(3,ii,np)%y = v(4,ii+1,np)%y
      v(3,ii,np)%z = v(4,ii+1,np)%z
      ! South
      !---------
      v(2,ii,1)%x  = v(1,ii+1,1)%x 
      v(2,ii,1)%y  = v(1,ii+1,1)%y
      v(2,ii,1)%z  = v(1,ii+1,1)%z
    end do

    do ii=2,np/2
      ! East
      !---------
      v(2,np,ii)%x = v(3,np,ii-1)%x
      v(2,np,ii)%y = v(3,np,ii-1)%y
      v(2,np,ii)%z = v(3,np,ii-1)%z
      ! West
      !---------
      v(1,1,ii)%x  = v(4,1,ii-1)%x
      v(1,1,ii)%y  = v(4,1,ii-1)%y
      v(1,1,ii)%z  = v(4,1,ii-1)%z
    end do

    do ii=np/2+1,np-1
      ! East
      !---------
      v(3,np,ii)%x = v(2,np,ii+1)%x 
      v(3,np,ii)%y = v(2,np,ii+1)%y
      v(3,np,ii)%z = v(2,np,ii+1)%z
      ! West
      !---------
      v(4,1,ii)%x  = v(1,1,ii+1)%x 
      v(4,1,ii)%y  = v(1,1,ii+1)%y
      v(4,1,ii)%z  = v(1,1,ii+1)%z
    end do

    ! End Routine 
    !-------------
    return
  end subroutine edgeVunpackVert
  !==================================================================


  !==================================================================
  subroutine edgeDGVunpack(edge,v,vlyr,kptr,desc)
    ! edgeDGVunpack:
    !
    ! Unpack edges from edge buffer into v...
    !============================================================
    use SE_Constants,only: np
    use SE_Options  ,only: north, south, east, west
    !
    ! Passed Variables
    !-------------------
    integer             ,intent(in ):: vlyr
    type(EdgeBuffer_t)  ,intent(in ):: edge
    real(kind=real_kind),intent(out):: v(0:np+1,0:np+1,vlyr)
    integer             ,intent(in ):: kptr
    type(EdgeDescriptor_t)          :: desc
    !
    ! Local Values
    !------------------
    integer :: is,ie,in,iw
    integer :: ii,kk

    threadsafe=.false.

    is=desc%getmapP(south)
    ie=desc%getmapP(east )
    in=desc%getmapP(north)
    iw=desc%getmapP(west )

    do kk=1,vlyr
    do ii=1,np
      v(ii  ,   0,kk)=edge%buf(kptr+kk,is+ii)
      v(np+1,ii  ,kk)=edge%buf(kptr+kk,ie+ii)
      v(ii  ,np+1,kk)=edge%buf(kptr+kk,in+ii)
      v(   0,ii  ,kk)=edge%buf(kptr+kk,iw+ii)
    end do
    end do

    ! End Routine 
    !-------------
    return
  end subroutine edgeDGVunpack
  !==================================================================


  !==================================================================
  subroutine edgeVunpackMAX(edge,v,vlyr,kptr,desc)
    ! edgeVunpackMIN/MAX:
    !
    ! Finds the Min/Max edges from edge buffer into v...
    !================================================================
    use SE_Constants,only: np
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !-------------------
    integer               ,intent(in   ):: vlyr
    type(EdgeBuffer_t)    ,intent(in   ):: edge
    real(kind=real_kind)  ,intent(inout):: v(np,np,vlyr)
    integer               ,intent(in   ):: kptr
    type(EdgeDescriptor_t),intent(in   ):: desc
    !
    ! Local Values
    !----------------
    integer :: is,ie,in,iw
    integer :: ii,kk,ll

    threadsafe=.false.

    is=desc%getmapP(south)
    ie=desc%getmapP(east )
    in=desc%getmapP(north)
    iw=desc%getmapP(west )

    do kk=1,vlyr
    do ii=1,np
      v(ii, 1,kk) = MAX(v(ii, 1,kk),edge%buf(kptr+kk,is+ii))
      v(np,ii,kk) = MAX(v(np,ii,kk),edge%buf(kptr+kk,ie+ii))
      v(ii,np,kk) = MAX(v(ii,np,kk),edge%buf(kptr+kk,in+ii))
      v( 1,ii,kk) = MAX(v( 1,ii,kk),edge%buf(kptr+kk,iw+ii))
    end do
    end do

    ! SWEST
    !----------
    do ll=swest,swest+max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(1,1,kk)=MAX(v(1,1,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! SEAST
    !----------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(np,1,kk)=MAX(v(np,1,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! NEAST
    !----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then
        do kk=1,vlyr
          v(np,np,kk)=MAX(v(np,np,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! NWEST
    !----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(1,np,kk)=MAX(v(1 ,np,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine edgeVunpackMAX
  !==================================================================


  !==================================================================
  subroutine edgeVunpackMIN(edge,v,vlyr,kptr,desc)
    ! edgeVunpackMIN:
    !
    !=================================================================
    use SE_Constants,only: np
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !-------------------
    integer               ,intent(in   ):: vlyr
    type (EdgeBuffer_t)   ,intent(in   ):: edge
    real (kind=real_kind) ,intent(inout):: v(np,np,vlyr)
    integer               ,intent(in   ):: kptr
    type(EdgeDescriptor_t),intent(in   ):: desc
    !
    ! Local Values
    !----------------
    integer :: is,ie,in,iw
    integer :: ii,kk,ll

    threadsafe=.false.

    is=desc%getmapP(south)
    ie=desc%getmapP(east )
    in=desc%getmapP(north)
    iw=desc%getmapP(west )

    do kk=1,vlyr
    do ii=1,np
      v(ii, 1,kk) = MIN(v(ii, 1,kk),edge%buf(kptr+kk,is+ii))
      v(np,ii,kk) = MIN(v(np,ii,kk),edge%buf(kptr+kk,ie+ii))
      v(ii,np,kk) = MIN(v(ii,np,kk),edge%buf(kptr+kk,in+ii))
      v( 1,ii,kk) = MIN(v( 1,ii,kk),edge%buf(kptr+kk,iw+ii))
    end do
    end do

    ! SWEST
    !-----------
    do ll=swest,swest+max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(1,1,kk)=MIN(v(1,1,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! SEAST
    !-----------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(np,1,kk)=MIN(v(np,1,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! NEAST
    !-----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(np,np,kk)=MIN(v(np,np,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! NWEST
    !-----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(1,np,kk)=MIN(v(1,np,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine edgeVunpackMIN
  !==================================================================


  !==================================================================
  subroutine LongEdgeVunpackMIN(edge,v,vlyr,kptr,desc)
    ! LongEdgeVunpackMIN:
    !
    ! Finds the Min edges from edge buffer into v...
    !==================================================================
    use SE_Constants,only: np
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !-------------------
    integer               ,intent(in   ):: vlyr
    type(LongEdgeBuffer_t),intent(in   ):: edge
    integer(kind=int_kind),intent(inout):: v(np,np,vlyr)
    integer               ,intent(in   ):: kptr
    type(EdgeDescriptor_t),intent(in   ):: desc
    !
    ! Local Values
    !----------------
    integer :: is,ie,in,iw
    integer :: ii,kk,ll

    threadsafe=.false.

    is=desc%getmapP(south)
    ie=desc%getmapP(east )
    in=desc%getmapP(north)
    iw=desc%getmapP(west )

    do kk=1,vlyr
    do ii=1,np
      v(ii, 1,kk) = MIN(v(ii, 1,kk),edge%buf(kptr+kk,is+ii))
      v(np,ii,kk) = MIN(v(np,ii,kk),edge%buf(kptr+kk,ie+ii))
      v(ii,np,kk) = MIN(v(ii,np,kk),edge%buf(kptr+kk,in+ii))
      v( 1,ii,kk) = MIN(v( 1,ii,kk),edge%buf(kptr+kk,iw+ii))
    end do
    end do

    ! SWEST
    !---------
    do ll=swest,swest+max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(1,1,kk)=MIN(v(1,1,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! SEAST
    !---------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(np,1,kk)=MIN(v(np,1,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! NEAST
    !---------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(np,np,kk)=MIN(v(np,np,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! NWEST
    !---------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(desc%getmapP(ll) /= -1) then 
        do kk=1,vlyr
          v(1,np,kk)=MIN(v(1,np,kk),edge%buf(kptr+kk,desc%getmapP(ll)+1))
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine LongEdgeVunpackMIN
  !==================================================================


  !==================================================================
  subroutine edgerotate(edge,vlyr,kptr,desc)
    ! edgerotate:
    !
    ! Rotate edges in buffer...
    !============================================================
    use SE_Constants,only: np
    !
    ! Passed Variables
    !----------------------
    type(EdgeBuffer_t)    :: edge      ! edge struct
    integer    ,intent(in):: vlyr      ! number of 2d vector fields to rotate
    integer    ,intent(in):: kptr      ! layer pointer into edge buffer
    type(EdgeDescriptor_t):: desc
    !
    ! Local Values
    !-------------------
    real(kind=real_kind),pointer:: R(:,:,:)
    real(kind=real_kind)        :: tmp1,tmp2
    integer                     :: irot,ia,nbr
    integer                     :: ii,kk,k1,k2

#ifdef _USEASSOCIATED
    if(associated(rot)) then
#else
    if(desc%use_rotation == 1) then
#endif
      do irot=1,SIZE(desc%rot)
        nbr =  desc%rot(irot)%nbr
        R   => desc%rot(irot)%R
        ia  =  desc%putmapP(nbr)

        ! If nbr direction is (1-4) => is an edge
        !----------------------------------------
        if(nbr <= 4) then
          !  Is an edge. Rotate it in place
          !----------------------------------
          do ii=1,np
          do kk=1,vlyr,2
            k1 = kptr + kk
            k2 = k1   + 1
            tmp1=R(1,1,ii)*edge%buf(k1,ia+ii) + R(1,2,ii)*edge%buf(k2,ia+ii)
            tmp2=R(2,1,ii)*edge%buf(k1,ia+ii) + R(2,2,ii)*edge%buf(k2,ia+ii)
            edge%buf(k1,ia+ii)=tmp1
            edge%buf(k2,ia+ii)=tmp2
          end do
          end do
        else
          ! Is an element corner point, but not a cube corner
          ! point, just rotate it in place.
          !----------------------------------------------------
          if(ia /= -1) then
            do kk=1,vlyr,2
              k1 = kptr + kk
              k2 = k1   + 1
              tmp1=R(1,1,1)*edge%buf(k1,ia+1) + R(1,2,1)*edge%buf(k2,ia+1)
              tmp2=R(2,1,1)*edge%buf(k1,ia+1) + R(2,2,1)*edge%buf(k2,ia+1)
              edge%buf(k1,ia+1)=tmp1
              edge%buf(k2,ia+1)=tmp2
            end do
          endif
        endif
      end do
    endif

    ! End Routine 
    !-------------
    return
  end subroutine edgerotate
  !==================================================================


  !==================================================================
   function buffermap(inum,facet) result(loc)
     ! buffermap:
     !
     ! buffermap translates element number, inum and
     ! element edge/corner, facet, into an edge buffer 
     ! memory location, loc.
     !==============================================================
     use SE_Constants,only: np
     !
     ! Passed Variables
     !--------------------
     integer,intent(in):: inum   
     integer,intent(in):: facet
     integer           :: loc

     if(facet>4) then
       if(inum == -1) then
         loc = inum
       else
         loc=(inum-1)*(4*np+4)+4*np+(facet-5)
       endif
     else
       loc=(inum-1)*(4*np+4)+np*(facet-1)
     endif

    ! End Function 
    !-------------
    return
  end function buffermap
  !==================================================================


  !==================================================================
  subroutine FreeGhostBufferTR(ghost) 
    !  FreeGhostBuffer:
    !  Author: Christoph Erath, Mark Taylor
    !  Freed an ghostpoints communication buffer
    !================================================================
    !
    ! Passed Variables
    !----------------------
    type(GhostBuffertr_t),intent(inout):: ghost

#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
    ghost%nbuf=0
    ghost%nlyr=0
    deallocate(ghost%buf    )
    deallocate(ghost%receive)
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
#endif

    ! End Routine 
    !-------------
    return
  end subroutine FreeGhostBufferTR 
  !==================================================================


  !==================================================================
  subroutine GhostVpackfull(edge,v,nc1,nc2,nc,vlyr,kptr,desc)
    ! GhostVpackfull:
    !
    !> @brief Pack edges of v into an edge buffer for boundary exchange.
    !
    !> This subroutine packs for one or more vertical layers into an edge 
    !! buffer. If the buffer associated with edge is not large enough to 
    !! hold all vertical layers you intent to pack, the method will 
    !! halt the program with a call to parallel_mod::endrun().
    !! @param[in] edge Ghost Buffer into which the data will be packed.
    !! This buffer must be previously allocated with initghostbufferfull().
    !! @param[in] v The data to be packed.
    !! @param[in] vlyr Number of vertical level coming into the subroutine
    !! for packing for input v.
    !! @param[in] kptr Vertical pointer to the place in the edge buffer where 
    !! data will be located.
    !=========================================================================
    use SE_Options,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !----------------------
    integer               ,intent(in):: vlyr
    type(Ghostbuffer3D_t)            :: edge
    real(kind=real_kind)  ,intent(in):: v(nc1:nc2,nc1:nc2,vlyr)
    integer               ,intent(in):: nc1,nc2,nc
    integer               ,intent(in):: kptr
    type(EdgeDescriptor_t),intent(in):: desc
    !
    ! Local Values
    !-----------------
    integer :: is,ie,in,iw
    integer :: ii,kk,ir,ll,ee

    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
      threadsafe=.true.
    endif

    ! Example convenction for buffer to the north:
    !  buf(:,,:,i,e) 
    !   each "edge" is a row of data (i=1,np) in the element 
    !     north most row of data goes into e=1
    !     next row of data goes into e=2
    !     .... 
    !     south most row of data goes into e=np
    !   We need to pack this way to preserve the orientation
    !   so the data can be unpacked correctly
    !
    ! note: we think of buf as dimensioned buf(k,is,i,e)
    ! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
    !----------------------------------------------------------
    is = desc%putmapP_ghost(south) 
    ie = desc%putmapP_ghost(east )  
    in = desc%putmapP_ghost(north) 
    iw = desc%putmapP_ghost(west )  

#if 1
    if(is > edge%nbuf) call endrun('error is=')
    if(ie > edge%nbuf) call endrun('error ie=')
    if(in > edge%nbuf) call endrun('error in=')
    if(iw > edge%nbuf) call endrun('error iw=')
    if(is <         1) call endrun('error is=0')
    if(ie <         1) call endrun('error ie=0')
    if(in <         1) call endrun('error in=0')
    if(iw <         1) call endrun('error iw=0')
#endif
!   print *,nc,is,ie,in,iw

    do kk=1,vlyr
    do ee=1,nc
    do ii=1,nc
      edge%buf(ii,ee,kptr+kk,is) = v(ii     ,ee     ,kk)
      edge%buf(ii,ee,kptr+kk,ie) = v(nc-ee+1,ii     ,kk)
      edge%buf(ii,ee,kptr+kk,in) = v(ii     ,nc-ee+1,kk)
      edge%buf(ii,ee,kptr+kk,iw) = v(ee     ,ii     ,kk)
    end do
    end do
    end do

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !  Check if the edge orientation of the recieving element is different
    !  if it is, swap the order of data in the edge
    !---------------------------------------------------------------------
    if(desc%reverse(south)) then
      is = desc%putmapP_ghost(south)
      do ee=1,nc
      do kk=1,vlyr
      do ii=1,nc
        ir = nc-ii+1
        edge%buf(ir,ee,kptr+kk,is)=v(ii,ee,kk)
      end do
      end do
      end do
    endif

    if(desc%reverse(east)) then
      ie = desc%putmapP_ghost(east)
      do ee=1,nc
      do kk=1,vlyr
      do ii=1,nc
        ir = nc-ii+1
        edge%buf(ir,ee,kptr+kk,ie)=v(nc-ee+1,ii,kk)
      end do
      end do
      end do
    endif

    if(desc%reverse(north)) then
      in = desc%putmapP_ghost(north)
      do ee=1,nc
      do kk=1,vlyr
      do ii=1,nc
        ir = nc-ii+1
        edge%buf(ir,ee,kptr+kk,in)=v(ii,nc-ee+1,kk)
      end do
      end do
      end do
    endif

    if(desc%reverse(west)) then
      iw = desc%putmapP_ghost(west)
      do ee=1,nc
      do kk=1,vlyr
      do ii=1,nc
        ir = nc-ii+1
        edge%buf(ir,ee,kptr+kk,iw)=v(ee,ii,kk)
      end do
      end do
      end do
    endif

    ! corners.  this is difficult because we dont know the orientaton
    ! of the corners, and this which (i,j) dimension maps to which dimension
    !------------------------------------------------------------------------

    ! SWEST
    !----------
    do ll=swest, swest+max_corner_elem-1
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do ee=1,nc
          edge%buf(:,ee,kptr+kk,desc%putmapP_ghost(ll))=v(1:nc,ee,kk)
        end do
        end do
      endif
    end do
    
    ! SEAST
    !----------
    do ll=swest+max_corner_elem, swest+2*max_corner_elem-1
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do ee=1,nc
          edge%buf(ee,:,kptr+kk,desc%putmapP_ghost(ll))=v(nc-ee+1,1:nc,kk)
        end do
        end do
      endif
    end do
    
    ! NEAST
    !----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do ee=1,nc
        do ii=1,nc
          edge%buf(ii,ee,kptr+kk,desc%putmapP_ghost(ll))=v(nc-ii+1,nc-ee+1,kk)
        end do
        end do
        end do
      endif
    end do
    
    ! NWEST
    !----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do ee=1,nc
          edge%buf(:,ee,kptr+kk,desc%putmapP_ghost(ll))=v(1:nc,nc-ee+1,kk)
        end do
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine GhostVpackfull
  !==================================================================


  !==================================================================
  subroutine GhostVunpackfull(edge,v,nc1,nc2,nc,vlyr,kptr,desc)
    ! GhostVunpackfull:
    !
    ! edgeVunpack:
    !
    ! Unpack edges from edge buffer into v...
    !===================================================================
    use SE_Options,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !------------------
    integer               ,intent(in   ):: vlyr
    type(Ghostbuffer3D_t) ,intent(in   ):: edge
    real(kind=real_kind)  ,intent(inout):: v(nc1:nc2,nc1:nc2,vlyr)
    integer               ,intent(in   ):: nc1,nc2,nc
    integer               ,intent(in   ):: kptr
    type(EdgeDescriptor_t)              :: desc
    !
    ! Local Values
    !---------------
    logical,parameter:: UseUnroll = .TRUE.
    integer:: is,ie,in,iw,ic
    logical:: reverse
    integer:: ii,kk,ll,ee

    ! make sure buffer is big enough:
    !----------------------------------
    if((nc2-nc1+1) <  3*nc) then
      call endrun("GhostVunpack:  insufficient ghost cell region")
    endif

    threadsafe=.false.

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east )  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west )  

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1,k)
    ! 2nd   row ('edge') goes in v(:,np+2,k)
    ! etc...
    !------------------------------------------
    do ee=1,nc
    do kk=1,vlyr
    do ii=1,nc
      v(ii   , 1-ee,kk) = edge%buf(ii,ee,kptr+kk,is)
      v(nc+ee,ii   ,kk) = edge%buf(ii,ee,kptr+kk,ie)
      v(ii   ,nc+ee,kk) = edge%buf(ii,ee,kptr+kk,in)
      v( 1-ee,ii   ,kk) = edge%buf(ii,ee,kptr+kk,iw)
    end do
    end do
    end do

    ! SWEST
    !------------
    do ll=swest,swest+max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if((ic /= -1).and.(ll.eq.swest)) then   
      !if(ic /= -1) then     ??? WTF!
        ! unpack swest corner, IGNORE all other corners
        !-----------------------------------------------
        reverse=desc%reverse(ll)
        if(reverse) then
          do kk=1,vlyr
            !v(0  ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do ee=1,nc
            do ii=1,nc
              v(1-ee,1-ii,kk)=edge%buf(ii,ee,kptr+kk,ic)
            end do
            end do
          end do
        else
          do kk=1,vlyr
            !v(0  ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do ee=1,nc
            do ii=1,nc
              v(1-ee,1-ii,kk)=edge%buf(ee,ii,kptr+kk,ic)
            end do
            end do
          end do
        endif
      endif
    end do
    
    ! SEAST
    !-------------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if((ic /= -1).and.(ll.eq.seast)) then
      !if(ic /= -1) then 
        ! unpack seast corner, IGNORE all other corners
        !-----------------------------------------------
        reverse=desc%reverse(ll)
        if(reverse) then
          do kk=1,vlyr
            !v(nc+1 ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do ee=1,nc
            do ii=1,nc
              v(nc+ii,1-ee,kk)=edge%buf(ee,ii,kptr+kk,ic)
            end do
            end do
          end do
        else
          do kk=1,vlyr
            !v(nc+1 ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do ee=1,nc
            do ii=1,nc
              v(nc+ii,1-ee,kk)=edge%buf(ii,ee,kptr+kk,ic)
            end do
            end do
          end do
        endif
      endif
    end do

    ! NEAST
    !-----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if((ic /= -1).and.(ll.eq.neast)) then
     !if(ic /= -1) then 
        ! unpack neast corner, IGNORE all other corners
        !-----------------------------------------------
        reverse=desc%reverse(ll)
        if(reverse) then
          do kk=1,vlyr
            !v(nc+1,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do ee=1,nc
            do ii=1,nc
              v(nc+ii,nc+ee,kk)=edge%buf(ee,ii,kptr+kk,ic)
            end do
            end do
          end do
        else
          do kk=1,vlyr
            !v(nc+1,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do ee=1,nc
            do ii=1,nc
              v(nc+ii,nc+ee,kk)=edge%buf(ii,ee,kptr+kk,ic)
            end do
            end do
          end do
        endif
      endif
    end do

    ! NWEST
    !----------------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if((ic /= -1).and.(ll.eq.nwest)) then
        ! unpack nwest corner, IGNORE all other corners
        !-----------------------------------------------
        reverse=desc%reverse(ll)
        if(reverse) then
          do kk=1,vlyr
            !v(0  ,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do ee=1,nc
            do ii=1,nc
              v(1-ii,nc+ee,kk)=edge%buf(ee,ii,kptr+kk,ic)
            end do
            end do
          end do
        else
          do kk=1,vlyr
            !v(0  ,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do ee=1,nc
            do ii=1,nc
              v(1-ii,nc+ee,kk)=edge%buf(ii,ee,kptr+kk,ic)
            end do
            end do
          end do
        endif
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine GhostVunpackfull
  !==================================================================


  !==================================================================
  subroutine GhostVpack_unoriented(edge,v,nc,vlyr,kptr,desc)
    ! GhostVpack_unoriented:
    !
    !> @brief Pack edges of v into an edge buffer for boundary exchange.
    !
    !> This subroutine packs for one or more vertical layers into an edge 
    !! buffer. If the buffer associated with edge is not large enough to 
    !! hold all vertical layers you intent to pack, the method will 
    !! halt the program with a call to parallel_mod::endrun().
    !! @param[in] edge Ghost Buffer into which the data will be packed.
    !! This buffer must be previously allocated with initghostbufferfull().
    !! @param[in] v The data to be packed.
    !! @param[in] vlyr Number of vertical level coming into the subroutine
    !! for packing for input v.
    !! @param[in] kptr Vertical pointer to the place in the edge buffer where 
    !! data will be located.
    !=======================================================================
    use SE_Options,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !-------------------
    integer               ,intent(in):: vlyr
    type(Ghostbuffer3D_t)            :: edge
    real(kind=real_kind)  ,intent(in):: v(nc,nc,vlyr)
    integer               ,intent(in):: nc
    integer               ,intent(in):: kptr
    type(EdgeDescriptor_t),intent(in):: desc
    !
    ! Local Values
    !-------------
    integer:: is,ie,in,iw
    integer:: ii,kk,ir,ll,ee,l_local

    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
      threadsafe=.true.
    endif

#if 0
    do ll=1,max_neigh_edges
      is = desc%putmapP_ghost(ll)
      if(is /= -1) then       
        do kk=1,vlyr
          edge%buf(:,:,kptr+kk,is) = v(:,:,kk)
        end do
      endif
    end do
#else
    do l_local=1,desc%actual_neigh_edges
      ll=desc%loc2buf(l_local)
      is = desc%putmapP_ghost(ll)
      do kk=1,vlyr
        edge%buf(:,:,kptr+kk,is) =v(:,:,kk)  
      end do
    end do
#endif

    ! End Routine 
    !-------------
    return
  end subroutine GhostVpack_unoriented
  !==================================================================


  !==================================================================
  subroutine GhostVunpack_unoriented(edge,v,nc,vlyr,kptr,desc)
    ! GhostVunpack_unoriented:
    !
    ! edgeVunpack:
    !
    ! Unpack edges from edge buffer into v...
    !===================================================================
    use SE_Options,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Varaibles
    !-------------------
    integer               ,intent(in   ):: vlyr
    type(Ghostbuffer3D_t) ,intent(in   ):: edge
    real(kind=real_kind)  ,intent(inout):: v(nc,nc,vlyr,*)
    integer               ,intent(in   ):: nc
    integer               ,intent(in   ):: kptr
    type(EdgeDescriptor_t)              :: desc
    !
    ! Local Values
    !---------------
    logical,parameter:: UseUnroll = .TRUE.
    integer:: l_local
    integer:: ii,kk,ll,ee,is

    threadsafe=.false.

#if 0
    do ll=1,max_neigh_edges
      is = desc%getmapP_ghost(ll)
      if(is  /= -1) then       
        l_local = ll
        do kk=1,vlyr
          v(:,:,kk,l_local) = edge%buf(:,:,kptr+kk,is) 
        end do
      endif
    end do
#endif

    do l_local=1,desc%actual_neigh_edges
      ll=desc%loc2buf(l_local)
      is = desc%getmapP_ghost(ll)
      do kk=1,vlyr
        v(:,:,kk,l_local ) = edge%buf(:,:,kptr+kk,is) 
      end do
    end do

    ! End Routine 
    !-------------
    return
  end subroutine GhostVunpack_unoriented
  !==================================================================


  !==================================================================
  subroutine initGhostBufferTR(ghost,nlyr,ntrac,nhc, npoints)
    ! initGhostBufferTR:
    !
    ! initGhostBuffer:
    ! Author: Christoph Erath
    ! create an Real based communication buffer
    ! npoints is the number of points on one side
    ! nhc is the deep of the ghost/halo zone
    !======================================================================
    !
    ! Passed Variables
    !------------------
    type(Ghostbuffertr_t),intent(out):: ghost
    integer              ,intent(in ):: nlyr,ntrac,nhc, npoints
    !
    ! Local Values
    !---------------
    integer :: nbuf

    ! make sure buffer is big enough:
    !--------------------------------
    if(nhc > npoints) then
      call endrun("intGhostBuffer:  halo region can not be larger then element size")
    endif
    if(ntrac < 1) then
      call endrun("intGhostBuffer:  you have to consider at least one tracer")
    endif

    ! sanity check for threading
    !-----------------------------
    if(omp_get_num_threads() > 1) then
      call endrun('ERROR: initGhostBuffer must be called before threaded region')
    endif

    nbuf=max_neigh_edges*nelemd
#if (! defined ELEMENT_OPENMP)
!$OMP MASTER
#endif
    ghost%nlyr=nlyr
    ghost%nbuf=nbuf
    allocate(ghost%buf(npoints,nhc,nlyr,ntrac,nbuf))
    ghost%buf=0

    allocate(ghost%receive(npoints,nhc,nlyr,ntrac,nbuf))
    ghost%receive=0
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif

    ! End Routine 
    !-------------
    return
  end subroutine initGhostBufferTR
  !==================================================================


  !==================================================================
  subroutine ghostVpack(edge,v,nhc,npoints,vlyr,ntrac,kptr,tn0,timelevels,desc)
    ! ghostVpack:
    !
    ! Christoph Erath
    !> Packs the halo zone from v
    ! =========================================
    ! NOTE: I have to give timelevels as argument, because element_mod is not compiled first
    ! and the array call has to be done in this way because of performance reasons!!!
    !===========================================================================
    use SE_Constants,only: ntrac_d
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !---------------------
    integer               ,intent(in):: vlyr
    integer               ,intent(in):: ntrac
    integer               ,intent(in):: nhc,npoints
    integer               ,intent(in):: kptr, tn0, timelevels
    type(Ghostbuffertr_t)            :: edge
    real(kind=real_kind)  ,intent(in):: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac_d,timelevels)
    type(EdgeDescriptor_t),intent(in):: desc
    !
    ! Local values
    !----------------
    integer :: is,ie,in,iw
    integer :: ii,jj,kk,ir,ll,itr

    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
      threadsafe=.true.
    endif

    ! Example convenction for buffer to the north:
    !  buf(:,,:,i,e) 
    !   each "edge" is a row of data (i=1,np) in the element 
    !     north most row of data goes into e=1
    !     next row of data goes into e=2
    !     .... 
    !     south most row of data goes into e=np
    !   We need to pack this way to preserve the orientation
    !   so the data can be unpacked correctly
    !
    ! note: we think of buf as dimensioned buf(k,is,i,e)
    ! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
    !---------------------------------------------------------------
    is = desc%putmapP_ghost(south) 
    ie = desc%putmapP_ghost(east )  
    in = desc%putmapP_ghost(north) 
    iw = desc%putmapP_ghost(west )  

!   print *,nc,is,ie,in,iw
    do itr=1,ntrac
    do kk=1,vlyr
    do ii=1,npoints
    do jj=1,nhc
      edge%buf(ii,jj,kptr+kk,itr,is) = v(ii          ,jj          ,kk,itr,tn0)
      edge%buf(ii,jj,kptr+kk,itr,ie) = v(npoints-jj+1,ii          ,kk,itr,tn0)
      edge%buf(ii,jj,kptr+kk,itr,in) = v(ii          ,npoints-jj+1,kk,itr,tn0)
      edge%buf(ii,jj,kptr+kk,itr,iw) = v(jj          ,ii          ,kk,itr,tn0)
    end do
    end do
    end do
    end do

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !  Check if the edge orientation of the recieving element is different
    !  if it is, swap the order of data in the edge
    !-----------------------------------------------------------------------------
    if(desc%reverse(south)) then
!     is = desc%putmapP_ghost(south)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        edge%buf(ir,jj,kptr+kk,itr,is)=v(ii,jj,kk,itr,tn0)
      end do
      end do
      end do
      end do
    endif

    if(desc%reverse(east)) then
!     ie = desc%putmapP_ghost(east)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        edge%buf(ir,jj,kptr+kk,itr,ie)=v(npoints-jj+1,ii,kk,itr,tn0)
      end do
      end do
      end do
      end do
    endif

    if(desc%reverse(north)) then
!     in = desc%putmapP_ghost(north)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        edge%buf(ir,jj,kptr+kk,itr,in)=v(ii,npoints-jj+1,kk,itr,tn0)
      end do
      end do
      end do
      end do
    endif

    if(desc%reverse(west)) then
!     iw = desc%putmapP_ghost(west)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        edge%buf(ir,jj,kptr+kk,itr,iw)=v(jj,ii,kk,itr,tn0)
      end do
      end do
      end do
      end do
    endif

    ! corners.  this is difficult because we dont know the orientaton
    ! of the corners, and this which (i,j) dimension maps to which dimension
    !-----------------------------------------------------------------------

    ! SWEST
    !------------
    do ll=swest,swest+max_corner_elem-1
      if(ll.ne.swest) call endrun('ERROR2: swest ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
          ! edge%buf(1,1,kptr+kk,desc%putmapP_ghost(ll))=v(1,1 ,kk)
          do ii=1,nhc
          do jj=1,nhc
            edge%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(ii,jj,kk,itr,tn0)
          end do
          end do
        end do
        end do
      endif
    end do
  
    ! SEAST
    !------------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(ll.ne.seast) call endrun('ERROR2: seast ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
          ! edge%buf(1,1,kptr+kk,desc%putmapP_ghost(ll))=v(nc ,1 ,kk)
          do ii=1,nhc
          do jj=1,nhc
            edge%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(npoints-ii+1,jj,kk,itr,tn0)
          end do
          end do
        end do
        end do
      endif
    end do
  
    ! NEAST
    !------------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(ll.ne.neast) call endrun('ERROR2: neast ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
          !edge%buf(1,1,kptr+kk,desc%putmapP_ghost(ll))=v(nc ,nc,kk)
          do ii=1,nhc
          do jj=1,nhc
            edge%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(npoints-ii+1,npoints-jj+1,kk,itr,tn0)
          end do
          end do
        end do
        end do
      endif
    end do
  
    ! NWEST
    !------------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(ll.ne.nwest) call endrun('ERROR2: nwest ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
          !edge%buf(1,1,kptr+kk,desc%putmapP_ghost(ll))=v(1  ,nc,kk)
          do ii=1,nhc
          do jj=1,nhc
            edge%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(ii,npoints-jj+1,kk,itr,tn0)
          end do
          end do
        end do
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine GhostVpack
  !==================================================================


  !==================================================================
  subroutine ghostVpackR(edge,v,nhc,npoints,vlyr,ntrac,kptr,desc)
    ! ghostVpackR:
    !
    ! Christoph Erath
    !> Packs the halo zone from v
    ! =========================================
    ! NOTE: I have to give timelevels as argument, because element_mod is not compiled first
    ! and the array call has to be done in this way because of performance reasons!!!
    !=====================================================================================
    use SE_Constants,only: ntrac_d
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !---------------------
    integer               ,intent(in):: vlyr
    integer               ,intent(in):: nhc,npoints
    integer               ,intent(in):: ntrac
    type(Ghostbuffertr_t)            :: edge
    real(kind=real_kind)  ,intent(in):: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac_d)
    integer               ,intent(in):: kptr
    type(EdgeDescriptor_t),intent(in):: desc
    !
    ! Local Values
    !------------------
    integer:: is,ie,in,iw
    integer:: ii,jj,kk,ir,ll,itr

    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
      threadsafe=.true.
    endif

    ! Example convenction for buffer to the north:
    !  buf(:,,:,i,e) 
    !   each "edge" is a row of data (i=1,np) in the element 
    !     north most row of data goes into e=1
    !     next row of data goes into e=2
    !     .... 
    !     south most row of data goes into e=np
    !   We need to pack this way to preserve the orientation
    !   so the data can be unpacked correctly
    !
    ! note: we think of buf as dimensioned buf(k,is,i,e)
    ! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
    !------------------------------------------------------------
    is = desc%putmapP_ghost(south) 
    ie = desc%putmapP_ghost(east )  
    in = desc%putmapP_ghost(north) 
    iw = desc%putmapP_ghost(west )  

!   print *,nc,is,ie,in,iw
    do itr=1,ntrac
    do kk=1,vlyr
    do ii=1,npoints
    do jj=1,nhc
      edge%buf(ii,jj,kptr+kk,itr,is) = v(ii          ,jj          ,kk,itr)
      edge%buf(ii,jj,kptr+kk,itr,ie) = v(npoints-jj+1,ii          ,kk,itr)
      edge%buf(ii,jj,kptr+kk,itr,in) = v(ii          ,npoints-jj+1,kk,itr)
      edge%buf(ii,jj,kptr+kk,itr,iw) = v(jj          ,ii          ,kk,itr)
    end do
    end do
    end do
    end do

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !  Check if the edge orientation of the recieving element is different
    !  if it is, swap the order of data in the edge
    !----------------------------------------------------------------------
    if(desc%reverse(south)) then
!     is = desc%putmapP_ghost(south)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        edge%buf(ir,jj,kptr+kk,itr,is)=v(ii,jj,kk,itr)
      end do
      end do
      end do
      end do
    endif

    if(desc%reverse(east)) then
!     ie = desc%putmapP_ghost(east)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        edge%buf(ir,jj,kptr+kk,itr,ie)=v(npoints-jj+1,ii,kk,itr)
      end do
      end do
      end do
      end do
    endif

    if(desc%reverse(north)) then
!     in = desc%putmapP_ghost(north)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        edge%buf(ir,jj,kptr+kk,itr,in)=v(ii,npoints-jj+1,kk,itr)
      end do
      end do
      end do
      end do
    endif

    if(desc%reverse(west)) then
!     iw = desc%putmapP_ghost(west)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        edge%buf(ir,jj,kptr+kk,itr,iw)=v(jj,ii,kk,itr)
      end do
      end do
      end do
      end do
    endif

    ! corners.  this is difficult because we dont know the orientaton
    ! of the corners, and this which (i,j) dimension maps to which dimension
    !--------------------------------------------------------------------------

    ! SWEST
    !----------
    do ll=swest,swest+max_corner_elem-1
      if(ll.ne.swest) call endrun('ERROR2: swest ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
          ! edge%buf(1,1,kptr+kk,desc%putmapP_ghost(ll))=v(1,1 ,kk)
          do ii=1,nhc
          do jj=1,nhc
            edge%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(ii,jj,kk,itr)
          end do
          end do
        end do
        end do
      endif
    end do
  
    ! SEAST
    !----------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(ll.ne.seast) call endrun('ERROR2: seast ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
          ! edge%buf(1,1,kptr+kk,desc%putmapP_ghost(ll))=v(nc ,1 ,kk)
          do ii=1,nhc
          do jj=1,nhc
            edge%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(npoints-ii+1,jj,kk,itr)
          end do
          end do
        end do
        end do
      endif
    end do
  
    ! NEAST
    !----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(ll.ne.neast) call endrun('ERROR2: neast ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
          !edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,nc,k)
          do ii=1,nhc
          do jj=1,nhc
            edge%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(npoints-ii+1,npoints-jj+1,kk,itr)
          end do
          end do
        end do
        end do
      endif
    end do
  
    ! NWEST
    !----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(ll.ne.nwest) call endrun('ERROR2: nwest ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
          !edge%buf(1,1,kptr+kk,desc%putmapP_ghost(ll))=v(1  ,nc,kk)
          do ii=1,nhc
          do jj=1,nhc
            edge%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(ii,npoints-jj+1,kk,itr)
          end do
          end do
        end do
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine GhostVpackR
  !==================================================================


  !==================================================================
  subroutine ghostVunpack(edge,v,nhc,npoints,vlyr,ntrac,kptr,tn0,timelevels,desc)
    ! ghostVunpack:
    !
    ! Christoph Erath
    !
    ! Unpack the halo zone into v
    ! ========================================
    ! NOTE: I have to give timelevels as argument, because element_mod is not compiled first
    ! and the array call has to be done in this way because of performance reasons!!!
    !===================================================================================
    use SE_Constants,only: ntrac_d
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !-------------------
    integer              ,intent(in   ):: vlyr
    integer              ,intent(in   ):: nhc,npoints
    integer              ,intent(in   ):: ntrac
    integer              ,intent(in   ):: kptr, tn0, timelevels
    type(Ghostbuffertr_t),intent(in   ):: edge
    real(kind=real_kind) ,intent(inout):: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac_d,timelevels)
    type(EdgeDescriptor_t)             :: desc
    !
    ! Local Values
    !-----------------
    logical,parameter   :: UseUnroll = .TRUE.
    real(kind=real_kind):: NaN=-1.0
    integer:: is,ie,in,iw,ic
    logical:: reverse
    integer:: ii,jj,kk,ll,itr

    NaN=sqrt(NaN)

    threadsafe=.false.

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east )  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west )  

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1,k)
    ! 2nd   row ('edge') goes in v(:,np+2,k)
    ! etc...
    !------------------------------------------
    do itr=1,ntrac
    do kk=1,vlyr
    do ii=1,npoints
    do jj=1,nhc
      v(ii        ,1-jj      ,kk,itr,tn0) = edge%buf(ii,jj,kptr+kk,itr,is)
      v(npoints+jj,ii        ,kk,itr,tn0) = edge%buf(ii,jj,kptr+kk,itr,ie)
      v(ii        ,npoints+jj,kk,itr,tn0) = edge%buf(ii,jj,kptr+kk,itr,in)
      v(1-jj      ,ii        ,kk,itr,tn0) = edge%buf(ii,jj,kptr+kk,itr,iw)
    end do
    end do
    end do
    end do

    ! SWEST
    !----------
    do ll=swest,swest+max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
            !v(0  ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(1-jj,1-ii,kk,itr,tn0)=edge%buf(ii,jj,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
            !v(0  ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(1-jj,1-ii,kk,itr,tn0)=edge%buf(jj,ii,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(1-ii,1-jj,kk,itr,tn0)=NaN
        end do
        end do
        end do
        end do    
      endif
    end do
  
    ! SEAST
    !----------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
            !v(nc+1 ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(npoints+ii,1-jj,kk,itr,tn0)=edge%buf(jj,ii,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
            !v(nc+1 ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(npoints+ii,1-jj,kk,itr,tn0)=edge%buf(ii,jj,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(npoints+ii,1-jj,kk,itr,tn0)=NaN
        end do
        end do
        end do
        end do  
      endif
    end do

    ! NEAST
    !----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
            !v(nc+1 ,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(npoints+ii,npoints+jj,kk,itr,tn0)=edge%buf(jj,ii,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
            !v(nc+1 ,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(npoints+ii,npoints+jj,kk,itr,tn0)=edge%buf(ii,jj,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr        
        do jj=1,nhc
        do ii=1,nhc
          v(npoints+ii,npoints+jj,kk,itr,tn0)=NaN
        end do
        end do
        end do
        end do    
      endif
    end do

    ! NWEST
    !----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
            !v(0  ,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(1-ii,npoints+jj,kk,itr,tn0)=edge%buf(jj,ii,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
            !v(0  ,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(1-ii,npoints+jj,kk,itr,tn0)=edge%buf(ii,jj,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(1-ii,npoints+jj,kk,itr,tn0)=NaN
        end do
        end do
        end do
        end do    
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine ghostVunpack
  !==================================================================


  !==================================================================
  subroutine ghostVunpackR(edge,v,nhc,npoints,vlyr,ntrac,kptr,desc)
    ! ghostVunpackR:
    !
    ! Christoph Erath
    !
    ! Unpack the halo zone into v
    ! ========================================
    ! NOTE: I have to give timelevels as argument, because element_mod is not compiled first
    ! and the array call has to be done in this way because of performance reasons!!!
    !==================================================================================
    use SE_Constants,only: ntrac_d
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !----------------------
    integer               ,intent(in   ):: vlyr
    integer               ,intent(in   ):: nhc,npoints
    integer               ,intent(in   ):: ntrac
    type (Ghostbuffertr_t),intent(in   ):: edge
    real (kind=real_kind) ,intent(inout):: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac_d)
    integer               ,intent(in   ):: kptr
    type (EdgeDescriptor_t)             :: desc
    !
    ! Local Values
    !-----------------
    logical,parameter   :: UseUnroll = .TRUE.
    real(kind=real_kind):: NaN=-1.0
    integer:: is,ie,in,iw,ic
    logical:: reverse
    integer:: ii,jj,kk,ll,itr
  
    NaN=sqrt(NaN)

    threadsafe=.false.

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east )  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west )  

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1,k)
    ! 2nd   row ('edge') goes in v(:,np+2,k)
    ! etc...
    !-------------------------------------------
    do itr=1,ntrac
    do kk=1,vlyr
    do ii=1,npoints
    do jj=1,nhc
      v(ii        ,1-jj      ,kk,itr) = edge%buf(ii,jj,kptr+kk,itr,is)
      v(npoints+jj,ii        ,kk,itr) = edge%buf(ii,jj,kptr+kk,itr,ie)
      v(ii        ,npoints+jj,kk,itr) = edge%buf(ii,jj,kptr+kk,itr,in)
      v(1-jj      ,ii        ,kk,itr) = edge%buf(ii,jj,kptr+kk,itr,iw)
    end do
    end do
    end do
    end do

    ! SWEST
    !----------
    do ll=swest,swest+max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
            !v(0  ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(1-jj,1-ii,kk,itr)=edge%buf(ii,jj,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
            !v(0  ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(1-jj,1-ii,kk,itr)=edge%buf(jj,ii,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(1-ii,1-jj,kk,itr)=NaN
        end do
        end do
        end do
        end do    
      endif
    end do
  
    ! SEAST
    !----------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
            !v(nc+1 ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(npoints+ii,1-jj,kk,itr)=edge%buf(jj,ii,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
            !v(nc+1 ,0 ,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(npoints+ii,1-jj,kk,itr)=edge%buf(ii,jj,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(npoints+ii,1-jj,kk,itr)=NaN
        end do
        end do
        end do
        end do  
      endif
    end do

    ! NEAST
    !----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
            !v(nc+1 ,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(npoints+ii,npoints+jj,kk,itr)=edge%buf(jj,ii,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
            !v(nc+1 ,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(npoints+ii,npoints+jj,kk,itr)=edge%buf(ii,jj,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr        
        do jj=1,nhc
        do ii=1,nhc
          v(npoints+ii,npoints+jj,kk,itr)=NaN
        end do
        end do
        end do
        end do    
      endif
    end do

    ! NWEST
    !----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
            !v(0  ,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(1-ii,npoints+jj,kk,itr)=edge%buf(jj,ii,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
            !v(0  ,nc+1,kk)=edge%buf(1,1,kptr+kk,desc%getmapP_ghost(ll))
            do jj=1,nhc
            do ii=1,nhc
              v(1-ii,npoints+jj,kk,itr)=edge%buf(ii,jj,kptr+kk,itr,ic)
            end do
            end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(1-ii,npoints+jj,kk,itr)=NaN
        end do
        end do
        end do
        end do    
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine ghostVunpackR
  !==================================================================


  !==================================================================
  subroutine ghostVpack2d(ghost,v,nhc, npoints,vlyr,ntrac,kptr,tn0,timelevels,desc)
    ! ghostVpack2d:
    !
    ! GHOSTVPACK2D
    ! AUTHOR: Christoph Erath 
    ! Pack edges of v into an ghost buffer for boundary exchange.
    !
    ! This subroutine packs for one vertical layers into an ghost
    ! buffer. It is for cartesian points (v is only two dimensional). 
    ! If the buffer associated with edge is not large enough to 
    ! hold all vertical layers you intent to pack, the method will 
    ! halt the program with a call to parallel_mod::endrun().
    ! INPUT: 
    ! - ghost Buffer into which the data will be packed.
    !   This buffer must be previously allocated with initGhostBuffer().
    ! - v The data to be packed.
    ! - nhc deep of ghost/halo zone
    ! - npoints number of points on on side
    ! - kptr Vertical pointer to the place in the edge buffer where 
    ! data will be located.
    !=================================================================================
    use SE_Constants,only: ntrac_d
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !--------------------
    integer               ,intent(in):: nhc,npoints
    integer               ,intent(in):: vlyr
    integer               ,intent(in):: ntrac
    integer               ,intent(in):: kptr, tn0, timelevels
    type(Ghostbuffertr_t)            :: ghost
    real(kind=real_kind)  ,intent(in):: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr, ntrac_d,timelevels)
    type(EdgeDescriptor_t),intent(in):: desc
    !
    ! Local Values
    !----------------
    integer:: itr,kk
    integer:: is,ie,in,iw
    integer:: ii,jj,ir,ll,ee

    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
  !$OMP BARRIER
#endif
      threadsafe=.true.
    endif

    ! Example convenction for buffer to the north:
    !  buf(:,,:,i,e) 
    !   each "edge" is a row of data (i=1,np) in the element 
    !     north most row of data goes into e=1
    !     next row of data goes into e=2
    !     .... 
    !     south most row of data goes into e=np
    !   We need to pack this way to preserve the orientation
    !   so the data can be unpacked correctly
    !
    ! note: we think of buf as dimensioned buf(k,is,i,e)
    ! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
    !----------------------------------------------------------

    is = desc%putmapP_ghost(south) 
    ie = desc%putmapP_ghost(east )  
    in = desc%putmapP_ghost(north) 
    iw = desc%putmapP_ghost(west )  

    do itr=1,ntrac
    do kk=1,vlyr
    do ii=1,npoints
    do jj=1,nhc
      ghost%buf(ii,jj,kptr+kk,itr,is) = v(ii        ,jj+1      ,kk,itr,tn0)
      ghost%buf(ii,jj,kptr+kk,itr,ie) = v(npoints-jj,ii        ,kk,itr,tn0)
      ghost%buf(ii,jj,kptr+kk,itr,in) = v(ii        ,npoints-jj,kk,itr,tn0)
      ghost%buf(ii,jj,kptr+kk,itr,iw) = v(jj+1      ,ii        ,kk,itr,tn0)
    end do
    end do
    end do
    end do

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !  Check if the edge orientation of the recieving element is different
    !  if it is, swap the order of data in the edge
    !------------------------------------------------------------------------
    if(desc%reverse(south)) then
  !   is = desc%putmapP_ghost(south)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,kptr+kk,itr,is)=v(ii,jj+1,kk,itr,tn0)
      end do
      end do
      end do
      end do
    endif

    if(desc%reverse(east)) then
  !   ie = desc%putmapP_ghost(east)
      do itr=1,ntrac
      do kk=1,vlyr  
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,kptr+kk,itr,ie)=v(npoints-jj,ii,kk,itr,tn0)
      end do
      end do
      end do
      end do
    endif

    if(desc%reverse(north)) then
  !   in = desc%putmapP_ghost(north)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,kptr+kk,itr,in)=v(ii,npoints-jj,kk,itr,tn0)
      end do
      end do
      end do
      end do
    endif

    if(desc%reverse(west)) then
  !   iw = desc%putmapP_ghost(west)
      do itr=1,ntrac
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,kptr+kk,itr,iw)=v(jj+1,ii,kk,itr,tn0)
      end do
      end do
      end do
      end do
    endif

    ! corners.  this is difficult because we dont know the orientaton
    ! of the corners, and this which (i,j) dimension maps to which dimension
    !------------------------------------------------------------------------

    ! SWEST
    !-----------
    do ll=swest,swest+max_corner_elem-1
      if(ll.ne.swest) call endrun('ERROR3: swest ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
        do ii=1,nhc
        do jj=1,nhc
          ghost%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(ii+1,jj+1,kk,itr,tn0)
        end do
        end do     
        end do
        end do        
      endif
    end do

    ! SEAST
    !-----------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(ll.ne.seast) call endrun('ERROR3: seast ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
        do ii=1,nhc             
        do jj=1,nhc
          ghost%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(npoints-ii,jj+1,kk,itr,tn0)
        end do
        end do
        end do
        end do             
      endif
    end do

    ! NEAST
    !-----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(ll.ne.neast) call endrun('ERROR3: neast ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
        do ii=1,nhc
        do jj=1,nhc
          ghost%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(npoints-ii,npoints-jj,kk,itr,tn0)           
        end do
        end do
        end do
        end do
      endif
    end do

    ! NWEST
    !-----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(ll.ne.nwest) call endrun('ERROR3: nwest ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do itr=1,ntrac
        do kk=1,vlyr
        do ii=1,nhc
        do jj=1,nhc
          ghost%buf(ii,jj,kptr+kk,itr,desc%putmapP_ghost(ll))=v(ii+1,npoints-jj,kk,itr,tn0)
        end do       
        end do
        end do
        end do      
      endif
    end do   

    ! End Routine 
    !-------------
    return
  end subroutine ghostVpack2d
  !==================================================================


  !==================================================================
  subroutine ghostVunpack2d(ghost,v,nhc,npoints,vlyr,ntrac,kptr,tn0,timelevels,desc)
    ! ghostVunpack2d:
    !
    ! GHOSTVUNPACK2D
    ! AUTHOR: Christoph Erath 
    ! Unpack ghost points from ghost buffer into v...
    ! It is for cartesian points (v is only two dimensional).
    ! INPUT SAME arguments as for GHOSTVPACK2d
    !=================================================================================
    use SE_Constants,only: ntrac_d
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    !  Passed Variables
    !--------------------
    integer              ,intent(in   ):: nhc,npoints
    integer              ,intent(in   ):: vlyr
    integer              ,intent(in   ):: ntrac
    integer              ,intent(in   ):: kptr, tn0,timelevels
    type(Ghostbuffertr_t),intent(in   ):: ghost
    real(kind=real_kind) ,intent(inout):: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr, ntrac_d,timelevels)
    type(EdgeDescriptor_t)             :: desc
    !
    ! Local Values
    !----------------
    logical,parameter   :: UseUnroll = .TRUE.
    real(kind=real_kind):: NaN=-1.0
    integer:: is,ie,in,iw,ic
    logical:: reverse
    integer:: ii,jj,ll,itr,kk
    
    NaN=sqrt(NaN)

    threadsafe=.false.

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east )  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west )  

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1)
    ! 2nd   row ('edge') goes in v(:,np+2)
    ! etc...
    !-----------------------------------------
    do itr=1,ntrac
    do kk=1,vlyr
    do ii=1,npoints
    do jj=1,nhc
      v(ii        ,1-jj      ,kk,itr,tn0) = ghost%buf(ii,jj,kptr+kk,itr,is)
      v(npoints+jj,ii        ,kk,itr,tn0) = ghost%buf(ii,jj,kptr+kk,itr,ie)
      v(ii        ,npoints+jj,kk,itr,tn0) = ghost%buf(ii,jj,kptr+kk,itr,in)
      v(1-jj      ,ii        ,kk,itr,tn0) = ghost%buf(ii,jj,kptr+kk,itr,iw)
    end do
    end do
    end do
    end do

    ! SWEST
    !------------
    do ll=swest,swest+max_corner_elem-1
      ic = desc%getmapP_ghost(ll)       
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,1-jj,kk,itr,tn0)=ghost%buf(jj,ii,kptr+kk,itr,ic)
          end do
          end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,1-jj,kk,itr,tn0)=ghost%buf(ii,jj,kptr+kk,itr,ic)
          end do
          end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(1-ii,1-jj,kk,itr,tn0)=NaN
        end do
        end do
        end do
        end do      
      endif
    end do
    
    ! SEAST
    !------------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,1-jj,kk,itr,tn0)=ghost%buf(jj,ii,kptr+kk,itr,ic)
          end do
          end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,1-jj,kk,itr,tn0)=ghost%buf(ii,jj,kptr+kk,itr,ic)
          end do
          end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(npoints+ii,1-jj,kk,itr,tn0)=NaN
        end do
        end do
        end do
        end do
      endif
    end do

    ! NEAST
    !------------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,npoints+jj,kk,itr,tn0)=ghost%buf(jj,ii,kptr+kk,itr,ic)
          end do
          end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,npoints+jj,kk,itr,tn0)=ghost%buf(ii,jj,kptr+kk,itr,ic)
          end do
          end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(npoints+ii,npoints+jj,kk,itr,tn0)=NaN
        end do
        end do 
        end do
        end do   
      endif
    end do

    ! NWEST
    !------------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do itr=1,ntrac
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,npoints+jj,kk,itr,tn0)=ghost%buf(jj,ii,kptr+kk,itr,ic)
          end do
          end do
          end do
          end do
        else
          do itr=1,ntrac
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,npoints+jj,kk,itr,tn0)=ghost%buf(ii,jj,kptr+kk,itr,ic)
          end do
          end do
          end do
          end do
        endif
      else
        do itr=1,ntrac
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(1-ii,npoints+jj,kk,itr,tn0)=NaN
        end do
        end do
        end do
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine ghostVunpack2d
  !==================================================================


  !==================================================================
  subroutine ghostVpack2d_level(ghost,v,kptr,nhc, npoints,vlyr,desc)
    ! ghostVpack2d_level:
    !
    ! GHOSTVPACK2D
    ! AUTHOR: Christoph Erath 
    ! Pack edges of v into an ghost buffer for boundary exchange.
    !
    ! This subroutine packs for one vertical layers into an ghost
    ! buffer. It is for cartesian points (v is only two dimensional). 
    ! If the buffer associated with edge is not large enough to 
    ! hold all vertical layers you intent to pack, the method will 
    ! halt the program with a call to parallel_mod::endrun().
    ! INPUT: 
    ! - ghost Buffer into which the data will be packed.
    !   This buffer must be previously allocated with initGhostBuffer().
    ! - v The data to be packed.
    ! - nhc deep of ghost/halo zone
    ! - npoints number of points on on side
    ! - kptr Vertical pointer to the place in the edge buffer where 
    ! data will be located.
    !=================================================================================
    use SE_Constants,only: ntrac_d
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !--------------------
    integer               ,intent(in):: kptr,nhc,npoints
    integer               ,intent(in):: vlyr
    type(Ghostbuffertr_t)            :: ghost
    real(kind=real_kind)  ,intent(in):: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr)
    type(EdgeDescriptor_t),intent(in):: desc
    !
    ! Local variables
    !--------------------
    integer :: itr,kk
    integer :: is,ie,in,iw
    integer :: ii,jj,ir,ll,ee

    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
  !$OMP BARRIER
#endif
      threadsafe=.true.
    endif

    ! Example convenction for buffer to the north:
    !  buf(:,,:,i,e) 
    !   each "edge" is a row of data (i=1,np) in the element 
    !     north most row of data goes into e=1
    !     next row of data goes into e=2
    !     .... 
    !     south most row of data goes into e=np
    !   We need to pack this way to preserve the orientation
    !   so the data can be unpacked correctly
    !
    ! note: we think of buf as dimensioned buf(k,is,i,e)
    ! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
    !-----------------------------------------------------------
    is = desc%putmapP_ghost(south) 
    ie = desc%putmapP_ghost(east )  
    in = desc%putmapP_ghost(north) 
    iw = desc%putmapP_ghost(west )  

    do kk=1,vlyr
    do ii=1,npoints
    do jj=1,nhc
      ghost%buf(ii,jj,kk,kptr,is) = v(ii        ,jj+1      ,kk)
      ghost%buf(ii,jj,kk,kptr,ie) = v(npoints-jj,ii        ,kk)
      ghost%buf(ii,jj,kk,kptr,in) = v(ii        ,npoints-jj,kk)
      ghost%buf(ii,jj,kk,kptr,iw) = v(jj+1      ,ii        ,kk)
    end do
    end do
    end do

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !  Check if the edge orientation of the recieving element is different
    !  if it is, swap the order of data in the edge
    !-------------------------------------------------------------------------
    if(desc%reverse(south)) then
  !   is = desc%putmapP_ghost(south)
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,kk,kptr,is)=v(ii,jj+1,kk)
      end do
      end do
      end do
    endif

    if(desc%reverse(east)) then
  !   ie = desc%putmapP_ghost(east)
      do kk=1,vlyr  
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,kk,kptr,ie)=v(npoints-jj,ii,kk)
      end do
      end do
      end do
    endif

    if(desc%reverse(north)) then
  !   in = desc%putmapP_ghost(north)
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,kk,kptr,in)=v(ii,npoints-jj,kk)
      end do
      end do
      end do
    endif

    if(desc%reverse(west)) then
  !   iw = desc%putmapP_ghost(west)
      do kk=1,vlyr
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,kk,kptr,iw)=v(jj+1,ii,kk)
      end do
      end do
      end do
    endif

    ! corners.  this is difficult because we dont know the orientaton
    ! of the corners, and this which (i,j) dimension maps to which dimension
    !--------------------------------------------------------------------------

    ! SWEST
    !-------------
    do ll=swest,swest+max_corner_elem-1
      if(ll.ne.swest) call endrun('ERROR3: swest ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do ii=1,nhc
        do jj=1,nhc
          ghost%buf(ii,jj,kk,kptr,desc%putmapP_ghost(ll))=v(ii+1,jj+1,kk)
        end do
        end do     
        end do
      endif
    end do
  
    ! SEAST
    !-------------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(ll.ne.seast) call endrun('ERROR3: seast ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do ii=1,nhc             
        do jj=1,nhc
          ghost%buf(ii,jj,kk,kptr,desc%putmapP_ghost(ll))=v(npoints-ii,jj+1,kk)
        end do
        end do
        end do
      endif
    end do

    ! NEAST
    !-------------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(ll.ne.neast) call endrun('ERROR3: neast ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do ii=1,nhc
        do jj=1,nhc
          ghost%buf(ii,jj,kk,kptr,desc%putmapP_ghost(ll))=v(npoints-ii,npoints-jj,kk)           
        end do
        end do
        end do
      endif
    end do

    ! NWEST
    !-------------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(ll.ne.nwest) call endrun('ERROR3: nwest ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do ii=1,nhc
        do jj=1,nhc
          ghost%buf(ii,jj,kk,kptr,desc%putmapP_ghost(ll))=v(ii+1,npoints-jj,kk)
        end do       
        end do
        end do
      endif
    end do   

    ! End Routine 
    !-------------
    return
  end subroutine ghostVpack2d_level
  !==================================================================


  !==================================================================
  subroutine ghostVunpack2d_level(ghost,v,kptr,nhc,npoints,vlyr,desc)
    ! ghostVunpack2d_level:
    !
    ! GHOSTVUNPACK2D
    ! AUTHOR: Christoph Erath 
    ! Unpack ghost points from ghost buffer into v...
    ! It is for cartesian points (v is only two dimensional).
    ! INPUT SAME arguments as for GHOSTVPACK2d
    !=================================================================================
    use SE_Constants,only: ntrac_d
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !--------------------
    integer              ,intent(in   ):: kptr,nhc,npoints,vlyr
    type(Ghostbuffertr_t),intent(in   ):: ghost
    real(kind=real_kind) ,intent(inout):: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr)
    type(EdgeDescriptor_t)             :: desc
    !
    ! Local Values
    !--------------
    logical,parameter   :: UseUnroll = .TRUE.
    real(kind=real_kind):: NaN=-1.0
    integer:: is,ie,in,iw,ic
    logical:: reverse
    integer:: ii,jj,ll,itr,kk
    
    NaN=sqrt(NaN)

    threadsafe=.false.

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east )  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west )  

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1)
    ! 2nd   row ('edge') goes in v(:,np+2)
    ! etc...
    !----------------------------------------------
    do kk=1,vlyr
    do ii=1,npoints
    do jj=1,nhc
      v(ii        ,1-jj      ,kk) = ghost%buf(ii,jj,kk,kptr,is)
      v(npoints+jj,ii        ,kk) = ghost%buf(ii,jj,kk,kptr,ie)
      v(ii        ,npoints+jj,kk) = ghost%buf(ii,jj,kk,kptr,in)
      v(1-jj      ,ii        ,kk) = ghost%buf(ii,jj,kk,kptr,iw)
    end do
    end do
    end do

    ! SWEST
    !-----------
    do ll=swest,swest+max_corner_elem-1
      ic = desc%getmapP_ghost(ll)       
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,1-jj,kk)=ghost%buf(jj,ii,kk,kptr,ic)
          end do
          end do
          end do
        else
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,1-jj,kk)=ghost%buf(ii,jj,kk,kptr,ic)
          end do
          end do
          end do
        endif
      else
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(1-ii,1-jj,kk)=NaN
        end do
        end do
        end do
      endif
    end do
    
    ! SEAST
    !-----------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,1-jj,kk)=ghost%buf(jj,ii,kk,kptr,ic)
          end do
          end do
          end do
        else
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,1-jj,kk)=ghost%buf(ii,jj,kk,kptr,ic)
          end do
          end do
          end do
        endif
      else
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(npoints+ii,1-jj,kk)=NaN
        end do
        end do
        end do
      endif
    end do

    ! NEAST
    !-----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,npoints+jj,kk)=ghost%buf(jj,ii,kk,kptr,ic)
          end do
          end do
          end do
        else
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,npoints+jj,kk)=ghost%buf(ii,jj,kk,kptr,ic)
          end do
          end do
          end do
        endif
      else
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(npoints+ii,npoints+jj,kk)=NaN
        end do
        end do 
        end do
      endif
    end do

    ! NWEST
    !-----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(reverse) then
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,npoints+jj,kk)=ghost%buf(jj,ii,kk,kptr,ic)
          end do
          end do
          end do
        else
          do kk=1,vlyr
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,npoints+jj,kk)=ghost%buf(ii,jj,kk,kptr,ic)
          end do
          end do
          end do
        endif
      else
        do kk=1,vlyr
        do jj=1,nhc
        do ii=1,nhc
          v(1-ii,npoints+jj,kk)=NaN
        end do
        end do
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine ghostVunpack2d_level
  !==================================================================


  !==================================================================
  subroutine ghostVpack2d_single(ghost,v,nhc, npoints,desc)
    ! ghostVpack2d_single:
    !
    ! GHOSTVPACK2D
    ! AUTHOR: Christoph Erath 
    ! Pack edges of v into an ghost buffer for boundary exchange.
    !
    ! This subroutine packs for one vertical layers into an ghost
    ! buffer. It is for cartesian points (v is only two dimensional). 
    ! If the buffer associated with edge is not large enough to 
    ! hold all vertical layers you intent to pack, the method will 
    ! halt the program with a call to parallel_mod::endrun().
    ! INPUT: 
    ! - ghost Buffer into which the data will be packed.
    !   This buffer must be previously allocated with initGhostBuffer().
    ! - v The data to be packed.
    ! - nhc deep of ghost/halo zone
    ! - npoints number of points on on side
    ! - kptr Vertical pointer to the place in the edge buffer where 
    ! data will be located.
    !=================================================================================
    use SE_Constants,only: ntrac_d
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !---------------------
    type(Ghostbuffertr_t)            :: ghost
    real(kind=real_kind)  ,intent(in):: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc)
    integer               ,intent(in):: nhc,npoints
    type(EdgeDescriptor_t),intent(in):: desc
    !
    ! Local Values
    !----------------
    integer :: itr,k
    integer :: is,ie,in,iw
    integer :: ii,jj,ir,ll,ee


    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
  !$OMP BARRIER
#endif
      threadsafe=.true.
    endif

    ! Example convenction for buffer to the north:
    !  buf(:,,:,i,e) 
    !   each "edge" is a row of data (i=1,np) in the element 
    !     north most row of data goes into e=1
    !     next row of data goes into e=2
    !     .... 
    !     south most row of data goes into e=np
    !   We need to pack this way to preserve the orientation
    !   so the data can be unpacked correctly
    !
    ! note: we think of buf as dimensioned buf(k,is,i,e)
    ! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
    !---------------------------------------------------------------
    is = desc%putmapP_ghost(south) 
    ie = desc%putmapP_ghost(east )  
    in = desc%putmapP_ghost(north) 
    iw = desc%putmapP_ghost(west )  

    do ii=1,npoints
    do jj=1,nhc
      ghost%buf(ii,jj,1,1,is)   = v( ii       , jj+1     )
      ghost%buf(ii,jj,1,1,ie)   = v(npoints-jj, ii       )
      ghost%buf(ii,jj,1,1,in)   = v(        ii,npoints-jj)
      ghost%buf(ii,jj,1,1,iw)   = v( jj+1     , ii       )
    end do
    end do

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !  Check if the edge orientation of the recieving element is different
    !  if it is, swap the order of data in the edge
    !-------------------------------------------------------
    if(desc%reverse(south)) then
!     is = desc%putmapP_ghost(south)
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,1,1,is)=v(ii,jj+1)
      end do
      end do
    endif

    if(desc%reverse(east)) then
!     ie = desc%putmapP_ghost(east) 
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,1,1,ie)=v(npoints-jj,ii)
      end do
      end do
    endif

    if(desc%reverse(north)) then
!     in = desc%putmapP_ghost(north)
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,1,1,in)=v(ii,npoints-jj)
      end do
      end do
    endif

    if(desc%reverse(west)) then
!     iw = desc%putmapP_ghost(west)
      do ii=1,npoints
      do jj=1,nhc
        ir = npoints-ii+1
        ghost%buf(ir,jj,1,1,iw)=v(jj+1,ii)
      end do
      end do
    endif

    ! corners.  this is difficult because we dont know the orientaton
    ! of the corners, and this which (i,j) dimension maps to which dimension
    !----------------------------------------------------------------------
 
    ! SWEST
    !--------
    do ll=swest,swest+max_corner_elem-1
      if(ll.ne.swest) call endrun('ERROR3: swest ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do ii=1,nhc
        do jj=1,nhc
          ghost%buf(ii,jj,1,1,desc%putmapP_ghost(ll))=v(ii+1,jj+1)
        end do
        end do           
      endif
    end do

    ! SEAST
    !--------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(ll.ne.seast) call endrun('ERROR3: seast ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do ii=1,nhc             
        do jj=1,nhc
          ghost%buf(ii,jj,1,1,desc%putmapP_ghost(ll))=v(npoints-ii,jj+1)
        end do
        end do         
      endif
    end do

    ! NEAST
    !--------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(ll.ne.neast) call endrun('ERROR3: neast ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do ii=1,nhc
        do jj=1,nhc
          ghost%buf(ii,jj,1,1,desc%putmapP_ghost(ll))=v(npoints-ii,npoints-jj)           
        end do
        end do
      endif
    end do

    ! NWEST
    !--------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(ll.ne.nwest) call endrun('ERROR3: nwest ghost cell update requires ne>0 cubed-sphere mesh')
      if(desc%putmapP_ghost(ll) /= -1) then
        do ii=1,nhc
        do jj=1,nhc
          ghost%buf(ii,jj,1,1,desc%putmapP_ghost(ll))=v(ii+1,npoints-jj)
        end do       
        end do     
      endif
    end do   

    ! End Routine 
    !-------------
    return
  end subroutine ghostVpack2d_single
  !==================================================================


  !==================================================================
  subroutine ghostVunpack2d_single(ghost,v,nhc,npoints,desc)
    ! ghostVunpack2d_single:
    !
    ! GHOSTVUNPACK2D
    ! AUTHOR: Christoph Erath 
    ! Unpack ghost points from ghost buffer into v...
    ! It is for cartesian points (v is only two dimensional).
    ! INPUT SAME arguments as for GHOSTVPACK2d
    !=======================================================================
    use SE_Constants,only: ntrac_d
    use SE_Options  ,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !--------------------
    type(Ghostbuffertr_t) ,intent(in   ):: ghost
    real(kind=real_kind)  ,intent(inout):: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc)
    integer               ,intent(in   ):: nhc,npoints
    type(EdgeDescriptor_t)              :: desc
    !
    ! Local Values
    !----------------
    logical,parameter   :: UseUnroll = .TRUE.
    real(kind=real_kind):: NaN=-1.0
    integer :: is,ie,in,iw,ic
    logical :: reverse
    integer :: ii,jj,ll, itr, kk
    
    NaN=sqrt(NaN)

    threadsafe=.false.

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east )  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west )  

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1)
    ! 2nd   row ('edge') goes in v(:,np+2)
    ! etc...
    !-----------------------------------------
    do ii=1,npoints
    do jj=1,nhc
      v(ii        ,1-jj      ) = ghost%buf(ii,jj,1,1,is)
      v(npoints+jj,ii        ) = ghost%buf(ii,jj,1,1,ie)
      v(ii        ,npoints+jj) = ghost%buf(ii,jj,1,1,in)
      v(1-jj      ,ii        ) = ghost%buf(ii,jj,1,1,iw)
    end do
    end do

    ! SWEST
    !-----------
    do ll=swest,swest+max_corner_elem-1
      ic = desc%getmapP_ghost(ll)       
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if (reverse) then
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,1-jj)=ghost%buf(jj,ii,1,1,ic)
          end do
          end do
        else
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,1-jj)=ghost%buf(ii,jj,1,1,ic)
          end do
          end do
        endif
      else
        do jj=1,nhc
        do ii=1,nhc
          v(1-ii,1-jj)=NaN
        end do
        end do     
      endif
    end do
    
    ! SEAST
    !----------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if (reverse) then
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,1-jj)=ghost%buf(jj,ii,1,1,ic)
          end do
          end do
        else
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,1-jj)=ghost%buf(ii,jj,1,1,ic)
          end do
          end do
        endif
      else
        do jj=1,nhc
        do ii=1,nhc
          v(npoints+ii,1-jj)=NaN
        end do
        end do
      endif
    end do

    ! NEAST
    !----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if (reverse) then
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,npoints+jj)=ghost%buf(jj,ii,1,1,ic)
          end do
          end do
        else
          do jj=1,nhc
          do ii=1,nhc
            v(npoints+ii,npoints+jj)=ghost%buf(ii,jj,1,1,ic)
          end do
          end do
        endif
      else
        do jj=1,nhc
        do ii=1,nhc
          v(npoints+ii,npoints+jj)=NaN
        end do
        end do  
      endif
    end do

    ! NWEST
    !----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if (reverse) then
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,npoints+jj)=ghost%buf(jj,ii,1,1,ic)
          end do
          end do
        else
          do jj=1,nhc
          do ii=1,nhc
            v(1-ii,npoints+jj)=ghost%buf(ii,jj,1,1,ic)
          end do
          end do
        endif
      else
        do jj=1,nhc
        do ii=1,nhc
          v(1-ii,npoints+jj)=NaN
        end do
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine ghostVunpack2d_single
  !==================================================================


  !==================================================================
  subroutine initGhostBuffer3d(ghost,nlyr,np,nhc_in)
    ! initGhostBuffer3d:
    ! Author: James Overfelt
    ! create an Real based communication buffer
    ! npoints is the number of points on one side
    ! nhc is the deep of the ghost/halo zone
    !===============================================================
    !
    ! Passed Variables
    !---------------------
    type(Ghostbuffer3d_t),intent(out):: ghost
    integer              ,intent(in ):: nlyr, np
    integer,optional     ,intent(in ):: nhc_in
    !
    ! Local Values
    !------------------
    integer:: nbuf,nhc,ii

    ! sanity check for threading
    !----------------------------
    if(omp_get_num_threads()>1) then
      call endrun('ERROR: initGhostBuffer must be called before threaded region')
    endif

    if(present(nhc_in)) then
      nhc=nhc_in
    else
      nhc = np-1
    endif

    nbuf=max_neigh_edges*nelemd
#if (! defined ELEMENT_OPENMP)
!$OMP MASTER
#endif
    ghost%nlyr      = nlyr
    ghost%nhc       = nhc
    ghost%np        = np
    ghost%nbuf      = nbuf
    ghost%elem_size = np*(nhc+1)

    allocate(ghost%buf    (np,(nhc+1),nlyr,nbuf))
    allocate(ghost%receive(np,(nhc+1),nlyr,nbuf))
    ghost%buf       =0
    ghost%receive   =0
#if (! defined ELEMENT_OPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif

    ! End Routine 
    !-------------
    return
  end subroutine initGhostBuffer3d
  !==================================================================


  !==================================================================
  subroutine ghostVpack3d(ghost, v, vlyr, kptr, desc)
    ! ghostVpack3d:
    !
    ! GHOSTVPACK3D
    ! AUTHOR: James Overfelt (from a subroutine of Christoph Erath, ghostvpack2D)
    ! Pack edges of v into an ghost buffer for boundary exchange.
    !
    ! This subroutine packs for many vertical layers into an ghost
    ! buffer. 
    ! If the buffer associated with edge is not large enough to 
    ! hold all vertical layers you intent to pack, the method will 
    ! halt the program with a call to parallel_mod::endrun().
    ! INPUT: 
    ! - ghost Buffer into which the data will be packed.
    !   This buffer must be previously allocated with initGhostBuffer().
    ! - v The data to be packed.
    ! - nhc deep of ghost/halo zone
    ! - npoints number of points on on side
    ! - kptr Vertical pointer to the place in the edge buffer where 
    ! data will be located.
    !=================================================================================
    use SE_Options,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !---------------------
    integer               ,intent(in):: vlyr
    type(Ghostbuffer3d_t)            :: ghost
    real(kind=real_kind)  ,intent(in):: v(ghost%np, ghost%np, vlyr)
    integer               ,intent(in):: kptr
    type(EdgeDescriptor_t),intent(in):: desc
    !
    ! Local Values
    !---------------
    integer:: nhc, np
    integer:: is,ie,in,iw
    integer:: ii,jj,kk,ir,ll,ee

    if(.not. threadsafe) then
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
      threadsafe=.true.
    endif

    ! Example convenction for buffer to the north:
    !  buf(:,,:,i,e) 
    !   each "edge" is a row of data (i=1,np) in the element 
    !     north most row of data goes into e=1
    !     next row of data goes into e=2
    !     .... 
    !     south most row of data goes into e=np
    !   We need to pack this way to preserve the orientation
    !   so the data can be unpacked correctly
    !
    ! note: we think of buf as dimensioned buf(k,is,i,e)
    ! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
    !-------------------------------------------------------------
    nhc= ghost%nhc
    np = ghost%np
    is = desc%putmapP_ghost(south) 
    ie = desc%putmapP_ghost(east )  
    in = desc%putmapP_ghost(north) 
    iw = desc%putmapP_ghost(west )  

    do kk=1,vlyr
    do jj=1,nhc
    do ii=1,np
      ghost%buf(ii,jj,kptr+kk,is)   = v(ii   ,jj+1 ,kk)
      ghost%buf(ii,jj,kptr+kk,ie)   = v(np-jj,ii   ,kk)
      ghost%buf(ii,jj,kptr+kk,in)   = v(ii   ,np-jj,kk)
      ghost%buf(ii,jj,kptr+kk,iw)   = v(jj+1 ,ii   ,kk)
    end do
    end do
    end do

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !  Check if the edge orientation of the recieving element is different
    !  if it is, swap the order of data in the edge
    !-----------------------------------------------------------------------
    if(desc%reverse(south)) then
      do kk=1,vlyr
      do jj=1,nhc
      do ii=1,np
        ir = np-ii+1
        ghost%buf(ir, jj, kptr+kk, is)=v(ii, jj+1, kk)
      end do
      end do
      end do
    endif

    if(desc%reverse(east)) then
      do kk=1,vlyr
      do jj=1,nhc
      do ii=1,np
        ir = np-ii+1
        ghost%buf(ir, jj, kptr+kk, ie)=v(np-jj, ii, kk)
      end do
      end do
      end do
    endif

    if(desc%reverse(north)) then
      do kk=1,vlyr
      do jj=1,nhc
      do ii=1,np
        ir = np-ii+1
        ghost%buf(ir, jj, kptr+kk, in)=v(ii, np-jj, kk)
      end do
      end do
      end do
    endif

    if(desc%reverse(west)) then
      do kk=1,vlyr
      do jj=1,nhc
      do ii=1,np
        ir = np-ii+1
        ghost%buf(ir, jj, kptr+kk, iw)=v(jj+1, ii, kk)
      end do
      end do
      end do
    endif

    ! corners.  this is difficult because we dont know the orientaton
    ! of the corners, and this which (i,j) dimension maps to which dimension
    !----------------------------------------------------------------------

    ! SWEST
    !----------
    do ll=swest, swest+max_corner_elem-1
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do jj=1,nhc+1
        do ii=1,nhc+1
          ghost%buf(ii,jj,kptr+kk,desc%putmapP_ghost(ll))=v(ii,jj,kk)
        end do
        end do   
        end do
      endif
    end do

    ! SEAST
    !----------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do jj=1,nhc+1
        do ii=1,nhc+1           
          ghost%buf(ii,jj,kptr+kk,desc%putmapP_ghost(ll))=v(np-ii+1,jj,kk)
        end do
        end do
        end do
      endif
    end do
    
    ! NEAST
    !----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do jj=1,nhc+1
        do ii=1,nhc+1
          ghost%buf(ii,jj,kptr+kk,desc%putmapP_ghost(ll))=v(np-ii+1,np-jj+1,kk)           
        end do
        end do
        end do
      endif
    end do
    
    ! NWEST
    !----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      if(desc%putmapP_ghost(ll) /= -1) then
        do kk=1,vlyr
        do jj=1,nhc+1
        do ii=1,nhc+1
          ghost%buf(ii,jj,kptr+kk,desc%putmapP_ghost(ll))=v(ii,np-jj+1,kk)
        end do
        end do
        end do
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine ghostVpack3d
  !==================================================================
    

  !==================================================================
  subroutine ghostVunpack3d(g, v, vlyr, kptr, desc, sw, se, nw, ne, mult)
    ! ghostVunpack3d:
    !
    ! GHOSTVUNPACK3D
    ! AUTHOR: James Overfelt (from a subroutine of Christoph Erath, ghostVunpack2d)
    ! Unpack ghost points from ghost buffer into v...
    ! It is for cartesian points (v is only two dimensional).
    ! INPUT SAME arguments as for GHOSTVPACK
    !=============================================================================
    use SE_Options,only: north,south,east,west,neast,nwest,seast,swest
    !
    ! Passed Variables
    !-------------------
    integer              ,intent(in   ):: vlyr
    type(Ghostbuffer3d_t),intent(in   ):: g
    real(kind=real_kind) ,intent(inout):: v (1-g%nhc:g%np+g%nhc, 1-g%nhc:g%np+g%nhc, vlyr)
    integer              ,intent(in   ):: kptr
    type (EdgeDescriptor_t)            :: desc
    real(kind=real_kind) ,intent(  out):: sw(1-g%nhc:1         ,1-g%nhc:1         ,vlyr,max_corner_elem-1)
    real(kind=real_kind) ,intent(  out):: se(   g%np:g%np+g%nhc,1-g%nhc:1         ,vlyr,max_corner_elem-1)
    real(kind=real_kind) ,intent(  out):: nw(1-g%nhc:1         ,   g%np:g%np+g%nhc,vlyr,max_corner_elem-1)
    real(kind=real_kind) ,intent(  out):: ne(   g%np:g%np+g%nhc,   g%np:g%np+g%nhc,vlyr,max_corner_elem-1)
    integer              ,intent(  out):: mult(5:8)
    !
    ! Local Values
    !---------------
    logical,parameter   :: UseUnroll = .TRUE.
    real(kind=real_kind):: NaN=-1.0
    integer:: nhc, np
    integer:: is,ie,in,iw,ic
    logical:: reverse
    integer:: ii,jj,kk,ll
    
    NaN=sqrt(NaN)

    threadsafe=.false.

    nhc = g%nhc
    np  = g%np

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east )  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west )  

    ! fill in optional values with NaN
    !-----------------------------------
    do kk=1,vlyr
    do jj=1,nhc
    do ii=1,nhc
      v( 1-ii, 1-jj,kk)=NaN
      v(np+ii, 1-jj,kk)=NaN
      v(np+ii,np+jj,kk)=NaN
      v( 1-ii,np+jj,kk)=NaN
    end do
    end do
    end do

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1)
    ! 2nd   row ('edge') goes in v(:,np+2)
    ! etc...
    !--------------------------------------
    do kk=1,vlyr
    do jj=1,nhc
    do ii=1,np
      v(   ii, 1-jj,kk) = g%buf(ii,jj,kptr+kk,is)
      v(np+jj,   ii,kk) = g%buf(ii,jj,kptr+kk,ie)
      v(   ii,np+jj,kk) = g%buf(ii,jj,kptr+kk,in)
      v( 1-jj,   ii,kk) = g%buf(ii,jj,kptr+kk,iw)
    end do
    end do
    end do

    ! four sides are always just one
    !----------------------------------
    mult(swest) = 0
    mult(seast) = 0
    mult(neast) = 0
    mult(nwest) = 0

    ! SWEST
    !-----------
    do ll=swest, swest+max_corner_elem-1
      ic = desc%getmapP_ghost(ll)       
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(mult(swest) .eq. 0) then
          if(reverse) then
            do kk=1,vlyr
            do jj=1,nhc
            do ii=1,nhc
              v(1-ii,1-jj,kk)=g%buf(jj+1,ii+1,kptr+kk,ic)
            end do
            end do
            end do
          else
            do kk=1,vlyr
            do jj=1,nhc
            do ii=1,nhc
              v(1-ii,1-jj,kk)=g%buf(ii+1,jj+1,kptr+kk,ic)
            end do
            end do
            end do
          endif
        else
          if(reverse) then
            do kk=1,vlyr
            do jj=0,nhc
            do ii=0,nhc
              sw(1-ii,1-jj,kk,mult(swest))=g%buf(jj+1,ii+1,kptr+kk,ic)
            end do
            end do
            end do
          else
            do kk=1,vlyr
            do jj=0,nhc
            do ii=0,nhc
              sw(1-ii,1-jj,kk,mult(swest))=g%buf(ii+1,jj+1,kptr+kk,ic)
            end do
            end do
            end do
          endif
        endif
        mult(swest) = mult(swest) + 1
      endif
    end do
    
    ! SEAST
    !-----------
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(mult(seast) .eq. 0) then
          if(reverse) then
            do kk=1,vlyr
            do jj=1,nhc
            do ii=1,nhc
              v(np+ii,1-jj,kk)=g%buf(jj+1,ii+1,kptr+kk,ic)
            end do
            end do
            end do
          else
            do kk=1,vlyr
            do jj=1,nhc
            do ii=1,nhc
              v(np+ii,1-jj,kk)=g%buf(ii+1,jj+1,kptr+kk,ic)
            end do
            end do
            end do
          endif
        else
          if(reverse) then
            do kk=1,vlyr
            do jj=0,nhc
            do ii=0,nhc
              se(np+ii,1-jj,kk,mult(seast))=g%buf(jj+1,ii+1,kptr+kk,ic)
            end do
            end do
            end do
          else
            do kk=1,vlyr
            do jj=0,nhc
            do ii=0,nhc
              se(np+ii,1-jj,kk,mult(seast))=g%buf(ii+1,jj+1,kptr+kk,ic)
            end do
            end do
            end do
          endif
        endif
        mult(seast) = mult(seast) + 1
      endif
    end do

    ! NEAST
    !-----------
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(mult(neast) .eq. 0) then
          if(reverse) then
            do kk=1,vlyr
            do jj=1,nhc
            do ii=1,nhc
              v(np+ii,np+jj,kk)=g%buf(jj+1,ii+1,kptr+kk,ic)
            end do
            end do
            end do
          else
            do kk=1,vlyr
            do jj=1,nhc
            do ii=1,nhc
              v(np+ii,np+jj,kk)=g%buf(ii+1,jj+1,kptr+kk,ic)
            end do
            end do
            end do
          endif
        else
          if(reverse) then
            do kk=1,vlyr
            do jj=0,nhc
            do ii=0,nhc
              ne(np+ii,np+jj,kk,mult(neast))=g%buf(jj+1,ii+1,kptr+kk,ic)
            end do
            end do
            end do
          else
            do kk=1,vlyr
            do jj=0,nhc
            do ii=0,nhc
              ne(np+ii,np+jj,kk,mult(neast))=g%buf(ii+1,jj+1,kptr+kk,ic)
            end do
            end do
            end do
          endif
        endif
        mult(neast) = mult(neast) + 1
      endif
    end do

    ! NWEST
    !-----------
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
      ic = desc%getmapP_ghost(ll)
      if(ic /= -1) then 
        reverse=desc%reverse(ll)
        if(mult(nwest) .eq. 0) then
          if(reverse) then
            do kk=1,vlyr
            do jj=1,nhc
            do ii=1,nhc
              v(1-ii,np+jj,kk)=g%buf(jj+1,ii+1,kptr+kk,ic)
            end do
            end do
            end do
          else
            do kk=1,vlyr
            do jj=1,nhc
            do ii=1,nhc
              v(1-ii,np+jj,kk)=g%buf(ii+1,jj+1,kptr+kk,ic)
            end do
            end do
            end do
          endif
        else
          if(reverse) then
            do kk=1,vlyr
            do jj=0,nhc
            do ii=0,nhc
              nw(1-ii,np+jj,kk,mult(nwest))=g%buf(jj+1,ii+1,kptr+kk,ic)
            end do
            end do
            end do
          else
            do kk=1,vlyr
            do jj=0,nhc
            do ii=0,nhc
              nw(1-ii,np+jj,kk,mult(nwest))=g%buf(ii+1,jj+1,kptr+kk,ic)
            end do
            end do
            end do
          endif
        endif
        mult(nwest) = mult(nwest) + 1
      endif
    end do

    ! End Routine 
    !-------------
    return
  end subroutine ghostVunpack3d
  !==================================================================

end module edge_mod


#ifndef HAVE_F2003_PTR_BND_REMAP
!==================================================================
subroutine remap_2D_ptr_buf(edge,nlyr,nbuf,src_array)
  ! remap_2D_ptr_buf:
  !
  ! subroutine to allow sharing edge buffers
  ! this has to be outside a module to allow us to (F77 style) access the same chunk 
  ! of memory with a different shape
  !
  ! some compilers dont allow the 'target' attribute to be used in a F77 subroutine
  ! such as cray.  if that is the case, try compiling with -DHAVE_F2003_PTR_BND_REMAP
  !================================================================================
  use edge_mod    ,only: EdgeBuffer_t ! _EXTERNAL    
  use SE_Constants,only: real_kind
  !
  ! Passed Variables
  !-------------------
  type(EdgeBuffer_t)         :: edge
  integer                    :: nlyr,nbuf
  real(kind=real_kind),target:: src_array(nlyr,nbuf)

  edge%buf  => src_array

  ! End Routine
  !----------------
  return
end subroutine remap_2D_ptr_buf
!==================================================================


!==================================================================
subroutine remap_2D_ptr_receive(edge,nlyr,nbuf,src_array)
  ! remap_2D_ptr_receive:
  ! 
  !=============================================================
  use edge_mod    ,only: EdgeBuffer_t ! _EXTERNAL   
  use SE_Constants,only: real_kind
  !
  ! Passed Variables
  !-------------------
  type(EdgeBuffer_t)         :: edge
  integer                    :: nlyr,nbuf
  real(kind=real_kind),target:: src_array(nlyr,nbuf)

  edge%receive  => src_array

  ! End Routine
  !----------------
  return
end subroutine remap_2D_ptr_receive
!==================================================================
#endif

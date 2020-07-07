#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module metagraph_mod

  ! Usefule Modules
  !---------------------
  use SE_Constants ,only: int_kind
  use err_exit     ,only: iulog
  use gridgraph_mod,only: gridvertex_t,gridedge_t,allocate_gridvertex_nbrs,assignment ( = )
!PFC  use pio_types    ! _EXTERNAL

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none 
  private 

  public :: MetaEdge_t
  public :: MetaVertex_t

  private:: copy_metaedge
! public :: findedge
! public :: MetaEdgeCount
  public :: LocalElemCount
  public :: edge_uses_vertex
  public :: PrintMetaEdge
  public :: PrintMetaVertex
  public :: initMetaGraph

  ! Publice Interfaces
  !----------------------
  interface assignment ( = )
    module procedure copy_metaedge
  end interface

  ! Type Definitions
  !------------------
  type MetaEdge_t
    type(GridEdge_t),pointer:: members(:)
    integer         ,pointer:: edgeptrP(:)
    integer         ,pointer:: edgeptrP_ghost(:)
    integer                 :: number
    integer                 :: type 
    integer                 :: wgtP       ! sum of lengths of all messages to pack for edges
    integer                 :: wgtP_ghost ! sum of lengths of all messages to pack for ghost cells
    integer                 :: HeadVertex ! processor number to send to
    integer                 :: TailVertex ! processor number to send from
    integer                 :: nmembers   ! number of messages to (un)pack (out)into this buffer
    integer                 :: padding    ! just to quite compiler
  end type MetaEdge_t

  type MetaVertex_t                         ! one for each processor
    integer                   :: number     ! USELESS just the local processor number
    integer                   :: nmembers   ! number of elements on this processor 
    type(GridVertex_t),pointer:: members(:) ! array of elements on this processor
    type(MetaEdge_t)  ,pointer:: edges(:)   ! description of messages to send/receive
    integer                   :: nedges     ! number of processors to communicate with (length of edges)
    integer                   :: padding    ! just to quite compiler
  end type MetaVertex_t


contains 

  ! =====================================
  ! copy vertex:
  ! copy device for overloading = sign.
  ! =====================================

  !===========================================================
  recursive subroutine copy_metaedge(edge2,edge1)
    ! copy_metaedge:
    !
    !========================================================
    !
    ! Passed Variables
    !------------------
    type(MetaEdge_t),intent(out):: edge2
    type(MetaEdge_t),intent(in ):: edge1
    !
    ! Local Values
    !--------------
    integer ii

    edge2%number    = edge1%number
    edge2%type      = edge1%type
    edge2%wgtP      = edge1%wgtP
    edge2%wgtP_ghost= edge1%wgtP_ghost
    edge2%nmembers  = edge1%nmembers

    if(associated(edge1%members)) then
      allocate(edge2%members(edge2%nmembers))
      do ii=1,edge2%nmembers
        edge2%members(ii)=edge1%members(ii)
      end do
    endif

    if(associated(edge1%edgeptrP)) then
      allocate(edge2%edgeptrP(edge2%nmembers))
      allocate(edge2%edgeptrP_ghost(edge2%nmembers))
      do ii=1,edge2%nmembers
        edge2%edgeptrP(ii)      = edge1%edgeptrP(ii)
        edge2%edgeptrP_ghost(ii)= edge1%edgeptrP_ghost(ii)
      end do
    endif

    edge2%HeadVertex = edge1%HeadVertex
    edge2%TailVertex = edge1%TailVertex

    ! End Routine
    !-------------
    return
  end subroutine copy_metaedge
  !===========================================================


  !===========================================================
! function findedge(mEdge,Edge) result(number)
!
!   type(MetaEdge_t), intent(inout) :: mEdge(:)
!   type(GridEdge_t), intent(in)    :: Edge
!   integer :: number
!
!   integer :: head,tail,exist
!
!   integer :: nedge
!   integer :: i
!
!   nedge=SIZE(mEdge)
!   number = 0
!   tail=Edge%tail%processor_number
!   head=Edge%head%processor_number
!
!   exist=0
!   do i=1,nedge
!      !       write(iulog,*)'mEdge(i)%number: ',mEdge(i)%number
!      if(mEdge(i)%number .ne. 0) then 
!         if  ((mEdge(i)%TailVertex==tail .and. mEdge(i)%HeadVertex==head) ) then
!            number=i
!         end if
!         exist=exist+1
!      endif
!   end do
!   if(number == 0) number = exist + 1
!
! end function findedge
  !===========================================================


  !===========================================================
! function MetaEdgeCount(Edge) result(nMedge)
!   implicit none
!
!   type (GridEdge_t),intent(in)  :: Edge(:)
!   integer                       :: nMedge
!
!   integer                       :: nedge,i,j,maxedges
!   integer                       :: head_processor_number,tail_processor_number
!   integer, allocatable          :: tmp(:,:)
!   logical                       :: found
!
!
!   nedge = SIZE(Edge)
!   maxedges=nedge
!
!   allocate(tmp(2,maxedges))
!   tmp=0
!   nMedge=0
!   do i=1,nedge
!      head_processor_number = Edge(i)%head%processor_number
!      tail_processor_number = Edge(i)%tail%processor_number
!      found = .FALSE.
!      do j=1,nMedge
!         if((tmp(1,j) .eq. head_processor_number).and. &
!              (tmp(2,j) .eq. tail_processor_number)) found=.TRUE.
!      enddo
!      if(.NOT. found) then 
!         nMedge=nMedge+1
!         tmp(1,nMedge) = head_processor_number
!         tmp(2,nMedge) = tail_processor_number
!      endif
!   enddo
!   !mem    write(iulog,*)'MetaEdgeCount: before call to deallocate(tmp)'
!   deallocate(tmp)
!
! end function MetaEdgeCount
  !===========================================================


  !===========================================================
  function LocalElemCount(Vertex) result(nelemd)
    ! LocalElemCount:
    !
    !=========================================================
    !
    ! Passed Variables
    !---------------------
    type(MetaVertex_t),intent(in):: Vertex
    integer                      :: nelemd

    nelemd=Vertex%nmembers

    ! End Function
    !--------------
    return
  end function LocalElemCount
  !===========================================================


  !===========================================================
  function edge_uses_vertex(Vertex,Edge) result(log)
    ! edge_uses_vertex:
    !
    !========================================================
    !
    ! Passed Variables
    !----------------------
    type(MetaVertex_t),intent(in):: Vertex
    type(MetaEdge_t)  ,intent(in):: Edge
    logical                      :: log
    !
    ! Local Values
    !-----------------
    integer:: number

    number = Vertex%number
    if((number == Edge%HeadVertex).or. &
       (number == Edge%TailVertex)     ) then 
      log = .TRUE.
    else
      log = .FALSE. 
    endif

    ! End Function
    !---------------
    return
  end function edge_uses_vertex
  !===========================================================


  !===========================================================
  subroutine PrintMetaEdge(Edge)
    ! PrintMetaEdge:
    !
    !=======================================================
    use gridgraph_mod,only: PrintGridEdge
    !
    ! Passed Variables
    !--------------------
    type(MetaEdge_t),intent(in):: Edge(:)
    !
    ! Local Values
    !--------------
    integer:: ii,nedge

    nedge = SIZE(Edge)
    do ii=1,nedge
      print *
      write(iulog,90)  Edge(ii)%number,Edge(ii)%type,Edge(ii)%wgtP,Edge(ii)%nmembers, &
                       Edge(ii)%TailVertex, Edge(ii)%HeadVertex
      if(associated(Edge(ii)%members)) then
        call PrintGridEdge(Edge(ii)%members)
      endif
    end do

90  format('METAEDGE #',I4,2x,'TYPE ',I1,2x,'WGT ',I4,2x,'NUM ',I6,2x,'Processors ',I4,' ---> ',I4)

    ! End Routine
    !--------------
    return
  end subroutine PrintMetaEdge
  !===========================================================


  !===========================================================
  subroutine PrintMetaVertex(Vertex)
    ! PrintMetaVertex:
    !
    !=========================================================
    use gridgraph_mod,only: PrintGridVertex
    !
    ! Passed Variables
    !---------------------
    type(MetaVertex_t),intent(in),target:: Vertex
    !
    ! Local Values
    !---------------
    integer:: ii,jj,npart

    ii=0  ! dummy that wrote this has ii undefined before use, so set it to 0 for now.

    write(iulog,*)
    write(iulog,90) ii
    write(iulog,95) Vertex%nmembers
    call PrintGridVertex(Vertex%members)
    write(iulog,96) Vertex%nedges
    if(associated(Vertex%edges)) then 
      do jj=1,Vertex%nedges
        write(iulog,97) Vertex%edges(jj)%number    ,Vertex%edges(jj)%type      , &
                        Vertex%edges(jj)%wgtP      ,Vertex%edges(jj)%HeadVertex, &
                        Vertex%edges(jj)%TailVertex
      end do
    endif

90  format('METAVERTEX #',I2,2x)
95  format(5x,I2,' Member Grid Vertices')
96  format(5x,I2,' Incident Meta Edges ')
97  format(10x,'METAEDGE #',I2,2x,'TYPE ',I1,2x,'WGT ',I4,2x,'Processors ',I2,' ---> ',I2)
98  format(5x,'GRIDVERTEX #',5x,'West',5x,'East',5x,'South',5x,'North')
99  format(10x,I2,7x,4(2x,I3,1x,'(',I1,')'))

    ! End Routine
    !------------
    return
  end subroutine PrintMetaVertex
  !===========================================================


  !===========================================================
  subroutine initMetaGraph(ThisProcessorNumber,MetaVertex,GridVertex,GridEdge)
    ! initMetaGraph:
    !
    !=============================================================================
!PFC    use ll_mod       ,only: root_t,LLSetEdgeCount,LLFree,LLInsertEdge,LLGetEdgeCount,LLFindEdge
    use gridgraph_mod,only: GridEdge_type, printGridVertex
    ! 
    ! Passed Variables
    !---------------------
    integer(kind=int_kind)   ,intent(in ):: ThisProcessorNumber
    type(MetaVertex_t)       ,intent(out):: MetaVertex
    type(GridVertex_t),target,intent(in ):: GridVertex(:)
    type(GridEdge_t)  ,target,intent(in ):: GridEdge(:)
    !
    ! Local Values
    !--------------
    !type (MetaEdge_t), allocatable :: MetaEdge(:)
    integer            :: nelem,nelem_edge, nedges  
    integer,allocatable:: icount(:)
    integer            :: ic,i0,j0,ii
    integer            :: npart
    integer            :: head_processor_number
    integer            :: tail_processor_number
    integer            :: nedge_active,enum
    logical            :: found
    integer            :: iTail, iHead, wgtP
!PFC    type(root_t)       :: mEdgeList ! root_t = C++ std::set<std::pair<int,int> >
    logical            :: Verbose = .FALSE.
    logical            :: Debug   = .FALSE.

    if(Debug) write(iulog,*)'initMetagraph: point #1'

    ! Set number of grid vertices and edges
    !-----------------------------------------
    nelem      = SIZE(GridVertex)
    nelem_edge = SIZE(GridEdge  )

!PFC    mEdgeList%number = ThisProcessorNumber
!PFC    NULLIFY(mEdgeList%first)
!PFC    call LLSetEdgeCount(0)

    do i0=1,nelem_edge
      tail_processor_number = GridEdge(i0)%tail%processor_number
      head_processor_number = GridEdge(i0)%head%processor_number
      if((tail_processor_number.eq.ThisProcessorNumber).or.  &
         (head_processor_number.eq.ThisProcessorNumber)      ) then 
!PFC        call LLInsertEdge(mEdgeList,tail_processor_number,head_processor_number,eNum)
      endif
    end do
!PFC    call LLGetEdgeCount(nedges)

    NULLIFY(MetaVertex%edges)
        
    allocate(MetaVertex%edges(nedges))

    ! Initalize the Meta Vertices to zero... probably should be done
    ! in a separate routine
    !----------------------------------------------------------------
    MetaVertex%nmembers= 0  
    MetaVertex%number  = 0
    MetaVertex%nedges  = 0
    if(Debug) write(iulog,*)'initMetagraph: point #2'

    ! Give some identity to the Meta_vertex
    !-------------------------------------------
    MetaVertex%number = ThisProcessorNumber
    if(Debug) write(iulog,*)'initMetagraph: point #3'

    !  Look through all the small_vertices and determine the number of
    !  member vertices
    !-------------------------------------------------------------------
    if(Debug) call PrintGridVertex(GridVertex)
    if(Debug) write(iulog,*)'initMetagraph: After call to PrintGridVertex point #3.1'
    if(Debug) write(iulog,*)'initMetaGraph: ThisProcessorNumber is ',ThisProcessorNumber

    ! count number of elements on this processor
    !-------------------------------------------
    do j0=1,nelem  
      if(GridVertex(j0)%processor_number.eq.ThisProcessorNumber) then
        MetaVertex%nmembers = MetaVertex%nmembers + 1
      endif
    end do

    !  Allocate space for the members of the MetaVertices
    !  *sloppy* programmer again uses undefined values...... ARG!!!
    !--------------------------------------------------------------
    if(Debug) write(iulog,*)'initMetagraph: point #4 '
!    if(Debug) write(iulog,*)'initMetagraph: point #4.1 i,MetaVertex%nmembers',i,MetaVertex%nmembers
    allocate(MetaVertex%members(MetaVertex%nmembers))

    do j0=1, MetaVertex%nmembers
      call allocate_gridvertex_nbrs(MetaVertex%members(j0))
    end do
    if(Debug) write(iulog,*)'initMetagraph: point #5'

    ! Set the identity of the members of the MetaVertices
    !-------------------------------------------------------
    ic=1
    do j0=1,nelem
      if(GridVertex(j0)%processor_number.eq.ThisProcessorNumber) then 
        MetaVertex%members(ic) = GridVertex(j0)
        ic=ic+1
      endif
    end do

    nedges = SIZE(MetaVertex%edges)
    if(Debug) write(iulog,*)'initMetagraph: point #6 nedges',nedges

    !  Zero out all the edge numbers ... this should probably be
    !  move to some initalization routine
    !---------------------------------------------------------------
    MetaVertex%edges%number    = 0
    MetaVertex%edges%nmembers  = 0
    MetaVertex%edges%wgtP      = 0
    MetaVertex%edges%wgtP_ghost= 0
    do i0=1,nedges
      NULLIFY(MetaVertex%edges(i0)%members)
    end do
    if(Debug) write(iulog,*)'initMetagraph: point #7'

    ! Insert all the grid edges into the Meta Edges
    !-----------------------------------------------
    do i0=1, nelem_edge
      !  Which Meta Edge does this grid edge belong
      !-----------------------------------------------
      head_processor_number = GridEdge(i0)%head%processor_number
      tail_processor_number = GridEdge(i0)%tail%processor_number
!PFC      call LLFindEdge(mEdgeList,tail_processor_number,head_processor_number,j0,found)
      if(found) then 
        !  Increment the number of grid edges contained in the grid edge
        !  and setup the pointers
        !-----------------------------------------------------------------
        if(Debug) write(iulog,*)'initMetagraph: point #8'
        ii=GridEdge(i0)%tail_face
        wgtP                           = Gridedge(i0)%tail%nbrs_wgt(ii)
        MetaVertex%edges(j0)%nmembers  = MetaVertex%edges(j0)%nmembers+1
        MetaVertex%edges(j0)%wgtP      = MetaVertex%edges(j0)%wgtP + wgtP
        MetaVertex%edges(j0)%wgtP_ghost= MetaVertex%edges(j0)%wgtP_ghost    &
                                       + Gridedge(i0)%tail%nbrs_wgt_ghost(ii)
        if(Debug) write(iulog,*)'initMetagraph: point #9'

        !  If this the first grid edge to be inserted into the Meta Edge
        !  do some more stuff
        !---------------------------------------------------------------
        if(MetaVertex%edges(j0)%nmembers.eq.1) then
          if(Debug) write(iulog,*)'initMetagraph: point #10'
          MetaVertex%edges(j0)%number = j0                          ! its identity
          MetaVertex%edges(j0)%type   = gridedge_type(GridEdge(i0)) ! Type of grid edge
          if(Debug) write(iulog,*)'initMetagraph: point #11'

          !  Setup the pointer to the head and tail of the Vertex
          !--------------------------------------------------------
          MetaVertex%edges(j0)%HeadVertex = head_processor_number
          MetaVertex%edges(j0)%TailVertex = tail_processor_number
          if(Debug) write(iulog,*)'initMetagraph: point #12'

          !  Determine the number of edges for the Meta_Vertex
          !  This is the number of processors to communicate with
          !-------------------------------------------------------
          MetaVertex%nedges =  MetaVertex%nedges + 1
          if(Debug) write(iulog,*)'initMetagraph: point #13'
        endif
      endif
    end do

    !  Allocate space for the member edges and edge index
    !-----------------------------------------------------
    do i0=1,nedges
      allocate(MetaVertex%edges(i0)%members       (MetaVertex%edges(i0)%nmembers))
      allocate(MetaVertex%edges(i0)%edgeptrP      (MetaVertex%edges(i0)%nmembers))
      allocate(MetaVertex%edges(i0)%edgeptrP_ghost(MetaVertex%edges(i0)%nmembers))
      MetaVertex%edges(i0)%edgeptrP      (:) = 0
      MetaVertex%edges(i0)%edgeptrP_ghost(:) = 0
    end do
    if(Debug) write(iulog,*)'initMetagraph: point #14'

    !  Insert the edges into the proper meta edges
    !-----------------------------------------------
    allocate(icount(nelem_edge))
    icount=1
    do i0=1,nelem_edge
      head_processor_number = GridEdge(i0)%head%processor_number
      tail_processor_number = GridEdge(i0)%tail%processor_number
!PFC      call LLFindEdge(mEdgeList,tail_processor_number,head_processor_number,j0,found)
      if(found) then 
        MetaVertex%edges(j0)%members(icount(j0)) = GridEdge(i0)
        if((icount(j0)+1).le.MetaVertex%edges(j0)%nmembers) then
          ii  =GridEdge(i0)%tail_face
          wgtP=Gridedge(i0)%tail%nbrs_wgt(ii)
          MetaVertex%edges(j0)%edgeptrP(icount(j0)+1) = MetaVertex%edges(j0)%edgeptrP(icount(j0)) + wgtP
          wgtP=Gridedge(i0)%tail%nbrs_wgt_ghost(ii)
          MetaVertex%edges(j0)%edgeptrP_ghost(icount(j0)+1) = MetaVertex%edges(j0)%edgeptrP_ghost(icount(j0)) + wgtP
        endif
        if(Debug) write(iulog,*)'initMetagraph: point #15'
        icount(j0)=icount(j0)+1
      endif
    end do
    deallocate(icount)
    if(Debug) write(iulog,*)'initMetagraph: point #16'

    if(Verbose) then
      print *
      write(iulog,*)"edge bundle list:(INITMETAGRAPH)"
      call PrintMetaEdge( MetaVertex%edges)
      write(iulog,*)'initmetagrap: Before last call to PrintMetaVertex'
      call PrintMetaVertex(MetaVertex)
    endif

!PFC    call LLFree(mEdgeList)

90  format('EDGE #',I2,2x,'TYPE ',I1,2x,'Processor Numbers ',I2,' ---> ',I2)
100 format(10x,I2,1x,'(',I1,') ---> ',I2,1x,'(',I1,')')

    ! End Routine
    !---------------
    return
  end subroutine initMetaGraph
  !===========================================================

end module metagraph_mod

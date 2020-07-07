module GridGraph_mod

  ! Useful modules
  !----------------
  use SE_Constants,only: real_kind
  use SE_Options  ,only: north, south, east, west, neast, nwest, seast, swest
  use err_exit    ,only: iulog

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private

  public :: num_neighbors
  public :: GridVertex_t
  public :: EdgeIndex_t
  public :: GridEdge_t

  public :: assignment ( = ) 
  private:: copy_gridedge
  private:: copy_edgeindex
  private:: copy_gridvertex

  public :: allocate_gridvertex_nbrs
  public :: deallocate_gridvertex_nbrs
  public :: FreeGraph
  public :: gridedge_search
  public :: gridedge_type
  public :: grid_edge_uses_vertex
  public :: PrintChecksum
  public :: CreateSubGridGraph
  public :: PrintGridEdge
  public :: set_GridVertex_number
  public :: PrintGridVertex
  public :: CheckGridNeighbors
  public :: initgridedge

  ! PFC: Hack to set global option values
  !--------------------------------------
  public :: init_gridgraph_mod
  private:: max_neigh_edges
  private:: nelemd
  private:: max_corner_elem
  integer:: max_neigh_edges
  integer:: nelemd
  integer:: max_corner_elem

  ! Parameter Values
  !------------------
  integer,parameter:: num_neighbors=8 ! for north, south, east, west, neast, nwest, seast, swest

  ! Type Definitions
  !--------------------
  type GridVertex_t
    integer,pointer:: nbrs          (:) => null() ! The numbers of the neighbor elements
    integer,pointer:: nbrs_face     (:) => null() ! The cube face number of the neighbor element (nbrs array)
    integer,pointer:: nbrs_wgt      (:) => null() ! The weights for edges defined by nbrs array
    integer,pointer:: nbrs_wgt_ghost(:) => null() ! The weights for edges defined by nbrs array
    integer        :: nbrs_ptr(num_neighbors + 1) !index into the nbrs array for each neighbor direction
    integer        :: face_number                 ! which face of the cube this vertex is on
    integer        :: number                      ! element number
    integer        :: processor_number            ! processor number
    integer        :: SpaceCurve                  ! index in Space-Filling curve
  end type GridVertex_t

  type EdgeIndex_t
    integer,pointer:: ixP(:) => null()
    integer,pointer:: iyP(:) => null()
  end type EdgeIndex_t

  type GridEdge_t
    integer                   :: head_face      ! needed if head vertex has shape (i.e. square)
    integer                   :: tail_face      ! needed if tail vertex has shape (i.e. square)
    integer                   :: head_dir       ! which of 8 neighbor directions is the head
    integer                   :: tail_dir       ! which of 8 neighbor directions is the tail
    type(GridVertex_t),pointer:: head => null() ! edge head vertex
    type(GridVertex_t),pointer:: tail => null() ! edge tail vertex
    logical                   :: reverse
  end type GridEdge_t
  
  ! Public Interfaces
  !----------------------------
  interface assignment ( = )
    module procedure copy_gridedge
    module procedure copy_edgeindex
    module procedure copy_gridvertex
  end interface

contains
  !==================================================================
  subroutine init_gridgraph_mod(I_SEopt)
    ! PFC: Hack to set global option values
    !--------------------------------------
    use SE_Options,only: SEoptions_t
    ! Passed Variables
    !------------------
    type(SEoptions_t),intent(in):: I_SEopt

    max_neigh_edges = I_SEopt%max_neigh_edges
    max_corner_elem = I_SEopt%max_corner_elem
    nelemd          = I_SEopt%nelemd

    ! End Routine
    !-------------
    return
  end subroutine init_gridgraph_mod
  !==================================================================


  !==================================================================
  subroutine allocate_gridvertex_nbrs(vertex, dim)
    ! allocate_gridvertex_nbrs:
    !
    !================================================================
    !
    ! Passed Variables
    !-------------------
    type(GridVertex_t),intent(inout):: vertex
    integer,optional  ,intent(in   ):: dim
    !
    ! Local Values
    !-----------------
    integer num

    if(present(dim)) then
      num = dim
    else
      num = max_neigh_edges
    endif

    allocate(vertex%nbrs          (num))
    allocate(vertex%nbrs_face     (num))
    allocate(vertex%nbrs_wgt      (num))
    allocate(vertex%nbrs_wgt_ghost(num))
 
    ! End Routine
    !-------------
    return
  end subroutine allocate_gridvertex_nbrs
  !==================================================================


  !==================================================================
  subroutine deallocate_gridvertex_nbrs(vertex)
    ! deallocate_gridvertex_nbrs:
    !
    !=================================================================
    !
    ! Passed Variables
    !-------------------
    type(GridVertex_t),intent(inout):: vertex

    deallocate(vertex%nbrs          )
    deallocate(vertex%nbrs_face     )
    deallocate(vertex%nbrs_wgt      )
    deallocate(vertex%nbrs_wgt_ghost)
 
    ! End Routine
    !-------------
    return
  end subroutine deallocate_gridvertex_nbrs
  !==================================================================


  !==================================================================
  recursive subroutine copy_gridedge(edge2, edge1)
    ! copy_gridedge:
    !
    ! copy edge:
    ! copy device for overloading = sign.
    !==================================================
    ! 
    ! Passed Variables
    !------------------
    type(GridEdge_t),intent(out):: edge2
    type(GridEdge_t),intent(in ):: edge1

    edge2%tail_face = edge1%tail_face
    edge2%head_face = edge1%head_face
    edge2%tail_dir  = edge1%tail_dir
    edge2%head_dir  = edge1%head_dir
    edge2%reverse   = edge1%reverse

    if(associated(edge1%tail)) then
      edge2%tail=>edge1%tail
    endif

    if(associated(edge1%head)) then
      edge2%head=>edge1%head
    endif

    ! End Routine
    !-------------
    return
  end subroutine copy_gridedge
  !==================================================================

  !==================================================================
  recursive subroutine copy_gridvertex(vertex2, vertex1)
    ! copy_gridvertex:
    !
    !===============================================================
    !
    ! Passed Variables
    !------------------
    type(GridVertex_t),intent(out):: vertex2
    type(GridVertex_t),intent(in ):: vertex1
    !
    ! Local Values
    !----------------
    integer ii,jj,nn
   
    nn = SIZE(vertex1%nbrs)

    if(associated(vertex2%nbrs)) then
      nullify(vertex2%nbrs)
    endif
    if(associated(vertex2%nbrs_face)) then
      nullify(vertex2%nbrs_face)
    endif
    if(associated(vertex2%nbrs_wgt)) then
      nullify(vertex2%nbrs_wgt)
    endif
    if(associated(vertex2%nbrs_wgt_ghost)) then
      nullify(vertex2%nbrs_wgt_ghost)
    endif

    call allocate_gridvertex_nbrs(vertex2)

    do ii=1,nn
      vertex2%nbrs          (ii) = vertex1%nbrs          (ii)
      vertex2%nbrs_face     (ii) = vertex1%nbrs_face     (ii)
      vertex2%nbrs_wgt      (ii) = vertex1%nbrs_wgt      (ii)
      vertex2%nbrs_wgt_ghost(ii) = vertex1%nbrs_wgt_ghost(ii)
    end do

    do ii=1, num_neighbors+1
      vertex2%nbrs_ptr(ii) = vertex1%nbrs_ptr(ii)
    end do

    vertex2%face_number      = vertex1%face_number
    vertex2%number           = vertex1%number
    vertex2%processor_number = vertex1%processor_number 
    vertex2%SpaceCurve       = vertex1%SpaceCurve 

    ! End Routine
    !-------------
    return
  end subroutine copy_gridvertex
  !==================================================================

 
  !==================================================================
  recursive subroutine copy_edgeindex(index2,index1)
    ! copy_edgeindex:
    !
    !==============================================================
    !
    ! Passed Variables
    !--------------------
    type(EdgeIndex_t),intent(out):: index2
    type(EdgeIndex_t),intent(in ):: index1

    if(associated(index1%iyP)) then 
      index2%iyP => index1%iyP
    endif

    if(associated(index1%ixP)) then 
      index2%ixP => index1%ixP
    endif
 
    ! End Routine
    !-------------
    return
  end subroutine copy_edgeindex
  !==================================================================


  !==================================================================
  subroutine FreeGraph(Vertex)
    ! FreeGraph:
    !
    !===============================================================
    !
    ! Passed Variables
    !---------------------
    type(GridVertex_t):: Vertex(:)
    !
    ! Local Values
    !----------------
    integer  ii,nelem

    nelem = SIZE(Vertex)

!JMD     do ii=1,nelem
!JMD        deallocate(Vertex(ii)%wgtV)
!JMD        deallocate(Vertex(ii)%wgtG)
!JMD        deallocate(Vertex(ii)%nbrs)
!JMD     end do

    ! End Routine
    !-------------
    return
  end subroutine FreeGraph
  !==================================================================


  !==================================================================
  function gridedge_search(nvert1, nvert2, edge) result(number)
    ! gridedge_search: search edge list for match
    !
    !=====================================================================
    !
    ! Passed Variables
    !---------------------
    integer         ,intent(in):: nvert1
    integer         ,intent(in):: nvert2
    type(GridEdge_t),intent(in):: edge(:)
    integer                    :: number
    !
    ! Local Values
    !---------------
    integer :: tmp
    integer :: head
    integer :: tail
    integer :: nedge
    integer :: ii

    nedge=SIZE(edge)
    tail =nvert1
    head =nvert2

    if(tail > head) then
      tmp  = tail
      tail = head
      head = tmp
    endif

    do ii=1,nedge
      if((edge(ii)%tail%number == tail).and. &
         (edge(ii)%head%number == head)      ) then
        number=ii
      endif
    end do

    ! End Function
    !-------------
    return
  end function gridedge_search
  !==================================================================


  !==================================================================
  function gridedge_type(edge) result(type)
    ! gridedge_type:
    !
    !=============================================================
    use SE_Options,only: INTERNAL_EDGE, EXTERNAL_EDGE
    !
    ! Passed Variables
    !------------------
    type(GridEdge_t),intent(in):: edge
    integer                    :: type

    if(edge%head%processor_number == edge%tail%processor_number) then
      type=INTERNAL_EDGE
    else
      type=EXTERNAL_EDGE
    endif

    ! End Function
    !-------------
    return
  end function gridedge_type
  !==================================================================


  !==================================================================
  function grid_edge_uses_vertex(Vertex,Edge) result(log)
    ! grid_edge_uses_vertex:
    !
    !========================================================
    !
    ! Passed Variables
    !-------------------
    type(GridVertex_t), intent(in) :: Vertex
    type(GridEdge_t),   intent(in) :: Edge
    logical :: log
    !
    ! Local Values
    !----------------
    integer number

    number = Vertex%number
    if((number == Edge%head%number).or. &
       (number == Edge%tail%number)     )then
      log = .TRUE.
    else
      log = .FALSE.
    endif

    ! End Function
    !-------------
    return
  end function grid_edge_uses_vertex
  !==================================================================


  !==================================================================
  subroutine PrintChecksum(TestPattern,Checksum)
    ! PrintChecksum:
    !
    !================================================================
    use SE_Constants,only: nlev, np
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind),target,intent(in):: TestPattern(:,:,:,:)
    real(kind=real_kind),target,intent(in)::    Checksum(:,:,:,:)
    !
    ! Local Values
    !--------------
    integer ii,kk,ix,iy

    print *
    write (iulog,*) 'checksums:'
    do ii=1,nelemd
      !  Lets start out only looking at the first element
      !--------------------------------------------------
      write(iulog,*)
      do kk=1,nlev
      do iy=1,np
      do ix=1,np
        write(iulog,*)INT(TestPattern(ix,iy,kk,ii))," checksum = ",INT(Checksum(ix,iy,kk,ii))
      end do
      end do
      end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine PrintChecksum
  !==================================================================


  !==================================================================
  subroutine CreateSubGridGraph(Vertex, SVertex, local2global)
    ! CreateSubGridGraph:
    !
    !===============================================================
    !
    ! Passed Variables
    !--------------------
    type(GridVertex_t),intent(in   ):: Vertex      (:)
    type(GridVertex_t),intent(inout):: SVertex     (:)
    integer           ,intent(in   ):: local2global(:)
    !
    ! Local Values
    !--------------
    logical,parameter  :: Debug = .FALSE.
    integer,allocatable:: global2local(:)
    integer            :: nelem,nelem_s,nn,ncount,cnt,pos,orig_start
    integer            :: inbr,ii,ig,jj,kk,new_pos
    
    nelem   = size( Vertex)
    nelem_s = size(SVertex) 

    if(Debug) write(iulog,*)'CreateSubGridGraph: point #1'
    allocate(global2local(nelem))
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #2'

    global2local(:) = 0
    do ii=1,nelem_s
      ig               = local2global(ii)
      global2local(ig) = ii
    end do
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #3'

    do ii=1,nelem_s
      ig = local2global(ii)
      if(Debug) write(iulog,*)'CreateSubGridGraph: point #4'
      call copy_gridvertex(SVertex(ii),Vertex(ig))  !svertex(ii) = vertex(ig)
      nn = SIZE(SVertex(ii)%nbrs(:))

      ! Apply the correction to the neighbors list to 
      ! reflect new subgraph numbers 
      !----------------------------------------------------
      if(Debug) write(iulog,*)'CreateSubGridGraph: point #5'
      orig_start = 1
      do jj=1,num_neighbors
        cnt    = Svertex(ii)%nbrs_ptr(jj+1) - orig_start  !number of neighbors for this direction
        ncount = 0
        do kk=1,cnt
          pos = orig_start + kk-1
          if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.1 ',        &
                                  'size(global2local) Svertex(ii)%nbrs(j) ', &
                                   size(global2local),Svertex(ii)%nbrs(pos)
          inbr = global2local(Svertex(ii)%nbrs(pos))
          if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.2'
          if(inbr .gt. 0) then 
            new_pos = Svertex(ii)%nbrs_ptr(jj) + ncount
            Svertex(ii)%nbrs          (new_pos) = inbr
            Svertex(ii)%nbrs_face     (new_pos) = Svertex(ii)%nbrs_face     (pos)
            Svertex(ii)%nbrs_wgt      (new_pos) = Svertex(ii)%nbrs_wgt      (pos)
            Svertex(ii)%nbrs_wgt_ghost(new_pos) = Svertex(ii)%nbrs_wgt_ghost(pos)
            ncount = ncount+1
          endif
        end do ! kk = 1, cnt

        !set neighbors ptr
        !--------------------
        orig_start                 =  Svertex(ii)%nbrs_ptr(jj+1)
        Svertex(ii)%nbrs_ptr(jj+1) =  Svertex(ii)%nbrs_ptr(jj) + ncount
        if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.3'
      end do ! jj=1,num_neighbors
      if(Debug) write(iulog,*)'CreateSubGridGraph: point #6'
      Svertex(ii)%number = ii
    end do ! ii=1,nelem_s

    if(Debug) write(iulog,*)'CreateSubGridGraph: point #7'
    deallocate(global2local)
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #8'

    ! End Routine
    !-------------
    return
  end subroutine CreateSubGridGraph
  !==================================================================


  !==================================================================
  subroutine PrintGridEdge(Edge)
    ! PrintGridEdge:
    !
    !===============================================================
    !
    ! Passed Variables
    !-------------------
    type(GridEdge_t),intent(in):: Edge(:)
    !
    ! Local Values
    !--------------
    integer jj,nedge,ii,wgtP

    nedge = SIZE(Edge)
    write(iulog,95)

    do jj=1,nedge
      ii=Edge(jj)%tail_face

      ! map to correct location - for now all on same nbr side have 
      ! same wgt, so take the first one
      !------------------------------------------------------------
      ii  = Edge(jj)%tail%nbrs_ptr(ii)
      wgtP= Edge(jj)%tail%nbrs_wgt(ii)       
      write(iulog,100) jj,                                                   &
               Edge(jj)%tail%number,Edge(jj)%tail_face,wgtP,                 &
               Edge(jj)%head%number,Edge(jj)%head_face,gridedge_type(Edge(jj))
    end do

  95 format(5x,'GRIDEDGE #',3x,'Tail (face)',5x,'Head (face)',3x,'Type')
 100 format(10x,I6,8x,I4,1x,'(',I1,')  --',I2,'--> ',I6,1x,'(',I1,')',5x,'[',I1,']')

    ! End Routine
    !-------------
    return
  end subroutine PrintGridEdge
  !==================================================================


  !==================================================================
  subroutine set_GridVertex_number(elem,number)
    ! set_GridVertex_neighbors: Set global element number for element elem
    !
    !==============================================================
    !
    ! Passed Variables
    !-------------------
    type(GridVertex_t):: elem
    integer           :: number
  
    elem%number=number

    ! End Routine
    !-------------
    return
  end subroutine set_GridVertex_number
  !==================================================================


  !==================================================================
  subroutine PrintGridVertex(Vertex)
    ! PrintGridVertex:
    !
    !===============================================================
    !
    ! Passed variables
    !------------------
    type(GridVertex_t),intent(in),target:: Vertex(:)
    !
    ! Local Values
    !-----------------
    integer,pointer:: np(:)
  
    integer:: n_west, n_east, n_south, n_north, n_swest, n_seast, n_nwest, n_neast
    integer:: w_west, w_east, w_south, w_north, w_swest, w_seast, w_nwest, w_neast
    integer:: nn, print_buf(90), nbr(8), jj, kk, start, cnt, nbrs_cnt(8)
    integer:: ii,nvert

    nbr   = (/ west, east, south, north, swest, seast, nwest, neast/)
    nvert = SIZE(Vertex)
        
    write(iulog,98)
    do ii=1,nvert
      print_buf(:) = 0
      nbrs_cnt (:) = 0
      cnt          = 1
      np =>  Vertex(ii)%nbrs_ptr  !alias
      do jj=1,num_neighbors
        nn           = np(nbr(jj)+1) - np(nbr(jj)) !num neigbors in that directions
        start        = np(nbr(jj))                 !start in array
        nbrs_cnt(jj) = nn
        do kk=1,nn
          print_buf(cnt  ) = Vertex(ii)%nbrs     (start+kk-1)
          print_buf(cnt+1) = Vertex(ii)%nbrs_wgt (start+kk-1)
          print_buf(cnt+2) = Vertex(ii)%nbrs_face(start+kk-1)
          cnt = cnt + 3
        end do
      end do ! jj=1,num_neighbors
      write(iulog,991) Vertex(ii)%number     , Vertex(ii)%processor_number, &
                       Vertex(ii)%face_number, print_buf(1:cnt-1)
      write(iulog,992) nbrs_cnt(1:8)
    end do ! ii=1,nvert

   98 format(5x,'GRIDVERTEX #',2x,'PART',2x,'DEG',4x,'W',9x,'E',9x, &
                'S',9x,'N',9x,'SW',9x,'SE',9x,'NW',9x,'NE')
  991 format(10x,I3,8x,I4,8x,I4,2x,30(1x,I4,1x,'(',I2,I2,')'))
  992 format(30x,'nbrs_cnt:', 2x,8(1x,I4))

    ! End Routine
    !-------------
    return
  end subroutine PrintGridVertex
  !==================================================================


  !==================================================================
  subroutine CheckGridNeighbors(Vertex)
    ! CheckGridNeighbors:
    !
    !================================================================
    !
    ! Passed Variables
    !---------------------
    type(GridVertex_t),intent(in):: Vertex(:)
    !
    ! Local Values
    !---------------
    integer ii,jj,kk,ll,mm,nnbrs,inbrs,nvert

    nvert = SIZE(Vertex)
    do ii=1,nvert
      nnbrs = SIZE(Vertex(ii)%nbrs)
      do jj=1,nnbrs
        inbrs = Vertex(ii)%nbrs(jj)
        if(inbrs > 0) then
          do kk=1,nnbrs
            if((inbrs .eq. Vertex(ii)%nbrs(kk)).and.(jj /= kk)) then
              write(iulog,*)'CheckGridNeighbors: ERROR identical neighbors detected  for Vertex ',ii
            endif
          end do
        endif
      end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine CheckGridNeighbors
  !==================================================================


  !==================================================================
  subroutine initgridedge(GridEdge,GridVertex)
    ! initgridedge:
    !
    !================================================================
    use err_exit  ,only: endrun
    !
    ! Passed Variables
    !--------------------
    type (GridEdge_t)  ,intent(inout)       :: GridEdge  (:)
    type (GridVertex_t),intent(in   ),target:: GridVertex(:)
    !
    ! Local Values
    !-----------------
    logical:: Verbose=.FALSE.
    integer:: nelem,nelem_edge,inbr
    integer:: mynbr_cnt, cnt, mystart, start
    integer:: ii,jj,kk,iptr,mm,nn,wgtV,wgtP

    nelem              = SIZE(GridVertex)
    nelem_edge         = SIZE(GridEdge)
    GridEdge(:)%reverse= .FALSE.
    iptr               = 1

    do jj=1,nelem
    do ii=1,num_neighbors    
      mynbr_cnt = GridVertex(jj)%nbrs_ptr(ii+1) &  ! length of neighbor location  
                 -GridVertex(jj)%nbrs_ptr(ii)
      mystart   = GridVertex(jj)%nbrs_ptr(ii) 
      do mm=0,mynbr_cnt-1
        ! Do this only if has a non-zero weight
        !-----------------------------------------
        if(GridVertex(jj)%nbrs_wgt(mystart + mm) .gt. 0) then
          if(nelem_edge < iptr) then
            call endrun('Error in initgridedge: Number of edges greater than expected.')
          endif
          GridEdge(iptr)%tail     => GridVertex(jj)
          GridEdge(iptr)%tail_face=  mystart+mm              ! needs to be mystart+m (location in array)
          GridEdge(iptr)%tail_dir =  ii*max_corner_elem + mm ! conversion needed for setcycle
          inbr                    =  GridVertex(jj)%nbrs(mystart+mm)
          GridEdge(iptr)%head     => GridVertex(inbr)

          ! Need this awful piece of code to determine
          ! which "face" of the neighbor element the
          ! edge links (i.e. the "head_face")
          !-------------------------------------------
          do kk=1,num_neighbors
            cnt   = GridVertex(inbr)%nbrs_ptr(kk+1) -GridVertex(inbr)%nbrs_ptr(kk)                     
            start = GridVertex(inbr)%nbrs_ptr(kk)
            do nn=0,cnt-1
              if(GridVertex(inbr)%nbrs(start+nn) == GridVertex(jj)%number) then
                GridEdge(iptr)%head_face= start+nn              ! needs to be start+n (location in array)
                GridEdge(iptr)%head_dir = kk*max_corner_elem+nn ! conversion (un-done in setcycle)
              endif
            end do
          end do
          iptr=iptr+1
        endif
      end do ! mm=0,mynbr_cnt-1
    end do !   ii=1,num_neighbors    
    end do !   jj=1,nelem

    if(nelem_edge+1 /= iptr) then
      call endrun('Error in initgridedge: Number of edges less than expected.')
    endif

    if(Verbose) then
      write(iulog,*) " "
      write(iulog,*) "element edge tail,head list: (TEST)"
      do ii=1,nelem_edge
        write(iulog,*) GridEdge(ii)%tail%number,GridEdge(ii)%head%number
      end do
      write(iulog,*) " "
      write(iulog,*) "element edge tail_face, head_face list: (TEST)"
      do ii=1,nelem_edge
        write(iulog,*) GridEdge(ii)%tail_face,GridEdge(ii)%head_face
      end do
    endif

    ! End Routine
    !-------------
    return
  end subroutine initgridedge
  !==================================================================

end module GridGraph_mod

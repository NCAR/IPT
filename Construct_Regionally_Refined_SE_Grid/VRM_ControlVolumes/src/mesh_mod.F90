#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module mesh_mod

  ! Useful modules
  !----------------
  use SE_Constants,only: real_kind, long_kind
  use SE_Constants,only: DD_PI
  use SE_Constants,only: MAX_FILE_LEN
  use netcdf                                       ! _EXTERNAL

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private 

  private:: MXSTLN
  private:: nfaces
  private:: nInnerElemEdge

  public :: MeshUseMeshFile
  private:: p_mesh_file_name
  private:: p_ncid
  private:: p_number_elements 
  private:: p_number_elements_per_face
  private:: p_number_blocks 
  private:: p_number_nodes 
  private:: p_number_dimensions 
  private:: p_number_neighbor_edges 
  private:: p_node_coordinates
  private:: p_connectivity
  private:: p_elem_block_ids

  private:: handle_error
  private:: open_mesh_file
  private:: close_mesh_file
  private:: get_number_of_dimensions
  private:: get_number_of_elements
  private:: get_number_of_nodes
  private:: get_number_of_element_blocks
  private:: get_number_of_elements_per_face
  private:: get_block_ids                    ! Not used will eliminate later
  private:: get_face_connectivity
  private:: get_node_multiplicity
  private:: get_node_coordinates
  private:: get_2D_sub_coordinate_indexes
  private:: mesh_connectivity
  private:: create_index_table
  private:: find_side_neighbors
  private:: smallest_diameter_element
  private:: cube_to_cube_coordinates
  private:: sphere_to_cube_coordinates
  private:: cube_face_element_centroids
  private:: initialize_space_filling_curve
  private:: find_corner_neighbors
  public :: MeshOpen
  public :: MeshClose  
  public :: MeshPrint          ! show the contents of the Mesh after it has been loaded into the module
  public :: MeshCubeTopology   ! called afer MeshOpen
  public :: MeshSetCoordinates ! called after MeshCubeTopology    
  public :: MeshCubeEdgeCount  ! called anytime afer MeshOpen
  public :: MeshCubeElemCount  ! called anytime afer MeshOpen
  private:: test_private_methods

  ! PFC: Hack to set global option values
  !--------------------------------------
  public :: init_mesh_mod
  private:: max_elements_attached_to_node
  private:: max_corner_elem
  integer:: max_elements_attached_to_node
  integer:: max_corner_elem

  ! Parameters
  !--------------
  integer,parameter:: MXSTLN          = 32
  integer,parameter:: nfaces          = 6  ! number of faces on the cube
  integer,parameter:: nInnerElemEdge  = 8  ! number of edges for an interior element

  ! Global Data
  !-------------
  logical                         :: MeshUseMeshFile = .false.
  character (len=MAX_FILE_LEN)    :: p_mesh_file_name
  integer                         :: p_ncid
  integer                         :: p_number_elements 
  integer                         :: p_number_elements_per_face
  integer                         :: p_number_blocks 
  integer                         :: p_number_nodes 
  integer                         :: p_number_dimensions 
  integer                         :: p_number_neighbor_edges 
  real(kind=real_kind),allocatable:: p_node_coordinates(:,:) 
  integer             ,allocatable:: p_connectivity    (:,:)
  integer                         :: p_elem_block_ids            ! Not used will eliminate later


contains
  !==================================================================
  subroutine init_mesh_mod(I_SEopt)
    ! PFC: Hack to set global option values
    !--------------------------------------
    use SE_Options,only: SEoptions_t
    ! Passed Variables
    !------------------
    type(SEoptions_t),intent(in):: I_SEopt

    max_elements_attached_to_node = I_SEopt%max_elements_attached_to_node
    max_corner_elem               = I_SEopt%max_corner_elem
 
    print *,' PFC: max_elements_attached_to_node =',max_elements_attached_to_node,I_SEopt%max_elements_attached_to_node
    print *,' PFC: max_corner_elem =',max_corner_elem,I_SEopt%max_corner_elem

    ! End Routine
    !-------------
    return
  end subroutine init_mesh_mod
  !==================================================================


  !======================================================================
  subroutine handle_error(status, file, line)
    !  subroutine handle_error:
    !==================================================================
    use err_exit,only: endrun
    !
    ! Passed Variables
    !-----------------
    integer         ,intent(in):: status
    character(len=*),intent(in):: file
    integer         ,intent(in):: line

    print *, file,':',line,': ',trim(nf90_strerror(status))
    call endrun("Terminating program due to netcdf error while obtaining mesh information, please see message in standard output.")

    ! End Routine
    !---------------
    return
  end subroutine handle_error
  !======================================================================

  
  !======================================================================
  subroutine open_mesh_file() 
    !  open_mesh_file:
    !
    ! Open the netcdf file containing the mesh.
    ! Assign the holder to the file to p_ncid so everyone else knows
    ! how to use it without passing the argument around.
    !======================================================================
    !
    ! Local Values
    !---------------
    integer:: status

    status = nf90_open(p_mesh_file_name, NF90_NOWRITE, p_ncid)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
 
    MeshUseMeshFile = .true. 

    ! End Routine
    !---------------
    return
  end subroutine open_mesh_file
  !======================================================================


  !======================================================================
  subroutine close_mesh_file() 
    ! close_mesh_file:
    !
    !======================================================================
    !
    ! Local Values
    !-----------------
    integer:: status
    
    status = nf90_close(p_ncid)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
    
    ! End Routine
    !---------------
    return
  end subroutine close_mesh_file
  !======================================================================


  !======================================================================
  function get_number_of_dimensions() result(number_dimensions)
    ! get_number_of_dimensions:
    !
    !======================================================================
    ! 
    ! Passed Variables
    !--------------------
    integer:: number_dimensions
    ! 
    ! Local Values
    !--------------
    integer:: status, number_of_dim_id

    ! Get the id of 'num_elem', if such dimension is 
    ! not there panic and quit :P
    !------------------------------------------------
    status = nf90_inq_dimid(p_ncid, "num_dim", number_of_dim_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! How many values for 'num_elem' are there?
    !----------------------------------------------
    status = nf90_inquire_dimension(p_ncid, number_of_dim_id, len = number_dimensions)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! End Function
    !---------------
    return
  end function get_number_of_dimensions
  !======================================================================


  !======================================================================
  function get_number_of_elements() result(number_elements)
    ! get_number_of_elements:
    !
    !======================================================================
    ! Passed Variables
    !----------------------
    integer:: number_elements 
    !
    ! Local Values
    !----------------
    integer:: status, number_of_elements_id

    ! Get the id of 'num_elem', if such dimension 
    ! is not there panic and quit :P
    !----------------------------------------------
    status = nf90_inq_dimid(p_ncid, "num_elem", number_of_elements_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! How many values for 'num_elem' are there?
    !-------------------------------------------
    status = nf90_inquire_dimension(p_ncid, number_of_elements_id, len = number_elements)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! End Function
    !---------------
    return
  end function get_number_of_elements
  !======================================================================


  !======================================================================
  function get_number_of_nodes() result(number_nodes)
    ! get_number_of_nodes:
    !
    !======================================================================
    ! 
    ! Passed Variables
    !------------------
    integer:: number_nodes
    !
    ! Local Values
    !----------------
    integer:: status, number_of_nodes_id

    ! Get the id of 'num_nodes', if such dimension 
    ! is not there panic and quit :P
    !-----------------------------------------------
    status = nf90_inq_dimid(p_ncid, "num_nodes", number_of_nodes_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! How many values for 'num_nodes' are there?
    !----------------------------------------------
    status = nf90_inquire_dimension(p_ncid, number_of_nodes_id, len = number_nodes)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! End Function
    !---------------
    return
  end function get_number_of_nodes
  !======================================================================


  !======================================================================
  function get_number_of_element_blocks() result(number_element_blocks)
    ! get_number_of_element_blocks:
    !
    !======================================================================
    use err_exit,only: endrun
    !
    ! Passed Variables
    !------------------
    integer:: number_element_blocks 
    !
    ! Local Values
    !---------------
    integer:: status, number_of_element_blocks_id
    
    ! Get the id of 'num_el_blk', if such dimension 
    ! is not there panic and quit :P
    !--------------------------------------------------
    status = nf90_inq_dimid(p_ncid, "num_el_blk", number_of_element_blocks_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! How many values for 'num_el_blk' are there?
    !--------------------------------------------
    status = nf90_inquire_dimension(p_ncid, number_of_element_blocks_id, len = number_element_blocks)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    if(number_element_blocks /= 1) then
      if(number_element_blocks /= 6  ) then
        call endrun('Reading cube-sphere from input file is not supported')
      else
        call endrun('Number of elements blocks not exactly 1 (sphere) or 6 (cube)')
      endif
    endif

    ! End Function
    !---------------
    return
  end function get_number_of_element_blocks
  !======================================================================


  !======================================================================
  function get_number_of_elements_per_face() result(number_elements_per_face)
    ! get_number_of_elements_per_face:
    !
    !======================================================================
    use err_exit,only: endrun
    ! 
    ! Passed Variables
    !--------------------
    integer:: number_elements_per_face
    !
    ! Local Values
    !--------------
    integer              :: face_num                ! For each of the face, we get the information
    character(len=MXSTLN):: element_type            ! Each face is composed of elements of certain type
    integer              :: number_elements_in_face ! How many elements in this face
    integer              :: num_nodes_per_elem      ! How many nodes in each element
    integer              :: number_of_attributes    ! How many attributes in the face
    integer              :: status, dimension_id

    if(p_number_blocks == 0)  then
      call endrun('get_number_of_elements_per_face called before MeshOpen')
    elseif(p_number_blocks == 1) then 
      ! we are in the presence of a sphere
      ! First we get sure the number of nodes per element is four
      !-----------------------------------------------------------------
      status = nf90_inq_dimid(p_ncid, "num_nod_per_el1", dimension_id)
      if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
      status = nf90_inquire_dimension(p_ncid, dimension_id, len =  num_nodes_per_elem)
      if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
      if (num_nodes_per_elem /= 4)  call endrun('Number of nodes per element is not four')

      ! now we check how many elements there are in the face
      !-------------------------------------------------------
      status = nf90_inq_dimid(p_ncid, "num_el_in_blk1", dimension_id)
      if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
      status = nf90_inquire_dimension(p_ncid, dimension_id, len = number_elements_in_face)
      if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
      number_elements_per_face =  number_elements_in_face
    elseif(p_number_blocks == 6) then 
      ! we are in the presence of a cube-sphere
      !----------------------------------------------
      call endrun('Reading a mesh for a cube-sphere is not supported')
    else
      call endrun('Number of elements blocks not exactly 1 (sphere) or 6 (cube)')
    endif

    ! End Function
    !---------------
    return
  end function get_number_of_elements_per_face
  !======================================================================


  !======================================================================
  function get_block_ids(idexo) result(block_ids)
    ! get_block_ids: This function is used to set the value of 
    !                p_elem_block_ids  but such variable is never used
    !======================================================================
    use err_exit,only: endrun
    !
    ! Passed Variables
    !------------------
    integer(kind=long_kind),intent(in):: idexo
    integer(kind=long_kind)           :: block_ids(p_number_blocks)

    block_ids = 0

    ! End Function
    !---------------
    return
  end function get_block_ids
  !======================================================================


  !======================================================================
  subroutine get_face_connectivity() 
    ! get_face_connectivity:
    !
    !======================================================================
    use err_exit,only: endrun
    !
    ! Local Values
    !--------------
    integer:: var_id, status

    status = nf90_inq_varid(p_ncid, "connect1", var_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
    status = nf90_get_var(p_ncid, var_id, p_connectivity)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
    
    ! End Routine
    !---------------
    return
  end subroutine get_face_connectivity
  !======================================================================


  !======================================================================
  subroutine get_node_multiplicity(node_multiplicity) 
    ! get_node_multiplicity:
    !
    !======================================================================
    use err_exit  ,only: endrun
    !
    ! Passed Variables
    !---------------------
    integer,intent(out) :: node_multiplicity(:)
    ! 
    ! Local Values
    !---------------
    integer:: node_num(4)
    integer:: kk, number_nodes

    node_multiplicity(:)= 0
    number_nodes        = SIZE(node_multiplicity)

    ! check this external buffer was allocated correctly
    !-----------------------------------------------------
    if(number_nodes /= p_number_nodes) then
      call endrun('Number of nodes does not matches size of node multiplicity array')
    endif

    ! for each node, we have for four other nodes
    !----------------------------------------------
    if((minval(p_connectivity) < 1).or.(number_nodes < maxval(p_connectivity))) then
      call endrun('get_node_multiplicity: Node number less than 1 or greater than max.')
    endif
    
    do kk=1,p_number_elements_per_face
      node_num                    = p_connectivity   (:,kk)
      node_multiplicity(node_num) = node_multiplicity(node_num) + 1
    end do

    if((    minval(node_multiplicity) < 3                        ).or. &
       (max_elements_attached_to_node < maxval(node_multiplicity))     ) then
      print *, 'minval(node_multiplicity)', minval(node_multiplicity)
      print *, 'maxval(node_multiplicity)', maxval(node_multiplicity),&
               ' and max_elements_attached_to_node ',max_elements_attached_to_node
      call endrun('get_node_multiplicity: Number of elements attached to node less than 3 or greater than maximum.')
    endif
    
    ! End Routine
    !---------------
    return
  end subroutine get_node_multiplicity
  !======================================================================


  !======================================================================
  subroutine get_node_coordinates ()
    ! get_node_coordinates:
    !
    !======================================================================
    use coordinate_systems_mod,only: cartesian3D_t
    use err_exit              ,only: endrun
    !
    ! Local Values
    !----------------
    integer:: var_id, status

    status = nf90_inq_varid(p_ncid, "coord", var_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
    status = nf90_get_var(p_ncid, var_id, p_node_coordinates)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
    
    ! End Routine
    !---------------
    return
  end subroutine get_node_coordinates
  !======================================================================


  !======================================================================
  subroutine get_2D_sub_coordinate_indexes(xx, yy, sgnx, sgny, face_no) 
    ! get_2D_sub_coordinate_indexes:
    !======================================================================
    !
    ! Passed Variables
    !-------------------
    integer, intent(out):: xx,yy
    integer, intent(out):: sgnx, sgny
    integer, intent(in ):: face_no

    if((face_no == 1).or.(face_no == 3)) then
      xx = 2
      yy = 3
    elseif((face_no == 2).or.(face_no == 4)) then
      xx = 1
      yy = 3
    else
      xx = 2
      yy = 1
    endif

    if((face_no == 1).or.(face_no == 4).or.(face_no == 5)) then
      sgnx =  1
      sgny =  1
    elseif((face_no == 2).or.(face_no == 3)) then
      sgnx = -1
      sgny =  1
    else  
      sgnx =  1
      sgny = -1
    endif
    
    ! End Routine
    !---------------
    return
  end subroutine get_2D_sub_coordinate_indexes
  !======================================================================
  

  !======================================================================
  subroutine  mesh_connectivity (connect) 
    ! mesh_connectivity: puts the transpose of p_connectivity into connect
    !
    !======================================================================
    use err_exit,only: endrun
    !
    ! Passed Variables
    !------------------
    integer,intent(out):: connect(p_number_elements,4)
    !
    ! Local Values
    !---------------
    integer kk, jj 

    if(0 == p_number_blocks) call endrun('mesh_connectivity called before MeshOpen')

    jj=0
    do kk=1, p_number_elements_per_face
      jj=jj+1
      connect(jj,:) = p_connectivity(:,kk)
    end do
      
    if(jj /= p_number_elements) then
      call endrun('mesh_connectivity: Number of elements in side sets not equal to total elements')
    endif
   
    if((minval(connect) < 1).or.(maxval(connect) > p_number_nodes)) then
      call endrun('mesh_connectivity: Node number out of bounds')
    endif
    
    ! End Routine
    !---------------
    return
  end subroutine mesh_connectivity 
  !======================================================================


  !======================================================================
  subroutine create_index_table(index_table, element_nodes)
    ! create_index_table:
    !
    ! this is needed to detremine side and corner neighbors
    !======================================================================
    use err_exit  ,only: endrun
    !
    ! Passed Variables
    !-------------------
    integer,allocatable,intent(inout):: index_table(:,:) 
    integer            ,intent(in   ):: element_nodes(p_number_elements, 4)
    !
    ! Local Values
    !----------------
    integer:: cnt, cnt_index, node
    integer:: kk, ll 

    ! Create an index table so that we can find neighbors on O(n)
    ! so for each node, we want to know which elements it is part of
    !----------------------------------------------------------------
    allocate(index_table(p_number_nodes, max_elements_attached_to_node + 1))
   
    ! the last column in the index table is a count of the number of elements
    !-----------------------------------------------------------------------------
    index_table = 0
    cnt_index =  max_elements_attached_to_node + 1
     
    do kk=1,p_number_elements
    do ll=1,4 
      node = element_nodes(kk, ll)          !the node
      cnt  = index_table(node, cnt_index)    !how many elements for that node already in table
      cnt  = cnt + 1                         !increment since we are adding an element
      if(cnt >  max_elements_attached_to_node) then
        call endrun('Found a node in too many elements.')
      endif
      index_table(node, cnt_index) = cnt  
      index_table(node, cnt      ) = kk           !put the element in the indextable
    end do
    end do
    
    ! End Routine
    !---------------
    return
  end subroutine create_index_table
  !======================================================================


  !======================================================================
  subroutine find_side_neighbors(GridVertex, normal_to_homme_ordering, &
                                 element_nodes, edge_wgt, index_table)
    ! find_side_neighbors:
    !
    ! find the element neighbors to the n,s,e,w and put them in GridVertex_t
    !  (only 1 neighbor to the n,s,e,w)
    !======================================================================
    use coordinate_systems_mod,only: cartesian3D_t
    use gridgraph_mod         ,only: GridVertex_t
    use err_exit              ,only: endrun
    !
    ! Passed Variables
    !------------------
    type (GridVertex_t),intent(inout):: GridVertex(:)
    integer            ,intent(in   ):: normal_to_homme_ordering(8)
    integer            ,intent(in   ):: element_nodes(p_number_elements, 4)
    integer            ,intent(in   ):: edge_wgt
    integer            ,intent(in   ):: index_table(:,:) 
    !
    ! Local Values
    !-------------
    integer:: i_node(2), my_node(2) 
    integer:: neighbor, direction, init_size
    integer:: jj,kk,ll,ii, mm
    integer:: i_elem, jump, end_i
    integer:: loc, node, cnt_index, cnt, a_count(2)
    logical:: found

    if(0 == p_number_blocks)  call endrun('find_side_neighbors called before MeshOpen')
   
    !the last column in the index table is a count of the number of elements
    !------------------------------------------------------------------------
    cnt_index =  max_elements_attached_to_node + 1
     
    !use index table to find neighbors for each element k 
    !----------------------------------------------------
    do kk=1,p_number_elements  

      ! set the side weights
      !----------------------------------------
      GridVertex(kk)%nbrs_wgt(1:4) = edge_wgt

      ! loop through the four sides
      !----------------------------
      do ll=1,4
        jump = normal_to_homme_ordering(ll)
        loc  = GridVertex(kk)%nbrs_ptr(jump)
        if(GridVertex(kk)%nbrs(loc) == 0) then  

          !if side is not set yet, then look for side element
          !---------------------------------------------------
          found      = .false.
          neighbor   = 0
          my_node(1) = element_nodes(kk, ll)
          a_count(1) = index_table  (my_node(1), cnt_index)
          my_node(2) = element_nodes(kk, mod(ll,4)+1)
          a_count(2) = index_table  (my_node(2), cnt_index)

          ! loop through the elements that are in the index table for each node
          ! and find the element number and direction of the side neighbor
          !---------------------------------------------------------------------
          do mm=1,2
            if(found) exit
            end_i = a_count(mm)
            do ii=1,end_i
              if(found) exit
              i_elem = index_table(my_node(mm),ii)

              ! k is the element we are setting sides for
              !-------------------------------------------
              if(i_elem /= kk) then 
                ! loop through each of i_elem's four sides
                !------------------------------------------
                do jj=1,4 
                  i_node(1) = element_nodes(i_elem,          jj)
                  i_node(2) = element_nodes(i_elem, mod(jj,4)+1)
                  if(((i_node(1) == my_node(2)).and.(i_node(2) == my_node(1))) .or. &
                     ((i_node(1) == my_node(1)).and.(i_node(2) == my_node(2)))      ) then
                    ! found a match
                    !----------------
                    neighbor  = i_elem
                    direction = jj
                    found     = .true.
                    exit
                  endif
                end do ! j loop
              endif
            end do ! i loop
          end do !m loop

          if(neighbor == 0) then
            call endrun('find_side_neighbor: Neighbor not found! Every side should have a neighbor.') 
          endif

          GridVertex(kk)%nbrs(loc)      = neighbor         
          jump                          = normal_to_homme_ordering(direction)
          loc                           = GridVertex(neighbor)%nbrs_ptr(jump)
          GridVertex(neighbor)%nbrs(loc)= kk
        endif
      enddo !  ll loop => 4 sides
    enddo ! k loop: each element
    
    do kk=1,p_number_elements
    do ll=1,4
      if( 0 == GridVertex(kk)%nbrs(ll)) then
        call endrun('Found one side of one element witout a neighbor.  Bummer!') 
      endif
    end do
    end do
    
    ! End Routine
    !---------------
    return
  end subroutine find_side_neighbors
  !======================================================================


  !======================================================================
  function smallest_diameter_element(element_nodes) result(min_diameter)
    ! smallest_diameter_element:
    !
    !======================================================================
    use err_exit,only: endrun
    !
    ! Passed Variables
    !---------------------
    integer,intent(in):: element_nodes(:,:)
    real              :: min_diameter
    !
    ! Local Values
    !--------------
    integer             :: node_numbers(4)
    real(kind=real_kind):: coordinates (4,3)
    real                :: xx(3), yy(3), rr(3), dd
    integer             :: ii, jj
    
    if(SIZE(element_nodes,dim=1) /= p_number_elements) then
      call endrun('smallest_diameter_element:Element count check failed in &
                    &exodus_mesh. Connectivity array length not equal to number of elements.')
    endif
    if( p_number_elements_per_face /= p_number_elements) then
      call endrun('smallest_diameter_element: Element count check failed in &
                    &exodus_mesh. Element array length not equal to sum of face.')
    endif
    
    min_diameter = 9999999.
    do ii=1, p_number_elements  
      node_numbers = element_nodes(ii,:)    
      coordinates  = p_node_coordinates(node_numbers,:)

      ! smallest side length
      !--------------------
      do jj=1,4
        xx = coordinates(jj         ,:)
        yy = coordinates(1+MOD(jj,4),:)
        rr = xx-yy
        dd  = dot_product(rr,rr)
        if(dd < min_diameter ) then
          min_diameter = dd
        endif
      end do

      ! smallest diameter length
      !--------------------------
      do jj=1,2
        xx = coordinates(jj         ,:)
        yy = coordinates(2+MOD(jj,4),:)
        rr = xx-yy
        dd  = dot_product(rr,rr)
        if(dd < min_diameter ) then
          min_diameter = dd
        endif
      end do
    end do
    min_diameter = SQRT(min_diameter)

    ! End Function
    !---------------
    return
  end function smallest_diameter_element
  !======================================================================
  

  !======================================================================
  subroutine cube_to_cube_coordinates (cube_coor, node_coor, face_number)
    ! cube_to_cube_coordinates:
    !
    !======================================================================
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind),intent(out):: cube_coor(4,2)
    real(kind=real_kind),intent(in ):: node_coor(4,3)
    integer             ,intent(in ):: face_number
    !
    ! Local Values
    !-------------
    real(kind=real_kind):: test_coor(4,2)
    integer             :: x_index, y_index, sgnx, sgny
    
    call get_2D_sub_coordinate_indexes(x_index, y_index, sgnx, sgny, face_number)
    cube_coor(:,1) = sgnx*node_coor(:,x_index)
    cube_coor(:,2) = sgny*node_coor(:,y_index)
    
    ! End Routine
    !---------------
    return
  end subroutine cube_to_cube_coordinates
  !======================================================================
  

  !======================================================================
  subroutine sphere_to_cube_coordinates (cube_coor, node_coor, face_number)
    ! sphere_to_cube_coordinates:
    !
    !======================================================================
    use coordinate_systems_mod,only: cartesian3D_t, cartesian2d_t, spherical_polar_t, &
                                     change_coordinates, sphere2cubedsphere
    !
    ! Passed Variables
    !--------------------
    real(kind=real_kind),intent(out):: cube_coor(4,2)
    real(kind=real_kind),intent(in ):: node_coor(4,3)
    integer             ,intent(in ):: face_number
    !
    ! Local Values
    !----------------
    type(cartesian2d_t):: cart(4)
    integer            :: ii
    
    do ii=1,4 
      cart(ii) = sphere2cubedsphere(change_coordinates(node_coor(ii,:)), face_number)
    end do
    cube_coor(:,1) = cart(:)%x
    cube_coor(:,2) = cart(:)%y
    
    ! End Routine
    !---------------
    return
  end subroutine sphere_to_cube_coordinates
  !======================================================================
  

  !======================================================================
  subroutine cube_face_element_centroids(centroids, face_numbers, element_nodes)
    ! cube_face_element_centroids:
    !
    !======================================================================
    use err_exit,only: endrun
    !
    ! Passed Variables
    !------------------
    real   ,intent(out):: centroids   (p_number_elements,2)
    integer,intent(in ):: face_numbers(p_number_elements)
    integer,intent(in ):: element_nodes(:,:)
    !
    ! Local Values
    !----------------
    real(kind=real_kind):: coordinates(4,3) 
    real(kind=real_kind):: cube_coor  (4,2) 
    integer             :: ii, node_numbers(4)
    
    if(0 == p_number_blocks)  call endrun('cube_face_element_centroids called before MeshOpen')
    if(SIZE(element_nodes,dim=1) /= p_number_elements) then
      call endrun('cube_face_element_centroids:Element count check failed in &
                    &exodus_mesh. Connectivity array length not equal to number of elements.')
    endif
    if( p_number_elements_per_face /= p_number_elements ) then
       call endrun('cube_face_element_centroids: Element count check failed in &
                     &exodus_mesh. Element array length not equal to sum of face.')
    endif
    
    do ii=1, p_number_elements  
      node_numbers = element_nodes(ii,:)    
      coordinates  = p_node_coordinates(node_numbers,:)
      if(6 == p_number_blocks) then
        call cube_to_cube_coordinates  (cube_coor, coordinates, face_numbers(ii))
      else
        call sphere_to_cube_coordinates(cube_coor, coordinates, face_numbers(ii))
      endif
      centroids(ii,:) = SUM(cube_coor,dim=1)/4.0
    end do
    
    ! End Routine
    !---------------
    return
  end subroutine cube_face_element_centroids
  !======================================================================
  

  !======================================================================
  subroutine initialize_space_filling_curve(GridVertex, element_nodes)
    ! initialize_space_filling_curve:
    !
    !======================================================================
    use gridgraph_mod ,only: GridVertex_t
    use err_exit      ,only: endrun
    use spacecurve_mod,only: GenspaceCurve
    !
    ! Passed Variables
    !--------------------
    type (GridVertex_t), intent(inout) :: GridVertex(:)
    integer            , intent(in)    :: element_nodes(:,:)
    !
    ! Local Values
    !--------------
    real               :: centroids   (p_number_elements,2)
    integer            :: face_numbers(p_number_elements)
    integer,allocatable:: Mesh2(:,:),Mesh2_map(:,:),sfcij(:,:)
    real               :: xx, yy, hh
    integer            :: ii, jj, i2, j2, ne, ne2
    integer            :: sfc_index, face, nelem
    
    if(SIZE(GridVertex) /= p_number_elements) then
      call endrun('initialize_space_filling_curve:Element count check failed &
                    &in exodus_mesh. Vertex array length not equal to number of elements.')
    endif
    if(SIZE(element_nodes,dim=1) /= p_number_elements) then
      call endrun('initialize_space_filling_curve:Element count check failed &
                    &in exodus_mesh. Connectivity array length not equal to number of elements.')
    endif
    
    face_numbers(:) = GridVertex(:)%face_number
    hh              = smallest_diameter_element(element_nodes)

    call cube_face_element_centroids(centroids, face_numbers, element_nodes)
    
    if (hh<.00001) call endrun('initialize_space_filling_curve: Unreasonably small element found. less than .00001')
    
    ne = CEILING(0.5*DD_PI/(hh/2));
    
    ! find the smallest ne2 which is a power of 2 and ne2>ne
    !-------------------------------------------------------
    ne2=2**ceiling( log(real(ne))/log(2d0) )
    if (ne2<ne) call endrun('initialize_space_filling_curve: Fatel SFC error')
    
    allocate(Mesh2    (ne2,ne2))
    allocate(Mesh2_map(ne2,ne2))
    allocate(sfcij(0:ne2*ne2,2))
    
    ! create a reverse index array for Mesh2
    ! j = Mesh2(i,j) 
    ! (i,j) = (sfcij(j,1),sfci(j,2)) 
    !------------------------------------------
    call GenspaceCurve(Mesh2)   ! SFC partition for ne2
    do j2=1,ne2
    do i2=1,ne2
      jj=Mesh2(i2,j2)
      sfcij(jj,1)=i2
      sfcij(jj,2)=j2
    end do
    end do
    
    GridVertex(:)%SpaceCurve=-1
    sfc_index   = 0
    do face = 1,nfaces
      ! associate every element on the ne x ne mesh (Mesh)
      ! with its closest element on the ne2 x ne2 mesh (Mesh2)
      ! Store this as a map from Mesh2 -> Mesh in Mesh2_map.
      ! elements in Mesh2 which are not mapped get assigned a value of 0
      !--------------------------------------------------------------------
      Mesh2_map=0
      do ii=1,p_number_elements
        if(face_numbers(ii) == face ) then
          xx = centroids(ii,1)
          yy = centroids(ii,2)
          ! map this element to an (i2,j2) element
          ! [ -DD_PI/4, DD_PI/4 ]  -> [ 0, ne2 ]
          !-------------------------------------------
          i2=nint( (0.5 + 2.0*xx/DD_PI)*ne2 + .5 )
          j2=nint( (0.5 + 2.0*yy/DD_PI)*ne2 + .5 )
          if((face == 4).or.(face == 6 )              ) i2 = ne2-i2+1
          if((face == 1).or.(face == 2).or.(face == 6)) j2 = ne2-j2+1
          if(i2 <   1) i2=1
          if(i2 > ne2) i2=ne2
          if(j2 <   1) j2=1
          if(j2 > ne2) j2=ne2
          Mesh2_map(i2,j2)=ii
        endif
      end do
       
       ! generate a SFC for Mesh with the same ordering as the 
       ! elements in Mesh2 which map to Mesh.
       !----------------------------------------------------------
       do jj=0,ne2*ne2-1
         i2=sfcij(jj,1)
         j2=sfcij(jj,2)
         ii=Mesh2_map(i2,j2)
         if(ii/=0) then
           ! (i2,j2) element maps to element
           GridVertex(ii)%SpaceCurve=sfc_index
           sfc_index=sfc_index+1
         endif
       end do
    end do ! face = 1,nfaces

    deallocate(Mesh2    )
    deallocate(Mesh2_map)
    deallocate(sfcij    )
    
    if(minval(GridVertex(:)%SpaceCurve) == -1) then
      do ii=1,p_number_elements
        if(-1==GridVertex(ii)%SpaceCurve) then
          write (*,*) " Error in projecting element ",ii," to space filling curve."
          write (*,*) " Face:",face_numbers (ii)
          write (*,*) " Centroid:",centroids(ii,:)
        endif
      end do
      call endrun('initialize_space_filling_curve: Vertex not on SpaceCurve')
    endif
    
    ! End Routine
    !---------------
    return
  end subroutine initialize_space_filling_curve
  !======================================================================
  

  !======================================================================
  subroutine find_corner_neighbors(GridVertex, normal_to_homme_ordering, &
                                   element_nodes, corner_wgt, index_table)
    ! find_corner_neighbors:
    !
    !======================================================================
    use err_exit     ,only: endrun
    use gridgraph_mod,only: GridVertex_t
    use SE_Options   ,only: north,south,east,west,neast,seast,nwest,swest
    !
    ! Passed Variables
    !------------------
    type (GridVertex_t),intent(inout):: GridVertex(:)
    integer            ,intent(in)   :: normal_to_homme_ordering(8)
    integer            ,intent(in)   :: element_nodes(p_number_elements, 4)
    integer            ,intent(in)   :: corner_wgt
    integer            ,intent(in)   :: index_table(:,:) 
    !
    ! Local Values
    !--------------
    integer:: node_elements (2*max_elements_attached_to_node)
    integer:: elem_neighbor (4*max_elements_attached_to_node)
    integer:: nbr_cnt(4)
    integer:: elem_nbr_start, start
    integer:: node, loc, cnt, cnt_index
    integer:: corner_array(max_corner_elem), orig_pos(max_corner_elem)
    integer:: face_array(max_corner_elem), a_corner_elems(max_corner_elem)
    integer:: corner_sides(2)
    integer:: side_elem, corner_elem, tmp_s 
    integer:: i0, j0, k0, ll, jj, kk

    !the last column in the index table is a count of the number of elements
    !-----------------------------------------------------------------------
    cnt_index =  max_elements_attached_to_node + 1

    ! loop through all elements
    !------------------------------
    do i0=1, p_number_elements
      node_elements(:) = 0
      elem_neighbor(:) = 0
      nbr_cnt      (:) = 0
      elem_nbr_start   = 0

      ! check each of the 4 nodes at the element corners
      !----------------------------------------------------
      do j0=1,4
        node = element_nodes(i0,j0)
        cnt  = index_table(node, cnt_index)
        if((cnt < 3).or.(max_elements_attached_to_node < cnt)) then
          call endrun('find_corner_neighbors: Number of elements attached to node less than 3 or greater than maximum.')
        endif
        node_elements(1:cnt) = index_table(node, 1:cnt)

        !now node_elements contains the element neighbors to that node - so grab the 
        ! corner neighbors - these are the ones that are not already side neighbors (or myself)
        !---------------------------------------------------------------------------------------
        k0 = 0
        do ll=1,cnt 
          if((           i0          /= node_elements(ll)).and. &
             (GridVertex(i0)%nbrs(1) /= node_elements(ll)).and. &
             (GridVertex(i0)%nbrs(2) /= node_elements(ll)).and. & ! etc ... 
             (GridVertex(i0)%nbrs(3) /= node_elements(ll)).and. &
             (GridVertex(i0)%nbrs(4) /= node_elements(ll))      ) then   
            k0 = k0 + 1
            elem_neighbor(elem_nbr_start + k0) = node_elements(ll)
          endif
        end do ! end of ll loop for multiplicity

        ! keep track of where we are starting in elem_neighbor for each corner j
        !  and how many neighbors in this corner
        !------------------------------------------------------------------------
        elem_nbr_start = elem_nbr_start + k0
        nbr_cnt(j0)    = k0
      end do ! end of j loop through 4 nodes

      !-------------------------------------------------------------------------
      ! now that we have done the 4 corners we can populate nbrs and nbrs_ptr 
      ! with the corners in the proper order (clockwise) in neighbors
      ! also we can add the corner weight
      !-------------------------------------------------------------------------

      ! loop through 4 corners
      !---------------------------
      do j0=5,8
        elem_nbr_start = 1

        !easiest to do the corner in ascending order - find loc
        !-------------------------------------------------------
        do jj = 5,8
          ll = normal_to_homme_ordering(jj)
          if(j0 == ll) then
            loc = jj
            exit
          endif
          elem_nbr_start = elem_nbr_start + nbr_cnt(jj-4)
        end do

        start                         = GridVertex(i0)%nbrs_ptr(j0)
        cnt                           = nbr_cnt(loc - 4)
        GridVertex(i0)%nbrs_ptr(j0+1) = start + cnt
          
        if(cnt > 0) then
          GridVertex(i0)%nbrs_wgt (start:(start+cnt-1)) = corner_wgt
          GridVertex(i0)%nbrs     (start:(start+cnt-1)) =                         &
                               elem_neighbor(elem_nbr_start:(elem_nbr_start+cnt-1))
          GridVertex(i0)%nbrs_face(start:(start+cnt-1)) =                                      &
                    GridVertex(elem_neighbor(elem_nbr_start:(elem_nbr_start+cnt-1)))%face_number
        endif

        ! within each corner neighbor, lets list the corners in clockwise order
        ! cnt is the number of neighbors in this corner j
        ! there can be at most max_corner element of these
        !---------------------------------------------------------------------------
        if(cnt > 1) then
          a_corner_elems = 0
          a_corner_elems = elem_neighbor(elem_nbr_start:(elem_nbr_start+cnt-1))

          !corner-sides(2) is clockwise of corner_side(1)
          !-----------------------------------------------
          corner_array= 0
          orig_pos    = 0
          select case (j0)
            case(neast)
                      corner_sides(1) = north
                      corner_sides(2) = east
            case(seast)
                      corner_sides(1) = east
                      corner_sides(2) = south
            case(swest)
                      corner_sides(1) = south
                      corner_sides(2) = west
            case(nwest)
                      corner_sides(1) = west
                      corner_sides(2) = north
          end select
             
          ! so the first element to list touches  corner_sides(1) element
          !---------------------------------------------------------------
          side_elem = GridVertex(i0)%nbrs(corner_sides(1))
             
          !loop though the corner elements(cnt) and see if any have a 
          ! side neighbor that = side_elem
          !------------------------------------------------------------
          do k0 = 1,cnt
            corner_elem = a_corner_elems(k0)

            ! number of sides to check
            !--------------------------
            do kk = 1,4
              loc   = GridVertex(corner_elem)%nbrs_ptr(kk)
              tmp_s = GridVertex(corner_elem)%nbrs    (loc)
              if(tmp_s == side_elem) then
                corner_array(1) = corner_elem
                orig_pos    (1) = k0
                exit
              endif
            end do
            if(corner_array(1)> 0) exit
          end do
          if(corner_array(1)==0) then
            print *, i0, cnt
            call endrun('find_corner_neighbors (1) : mistake finding corner neighbor order')
          endif

          ! if cnt == 2, we are done (we know the order of neighbors)
          !-----------------------------------------------------------
          if(cnt ==2) then
            if(corner_array(1) ==  a_corner_elems(1)) then
              corner_array(2) =  a_corner_elems(2)
              orig_pos    (2) = 2
            else
              corner_array(2) =  a_corner_elems(1)
              orig_pos    (2) = 1
            endif
          else ! cnt = 3 or 4
            ! find which corner element borders corner_sides(2)
            !--------------------------------------------------
            side_elem = GridVertex(i0)%nbrs(corner_sides(2))
            do k0 = 1,cnt
              corner_elem = a_corner_elems(k0)
              do kk = 1,4
                loc   = GridVertex(corner_elem)%nbrs_ptr(kk)
                tmp_s = GridVertex(corner_elem)%nbrs    (loc)
                if(tmp_s == side_elem) then
                  corner_array(4) = corner_elem
                  orig_pos    (4) = k0
                  exit
                endif
              enddo
              if(corner_array(4)> 0) exit
            enddo
            if((corner_array(4)==0).or.(corner_array(4) == corner_array(1))) then
              print *, i0, cnt
              call endrun('find_corner_neighbors (2) : mistake finding corner neighbor order')
            endif
                
            !now if cnt = 3 then we are done
            !--------------------------------
            if(cnt ==3) then
              corner_array(3) = corner_array(4)
              orig_pos(3)     = orig_pos(4) 
                   
              ! find the "middle" element
              !------------------------------
              do k0 = 1,cnt
                if((k0 /= orig_pos(1)).and.(k0 /= orig_pos(3))) then
                  corner_array(2) = a_corner_elems(k0)
                  orig_pos    (2) = k0
                  exit
                endif
              end do
            else  !cnt = 4 
              ! which of the two unassigned elements borders the element in
              ! corner_array(1) => put in corner_array(2)
              !---------------------------------------------------------------
              side_elem = corner_array(1)
              do k0 = 1,cnt
                corner_elem = a_corner_elems(k0)
                if((corner_elem == corner_array(4)).or. &
                   (corner_elem == corner_array(1))     ) then
                  cycle
                else
                  ! check each side
                  !-----------------
                  do kk = 1,4
                    loc   = GridVertex(corner_elem)%nbrs_ptr(kk)
                    tmp_s = GridVertex(corner_elem)%nbrs    (loc)
                    if(tmp_s == side_elem) then
                      corner_array(2) = corner_elem
                      orig_pos    (2) = k0
                      exit
                    endif
                  end do
                endif
                if(corner_array(2)> 0) exit
              end do
              !now put the remaining one in pos 3
              !------------------------------------
              do k0 = 1,cnt
                corner_elem = a_corner_elems(k0)
                if((corner_elem /= corner_array(4)).and. &
                   (corner_elem /= corner_array(2)).and. &
                   (corner_elem /= corner_array(1))      )then
                  corner_array(3) = corner_elem
                  orig_pos    (3) = k0
                  exit
                endif
              end do
            endif ! end of cnt=4
          endif! end of not cnt=2

          ! now re-set the elements in this corner
          !----------------------------------------
          GridVertex(i0)%nbrs(start:(start+cnt-1)) = corner_array(1:cnt)

          ! nbrs_wgt are the same - nothing to do fix neighbors face
          !----------------------------------------------------------
          do k0 = 1,cnt
            face_array(k0) = GridVertex(i0)%nbrs_face(start+orig_pos(k0)-1)
          end do
          GridVertex(i0)%nbrs_face(start:(start+cnt-1)) = face_array(1:cnt)
        endif !end of cnt > 1 loop for corners
      end do !j loop through each corner
    end do ! end of i loop through elements
    
    ! End Routine
    !---------------
    return
  end subroutine find_corner_neighbors
  !======================================================================


  !======================================================================
!PFC  subroutine MeshOpen(mesh_file_name, par) 
  subroutine MeshOpen(mesh_file_name) 
    ! MeshOpen:
    !
    !======================================================================
!PFC    use parallel_mod,only: parallel_t
    use err_exit    ,only: endrun, iulog
    use SE_Constants,only: real_kind
    !
    ! Passed Variables
    !------------------
    character (len=*),intent(in):: mesh_file_name
!PFC    type (parallel_t),intent(in):: par
    !
    ! Local Values
    !------------------
    integer            :: ncid
    integer,allocatable:: node_multiplicity(:)
    integer            :: kk

    p_mesh_file_name = mesh_file_name
    call open_mesh_file()
   
    p_number_elements   = get_number_of_elements      ()
    p_number_nodes      = get_number_of_nodes         ()
    p_number_blocks     = get_number_of_element_blocks()
    p_number_dimensions = get_number_of_dimensions    ()

    if(p_number_dimensions /= 3) then
      call endrun('The number of dimensions must be 3, otherwise the mesh algorithms will not work')
    endif

    ! Only spheres are allowed in input files.
    !------------------------------------------
!PFC    if(par%masterproc) then
      if(p_number_blocks == 1) then
        write(iulog,*) "Since the mesh file has only one block, it is assumed to be a sphere."
      endif
!PFC    endif
    if(p_number_blocks /= 1) then
      call endrun('Number of elements blocks not exactly 1 (sphere)')
    endif

    ! Because all elements are in one face, this value must match  p_number_elements
    !--------------------------------------------------------------------------------
    p_number_elements_per_face = get_number_of_elements_per_face()

    if(p_number_elements /= p_number_elements_per_face) then
      call endrun('The value of the total number of elements does not match all the elements found in face 1')
    endif

    allocate(p_connectivity(4,p_number_elements_per_face) )
    p_connectivity(:,:)=0

    ! extract the connectivity from the netcdf file
    !-----------------------------------------------
    call get_face_connectivity()
    
    allocate(node_multiplicity(p_number_nodes))
    call get_node_multiplicity(node_multiplicity) 
   
    ! tricky:  For each node with multiplicity n, there are n(n-1) neighbor links
    ! created.  But this counts each edge twice, so:  n(n-1) -n
    ! Should be the same as SUM(SIZE(GridVertex(i)%nbrs(j)%n),i=1:p_number_elements,j=1:8)
    ! p_number_neighbor_edges = dot_product(mult,mult) - 2*sum(mult)
    !----------------------------------------------------------------------------------------
    p_number_neighbor_edges = 0
    do kk=1,p_number_nodes
      p_number_neighbor_edges = p_number_neighbor_edges  &
                              + node_multiplicity(kk)*(node_multiplicity(kk)-2)
    end do
    deallocate(node_multiplicity)
    
    ! allocate the space for the coordinates, this is used in many functions
    !-------------------------------------------------------------------------
    allocate(p_node_coordinates(p_number_nodes, p_number_dimensions))
    call get_node_coordinates()

    if(p_number_elements_per_face /= p_number_elements) then
      call endrun('MeshOpen: Total number of elements not equal to the number of elements on face 1!')
    endif
    
    ! End Routine
    !---------------
    return
  end subroutine MeshOpen
  !======================================================================


  !======================================================================
  subroutine MeshClose
    ! MeshClose:
    !
    ! This routine acts as a destructor cleaning the memory allocated in MeshOpen 
    ! which acts as a constructor allocated dynamical memory for the nodes coordinates.
    !======================================================================

    ! release memory
    !----------------
    deallocate(p_node_coordinates)
    deallocate(p_connectivity)

    ! let the file go
    !----------------
    call close_mesh_file ()
    
    ! End Routine
    !---------------
    return
  end subroutine MeshClose
  !======================================================================


  !======================================================================
!PFC  subroutine MeshPrint(par)
  subroutine MeshPrint
    ! MeshPrint:
    !
    !======================================================================
!PFC    use parallel_mod, only : parallel_t
    use err_exit,only: endrun
    !
    ! Passed Variables
    !------------------
!PFC    type (parallel_t),intent(in):: par

!PFC    if(par%masterproc) then
      print *, 'This are the values for file ', trim(p_mesh_file_name)
      print *, 'The value for the number of dimensions (num_dim) is ', p_number_dimensions
      print *, 'The number of elements in the mesh file is ', p_number_elements
      print *, 'The number of nodes in the mesh file is ', p_number_nodes
      print *, 'The number of blocks in the mesh file is ',  p_number_blocks
      print *, 'The number of elements in the face 1 (sphere) is ',  p_number_elements_per_face
      if(p_number_elements == p_number_elements) then
        print *, 'The value of the total number of elements does match all the elements found in face 1 (the only face)' 
      else
        print *, 'The value of the total number of elements does not match all the elements found in face 1'
        print *, 'This message should not be appearing, there is something wrong in the code'
      endif
      print *, 'The number of neighbor edges ', p_number_neighbor_edges
      !print *, 'The node connectivity are (compare with ncdump -v connect1) ', p_connectivity
      !print *, ' ========================================================='
      !print *, 'The node coordinates are (compare with ncdump -v coord) ', p_node_coordinates
!PFC    endif
    
    ! End Routine
    !---------------
    return
  end subroutine MeshPrint
  !======================================================================


  !======================================================================
  subroutine MeshCubeTopology(GridEdge, GridVertex)
    ! MeshCubeTopology:
    !
    !======================================================================
    use err_exit              ,only: endrun
    use SE_Constants          ,only: np
    use SE_Options            ,only: north,south,east,west,neast,seast,swest,nwest
    use coordinate_systems_mod,only: cartesian3D_t,cube_face_number_from_cart
    use coordinate_systems_mod,only: cube_face_number_from_sphere
    use gridgraph_mod         ,only: GridVertex_t
    use gridgraph_mod         ,only: GridEdge_t
    use gridgraph_mod         ,only: initgridedge, num_neighbors
    use cube_mod              ,only: CubeSetupEdgeIndex
    !
    ! Passed Variables
    !--------------------
    type (GridEdge_t),  intent(inout):: GridEdge(:)
    type (GridVertex_t),intent(inout):: GridVertex(:)
    !
    ! Local Values
    !---------------
    real(kind=real_kind):: coordinates(4,3) 
    real(kind=real_kind):: centroid(3)
    type (cartesian3D_t):: face_center
    integer             :: element_nodes(p_number_elements, 4)
    integer             :: EdgeWgtP,CornerWgt
    integer             :: normal_to_homme_ordering(8)
    integer             :: node_numbers(4)
    integer, allocatable:: index_table(:,:) 

    integer ii,jj,kk,ll,mm,loc

    normal_to_homme_ordering(1) = south
    normal_to_homme_ordering(2) =  east
    normal_to_homme_ordering(3) = north 
    normal_to_homme_ordering(4) =  west
    normal_to_homme_ordering(5) = swest
    normal_to_homme_ordering(6) = seast
    normal_to_homme_ordering(7) = neast
    normal_to_homme_ordering(8) = nwest

    if(SIZE(GridVertex) /= p_number_elements) then
      call endrun('MeshCubeTopology: Element count check failed in exodus_mesh. &
                    &Vertex array length not equal to number of elements.')
    endif
    if(p_number_elements_per_face /= p_number_elements) then
      call endrun('MeshCubeTopology: Element count check failed in exodus_mesh. &
                    &Element array length not equal to sum of face.')
    endif

    EdgeWgtP = np
    CornerWgt= 1

    call mesh_connectivity (element_nodes)

    do ii=1, p_number_elements  
      GridVertex(ii)%number           = ii
      GridVertex(ii)%face_number      = 0
      GridVertex(ii)%processor_number = 0
      GridVertex(ii)%SpaceCurve       = 0
      GridVertex(ii)%nbrs          (:)= 0
      GridVertex(ii)%nbrs_face     (:)= 0
      GridVertex(ii)%nbrs_wgt      (:)= 0
      GridVertex(ii)%nbrs_wgt_ghost(:)= 1

      ! each elements has one side neighbor (first 4)
      !-----------------------------------------------
      GridVertex(ii)%nbrs_ptr(1) = 1
      GridVertex(ii)%nbrs_ptr(2) = 2
      GridVertex(ii)%nbrs_ptr(3) = 3
      GridVertex(ii)%nbrs_ptr(4) = 4

      ! don't know about corners yet
      !--------------------------------
      GridVertex(ii)%nbrs_ptr(5:num_neighbors+1) = 5
    end do ! ii=1, p_number_elements  

    ! create index table to find neighbors
    !---------------------------------------
    call create_index_table(index_table, element_nodes)

    ! side neighbors 
    !---------------
    call find_side_neighbors(GridVertex, normal_to_homme_ordering, &
                             element_nodes, EdgeWgtP, index_table)
   
    ! set vertex faces
    !-------------------
    do ii=1, p_number_elements
      node_numbers  = element_nodes     (ii,:)    
      coordinates   = p_node_coordinates(node_numbers,:)
      centroid      = SUM(coordinates, dim=1)/4.0
      face_center%x = centroid(1)
      face_center%y = centroid(2)
      face_center%z = centroid(3)
      GridVertex(ii)%face_number = cube_face_number_from_cart(face_center)
    end do

    ! set side neighbor faces
    !-------------------------
    do ii=1, p_number_elements
    do jj=1,4 ! look at each side
      kk  = normal_to_homme_ordering(jj)
      loc = GridVertex(ii)%nbrs_ptr(kk)
      ll  = GridVertex(ii)%nbrs    (loc)
      GridVertex(ii)%nbrs_face(loc) = GridVertex(ll)%face_number
    end do
    end do

    ! find corner neighbor and faces (weights added also)
    !------------------------------------------------------
    call find_corner_neighbors(GridVertex, normal_to_homme_ordering, &
                               element_nodes, CornerWgt, index_table)

    ! done with the index table
    !--------------------------
    deallocate(index_table)
  
    call initgridedge(GridEdge,GridVertex) 
    do ii=1,SIZE(GridEdge)
      call CubeSetupEdgeIndex(GridEdge(ii)) 
    end do

    call initialize_space_filling_curve(GridVertex, element_nodes)
    
    ! End Routine
    !---------------
    return
  end subroutine MeshCubeTopology
  !======================================================================


  !======================================================================
  subroutine MeshSetCoordinates(elem)
    ! MeshSetCoordinates:
    !
    !======================================================================
    use element_mod           ,only: element_t
    use err_exit              ,only: endrun
    use coordinate_systems_mod,only: cartesian3D_t, cartesian2d_t, spherical_polar_t, &
                                     change_coordinates, sphere2cubedsphere
    !
    ! Passed Variables
    !--------------------
    type(element_t),intent(inout):: elem(:)
    !
    ! Local Values
    !--------------
    integer             :: connectivity(p_number_elements,4)
    integer             :: node_multiplicity(p_number_nodes)
    integer             :: number
    integer             :: node_num(4)
    real(kind=real_kind):: coordinates(4,3)
    real(kind=real_kind):: cube_coor  (4,2)
    real(kind=real_kind):: x_double          
    real                :: x_real            
    type(cartesian2d_t) :: cart2

    integer face_no, ii, kk, ll
    
    connectivity     = 0
    node_multiplicity= 0
    call mesh_connectivity (connectivity)

    do kk=1,p_number_elements
      node_num                       = connectivity(kk,:)
      node_multiplicity(node_num(:)) = node_multiplicity(node_num(:)) + 1
    end do

    do kk=1,SIZE(elem) 
      number      = elem(kk)%vertex%number
      face_no     = elem(kk)%vertex%face_number
      node_num    = connectivity(number,:)
      coordinates = p_node_coordinates(node_num,:)
      if(6 == p_number_blocks) then
        call cube_to_cube_coordinates  (cube_coor, coordinates, face_no)
      else
        call sphere_to_cube_coordinates(cube_coor, coordinates, face_no)
      endif
      elem(kk)%node_numbers         = node_num 
      elem(kk)%node_multiplicity(:) = node_multiplicity(node_num(:))
      elem(kk)%corners(:)%x         = cube_coor(:,1)
      elem(kk)%corners(:)%y         = cube_coor(:,2)
    end do
    
    ! End Routine
    !---------------
    return
  end subroutine MeshSetCoordinates
  !======================================================================


  !======================================================================
  function MeshCubeEdgeCount() result(nedge)
    ! MeshCubeEdgeCount:
    !
    !======================================================================
    use err_exit,only: endrun
    !
    ! Passed Variables
    !-------------------
    integer:: nedge
   
    if(0 == p_number_blocks)  call endrun('MeshCubeEdgeCount called before MeshOpenMesh')
    if(MeshUseMeshFile) then
      ! should be the same as SUM(SIZE(GridVertex(i)%nbrs(j)%n),i=1:p_number_elements,j=1:nInnerElemEdge)
      ! the total number of neighbors.
      !----------------------------------------------------------------------------------
      nedge = p_number_neighbor_edges
    else
      call endrun('Error in MeshCubeEdgeCount: Should not call for non-exodus mesh file.')
    endif

    ! End Function
    !---------------
    return
  end function MeshCubeEdgeCount
  !======================================================================


  !======================================================================
  function MeshCubeElemCount()  result(nelem)
    ! MeshCubeElemCount:
    !
    !======================================================================
    use err_exit,only: endrun
    !
    ! Passed Variables
    !--------------------
    integer:: nelem

    if(0 == p_number_blocks)  call endrun('MeshCubeElemCount called before MeshOpenMesh')
    if(MeshUseMeshFile) then
      nelem = p_number_elements
    else
      call endrun('Error in MeshCubeElemCount: Should not call for non-exodus mesh file.')
    end if

    ! End Function
    !---------------
    return
  end function MeshCubeElemCount
  !======================================================================


  !======================================================================
  subroutine test_private_methods
    ! test_private_methods:
    !
    !======================================================================
    !
    ! Local Values
    !-----------------
    integer:: element_nodes(p_number_elements, 4)

    call mesh_connectivity (element_nodes)
    
    ! End Routine
    !---------------
    return
  end subroutine test_private_methods
  !======================================================================

end module mesh_mod



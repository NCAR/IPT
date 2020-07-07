module SE_Element_mod
!===================================================================================
! 
! Purpose: Provide an datastructure (SEelem_t) for cubed-sphere Spectral 
!          Element (SE) gridded data and a routine to initialize its values
!          for the given options.
!
! Author: Patrick Callaghan
!
! Description:
!     This module initializes the SEelem_t datastructure using the step by 
!     step process used in the 'cesm1_3_beta06' version of the model. The 
!     subset of needed values to support SE computations is extracted and 
!     stored in trimmed down datastructure. Modules extracted from the 
!     SE dycore were only modified to remove dependence on MPI and unrelated 
!     values used to solve the primitive equations. They were also modified 
!     to try and introduce some form of sanity to the formatting. The use of 
!     these modules is isolated to the Init_BASEelem() routine, which follows 
!     the step by step calls in the initialization of CAM-SE. The base element
!     values are augmented with additional information needed to implement
!     transformations between the Sphere and the reference element 
!     (e.g. transformation of 2nd rank tensor, Christoffel symbols, etc..)
!
! VERSION: 1.0 
!
! Revisions: 
!   Feb 2015 - Original Version
!
! Usage: 
!    Initialization: 
!    --------------- 
!      Init_SEelem(SEopt,SEelem,SEcv)
!                 Given the grid options specifying resolution and related
!                 parameters (SEopt), return the iniitalized spectral element 
!                 grid datastructure(SEelem). If the optional control volume
!                 datastructure is present, then the control volume values
!                 are calculated and returned as well.
!
!      Init_BASEelem(SEopt,elem)  :  PRIVATE routine
!                 This routine follows the cesm1_3_beta06 revision of the SE 
!                 dycore to initialize the elem data structure from that model.
!                 After all of the init code is complete, the values needed 
!                 are copied to initialize SEelem which contain only the essential 
!                 components needed.
!
!    Grid Validation:
!    ----------------
!      SE_Get_Ncol(ncol) 
!                 Routine to return the number of SE grid columns (ncol) after
!                 initialization.
!        
!      SE_Get_Nelem(nelem) 
!                 Routine to return the number of SE elements (nelem) after
!                 initialization.
!        
!      SE_Get_Face(SEelem,FaceNumber)
!                 Given the SEelem grid datastructure, return an (np,np,nelem)
!                 array in which each point is set the face number it is located
!                 on after initialization.
!        
!      SE_Get_Number(SEelem,Number)
!                 Given the SEelem grid datastructure, return an (np,np,nelem)
!                 array in which each point is set the element number it is 
!                 located in after initialization.
!
! TODO:
!-------
!
!===================================================================================
  ! Useful modules
  !----------------
  use err_exit       ,only: iulog,endrun
  use quadrature_mod ,only: quadrature_t,gausslobatto
  use SE_Constants   ,only: real_kind,log_kind,int_kind,nlev,         &
                            np,nc,npsq,DD_PI,MAX_FILE_LEN,rrearth
  use SE_Options     ,only: south,east,north,west,swest,SFCURVE

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private
  save

  public :: SEedge_Desc_t
  private:: SEindex_t
  private:: SEspherical_t
  public :: SEelem_t

  public :: Init_SEelem
  private:: Init_BASEelem
  public :: SE_Get_Ncol
  public :: SE_Get_Nelem
  public :: SE_Get_Face
  public :: SE_Get_Number

  private:: SEinitialized
  private:: SEncol
  private:: SEnp
  private:: SEnlev
  private:: SEnelem

  ! Type Definitions
  !-------------------
  type SEedge_Desc_t
    integer                       :: actual_neigh_edges
    integer(kind=int_kind),pointer:: putmapP(:)       => null()
    integer(kind=int_kind),pointer:: getmapP(:)       => null()
    integer(kind=int_kind),pointer:: globalID(:)      => null()
    logical(kind=log_kind),pointer:: reverse(:)       => null()
  end type SEedge_Desc_t

  type SEindex_t
    integer(kind=int_kind):: ia(npsq),ja(npsq)
    integer(kind=int_kind):: is,ie
    integer(kind=int_kind):: NumUniquePts
    integer(kind=int_kind):: UniquePtOffset
  end type SEindex_t

  type SEcart2D_t
    real(real_kind):: x             ! x coordinate
    real(real_kind):: y             ! y coordinate
  end type SEcart2D_t

  type SEcart3D_t
    real(real_kind):: x             ! x coordinate
    real(real_kind):: y             ! y coordinate
    real(real_kind):: z             ! z coordinate
  end type SEcart3D_t

  type SEspherical_t
    real(real_kind):: r             ! radius
    real(real_kind):: lon           ! longitude
    real(real_kind):: lat           ! latitude
  end type SEspherical_t

  type SEelem_t
    ! Base elemet values
    !--------------------
    integer             :: LocalId
    integer             :: GlobalId
    integer             :: face_number
    integer             :: number
    type(SEspherical_t) :: spherep  (np,np)
    type(SEcart2D_t   ) :: cartp    (np,np) ! gnomonic coords of GLL points 
    type(SEcart2D_t   ) :: corners  (4)         ! gnomonic coords of element corners
    type(SEcart3D_t   ) :: corners3D(4)
    real(kind=real_kind):: spheremp (np,np) ! mass matrix on v and p grid
    real(kind=real_kind):: rspheremp(np,np) ! inverse mass matrix on v and p grid
    real(kind=real_kind):: metinv(2,2,np,np)! metric tensor on velocity and pressure grid
    real(kind=real_kind):: metdet   (np,np) ! g=SQRT(det(g_ij)) on velocity/pressure grid
    real(kind=real_kind):: rmetdet  (np,np) ! 1/metdet on velocity/pressure grid
    real(kind=real_kind):: mp       (np,np) ! mass matrix on v and p grid 
    real(kind=real_kind):: fcor     (np,np) ! Coriolis factor
    real(kind=real_kind):: Dinv (2,2,np,np) ! Map vector field on sphere to covariant v on cube
    real(kind=real_kind):: D    (2,2,np,np) ! Map covariant v on cube to vector field on sphere
    real(kind=real_kind):: variable_hyperviscosity(np,np) ! hyperviscosity based on above
    real(kind=real_kind):: tensorVisc(2,2,np,np)          !og, matrix V for tensor viscosity
    real(kind=real_kind):: vec_sphere2cart(np,np,3,2)     ! Conv vectors spherical-rectangular
    type(SEindex_t)     :: idxP
    type(SEedge_Desc_t) :: desc
  end type SEelem_t

  ! Global Data
  !-------------
  logical                          :: SEinitialized = .false.
  integer                          :: SEncol
  integer                          :: SEnp
  integer                          :: SEnlev
  integer                          :: SEnelem
  character(len=MAX_FILE_LEN)      :: SEmesh_file

contains

  !=============================================================================
  subroutine Init_SEelem(IO_SEopt,O_SEelem,O_SEcv)
    !
    ! Init_SEelem: Given the grid options specifying resolution and related
    !              parameters, the routine BASE elem data structure found in 
    !              the cesm1_3_beta06 revision of the SE dycore is initialized. 
    !
    !              After all of the init code is complete, the values needed 
    !              are copied to initialize SEelem  which contain
    !              only the essential components needed for the differential 
    !              operators.
    !===========================================================================
    use SE_Options ,only: SEoptions_t
    use element_mod,only: element_t
    use SE_ControlVolume_mod,only: SEctrlvol_t, Calc_ControlVolumes
    ! 
    ! Passed Variables
    !------------------
    type(SEoptions_t)            ,intent(inout):: IO_SEopt
    type(SEelem_t   ),allocatable,intent(out  ):: O_SEelem(:)
    type(SEctrlvol_t),optional   ,intent(out  ):: O_SEcv
    !
    ! Local Values
    !--------------
    type(element_t),target,allocatable:: elem(:)
    integer:: ie,ii,jj
    
    ! Initialize elem() data structure from CAM-SE base code
    !----------------------------------------------------------
    call Init_BASEelem(IO_SEopt,elem)


    !===========================================================
    ! All CAM-SE initializations are now complete. Allocate and
    ! set global values needed for the stand-alone SE module.
    !===========================================================
    SEncol = IO_SEopt%GlobalUniqueCols
    SEnp   = np
    SEnlev = nlev
    SEnelem= IO_SEopt%nelem
    write(iulog,*) ' SEncol =',SEncol
    write(iulog,*) ' SEnp   =',SEnp
    write(iulog,*) ' SEnlev =',SEnlev
    write(iulog,*) ' SEnelem=',SEnelem

    allocate(O_SEelem(IO_SEopt%nelem))
    do ie=1,IO_SEopt%nelem
      allocate(O_SEelem(ie)%desc%putmapP      (IO_SEopt%max_neigh_edges))
      allocate(O_SEelem(ie)%desc%getmapP      (IO_SEopt%max_neigh_edges))
!XX      allocate(O_SEelem(ie)%desc%putmapP_ghost(IO_SEopt%max_neigh_edges))
!XX      allocate(O_SEelem(ie)%desc%getmapP_ghost(IO_SEopt%max_neigh_edges))
      allocate(O_SEelem(ie)%desc%reverse      (IO_SEopt%max_neigh_edges))
      allocate(O_SEelem(ie)%desc%globalID     (IO_SEopt%max_neigh_edges))
!XX      allocate(O_SEelem(ie)%desc%loc2buf      (IO_SEopt%max_neigh_edges))
      do ii=1,IO_SEopt%max_neigh_edges
!XX        O_SEelem(ie)%desc%loc2buf (ii)=ii
        O_SEelem(ie)%desc%globalID(ii)=-1
      end do
    end do
    write(iulog,*) ' O_SEelem allocated'

    do ie=1,IO_SEopt%nelem
      O_SEelem(ie)%spheremp (:,:)          = elem(ie)%spheremp (:,:)
      O_SEelem(ie)%rspheremp(:,:)          = elem(ie)%rspheremp(:,:)
      O_SEelem(ie)%metinv(:,:,:,:)         = elem(ie)%metinv(:,:,:,:)
      O_SEelem(ie)%metdet   (:,:)          = elem(ie)%metdet   (:,:)
      O_SEelem(ie)%rmetdet  (:,:)          = elem(ie)%rmetdet  (:,:)
      O_SEelem(ie)%mp       (:,:)          = elem(ie)%mp       (:,:)
      O_SEelem(ie)%spherep  (:,:)%r        = elem(ie)%spherep  (:,:)%r
      O_SEelem(ie)%spherep  (:,:)%lat      = elem(ie)%spherep  (:,:)%lat
      O_SEelem(ie)%spherep  (:,:)%lon      = elem(ie)%spherep  (:,:)%lon
      O_SEelem(ie)%cartp    (:,:)%x        = elem(ie)%cartp    (:,:)%x
      O_SEelem(ie)%cartp    (:,:)%y        = elem(ie)%cartp    (:,:)%y
      O_SEelem(ie)%corners  (:)%x          = elem(ie)%corners  (:)%x
      O_SEelem(ie)%corners  (:)%y          = elem(ie)%corners  (:)%y
      O_SEelem(ie)%corners3D(:)%x          = elem(ie)%corners3D(:)%x
      O_SEelem(ie)%corners3D(:)%y          = elem(ie)%corners3D(:)%y
      O_SEelem(ie)%corners3D(:)%z          = elem(ie)%corners3D(:)%z
      O_SEelem(ie)%fcor     (:,:)          = elem(ie)%fcor     (:,:)
      O_SEelem(ie)%Dinv     (:,:,:,:)      = elem(ie)%Dinv     (:,:,:,:)
      O_SEelem(ie)%D        (:,:,:,:)      = elem(ie)%D        (:,:,:,:)
      O_SEelem(ie)%vec_sphere2cart(:,:,:,:)= elem(ie)%vec_sphere2cart(:,:,:,:)
      O_SEelem(ie)%variable_hyperviscosity(:,:)= elem(ie)%variable_hyperviscosity(:,:)
      O_SEelem(ie)%tensorVisc(:,:,:,:)         = elem(ie)%tensorVisc(:,:,:,:)         
!XX      O_SEelem(ie)%desc%use_rotation      = elem(ie)%desc%use_rotation
!XX      O_SEelem(ie)%desc%padding           = elem(ie)%desc%padding
      O_SEelem(ie)%desc%actual_neigh_edges= elem(ie)%desc%actual_neigh_edges
      O_SEelem(ie)%desc%putmapP      (:)  = elem(ie)%desc%putmapP      (:)
      O_SEelem(ie)%desc%getmapP      (:)  = elem(ie)%desc%getmapP      (:)
!XX      O_SEelem(ie)%desc%putmapP_ghost(:)  = elem(ie)%desc%putmapP_ghost(:)
!XX      O_SEelem(ie)%desc%getmapP_ghost(:)  = elem(ie)%desc%getmapP_ghost(:)
      O_SEelem(ie)%desc%reverse      (:)  = elem(ie)%desc%reverse      (:)
      O_SEelem(ie)%desc%globalID     (:)  = elem(ie)%desc%globalID     (:)
!XX      O_SEelem(ie)%desc%loc2buf      (:)  = elem(ie)%desc%loc2buf      (:)
!XX      if(associated(elem(ie)%desc%rot)) then
!XX        nrot = SIZE(elem(ie)%desc%rot)
!XX        allocate(O_SEelem(ie)%desc%rot(nrot))
!XX        do ii=1,nrot
!XX          O_SEelem(ie)%desc%rot(ii)%nbr        = elem(ie)%desc%rot(ii)%nbr
!XX          O_SEelem(ie)%desc%rot(ii)%reverse    = elem(ie)%desc%rot(ii)%reverse
!XX          npt = SIZE(elem(ie)%desc%rot(ii)%R,dim=3)
!XX          allocate(O_SEelem(ie)%desc%rot(ii)%R(2,2,npt))
!XX          O_SEelem(ie)%desc%rot(ii)%R(:,:,:)   = elem(ie)%desc%rot(ii)%R(:,:,:)
!XX        end do
!XX      endif
      O_SEelem(ie)%idxP%ia(:)          = elem(ie)%idxP%ia(:)
      O_SEelem(ie)%idxP%ja(:)          = elem(ie)%idxP%ja(:)
      O_SEelem(ie)%idxP%is             = elem(ie)%idxP%is
      O_SEelem(ie)%idxP%ie             = elem(ie)%idxP%ie
      O_SEelem(ie)%idxP%NumUniquePts   = elem(ie)%idxP%NumUniquePts
      O_SEelem(ie)%idxP%UniquePtOffset = elem(ie)%idxP%UniquePtOffset
      O_SEelem(ie)%face_number         = elem(ie)%vertex%face_number
      O_SEelem(ie)%number              = elem(ie)%vertex%number
      O_SEelem(ie)%LocalId             = elem(ie)%LocalId
      O_SEelem(ie)%GlobalId            = elem(ie)%GlobalId
    end do

    ! Initialization is done, set the flag
    !-------------------------------------
    write(iulog,*) ' O_SEelem initialized'
    SEinitialized = .true.

    ! Calculate Control volume information for current grid
    !-------------------------------------------------------
    if(present(O_SEcv)) then
      call Calc_ControlVolumes(elem,IO_SEopt,O_SEcv)
    endif

    ! Clean up: no need for elem any more.
    !---------------------------------------
    do ie=1,IO_SEopt%nelem
      deallocate(elem(ie)%desc%putmapP      )
      deallocate(elem(ie)%desc%getmapP      )
      deallocate(elem(ie)%desc%putmapP_ghost)
      deallocate(elem(ie)%desc%getmapP_ghost)
      deallocate(elem(ie)%desc%reverse      )
      deallocate(elem(ie)%desc%globalID     )
      deallocate(elem(ie)%desc%loc2buf      )
    end do
    deallocate(elem)

    ! End Routine
    !--------------
    return
  end subroutine Init_SEelem
  !==================================================================


  !=============================================================================
  subroutine Init_BASEelem(IO_SEopt,elem)
    !
    ! Init_BASEelem: Given the grid options specifying resolution and related
    !                parameters, the routine goes thru the initialization step
    !                for the elem/deriv datastructures found in the cesm1_3_beta06
    !                revision of the SE dycore. In general, the original coding
    !                is followed as closely as possible. Ommisions focues on 
    !                initializations of unneeded values related to the primitive
    !                equations and all MPI realted coding. 
    !                After all of the init code is complete, the elem() data 
    !                structure is returned.
    !===========================================================================
    use SE_Constants   ,only: real_kind,longdouble_kind,rearth,rrearth, &
                              DD_PI,np,nc,nep
    use SE_Options     ,only: SEoptions_t
    use quadrature_mod ,only: quadrature_t,gausslobatto, &
                              test_gausslobatto,test_gauss
    use element_mod    ,only: element_t,allocate_element_desc,init_element_mod
    use gridgraph_mod  ,only: GridVertex_t,GridEdge_t,                             &
                              allocate_gridvertex_nbrs,deallocate_gridvertex_nbrs, &
                              init_gridgraph_mod
    use cube_mod       ,only: CubeElemCount,CubeEdgeCount,CubeTopology,           &
                              set_corner_coordinates,assign_node_numbers_to_elem, &
                              cube_init_atomic,rotation_init_atomic,init_cube_mod
    use mesh_mod       ,only: MeshUseMeshFile,MeshCubeElemCount,  &
                              MeshCubeEdgeCount,MeshCubeTopology, &
                              MeshSetCoordinates,MeshOpen,init_mesh_mod
    use spacecurve_mod ,only: genspacepart,init_spacecurve_mod
    use mass_matrix_mod,only: mass_matrix,init_mass_matrix_mod
    use dof_mod        ,only: global_dof,CreateUniqueIndex,SetElemOffset,init_dof_mod
    use derivative_mod ,only: derivative_t,derivinit,init_derivative_mod
    use edge_mod       ,only: initEdgeBuffer,edgeVpack,edgeVunpack,FreeEdgeBuffer, &
                              EdgeBuffer_t,init_edge_mod
    use coordinate_systems_mod,only: spherical_polar_t
    ! 
    ! Passed Variables
    !------------------
    type(SEoptions_t)                   ,intent(inout):: IO_SEopt
    type(element_t  ),target,allocatable,intent(out  ):: elem(:)
    !
    ! Local Values
    !--------------
    integer             jj,ii
    type(GridVertex_t),target,allocatable:: GridVertex(:)
    type(GridEdge_t)  ,target,allocatable:: Gridedge(:)
    type(derivative_t)                   :: deriv
    type(EdgeBuffer_t)                   :: edgebuf
    
    integer,allocatable::edgeptrP      (:)
    integer,allocatable::edgeptrP_ghost(:)
    integer,allocatable::Global2Local  (:)

    integer icount,i0,il,ig
    integer face
    integer loc
    integer direction
    integer nelem_edge
    integer ie
    integer l1,l2,l1id,l2id

    type(quadrature_t)  :: gp
    real(kind=real_kind):: area

    real(kind=real_kind),allocatable:: Htest(:,:,:)
    real(kind=real_kind):: I_sphere,Darea
    real(kind=real_kind):: min_area,max_area,max_ratio
    real(kind=real_kind):: min_max_eig,max_max_eig,min_min_eig,max_min_eig
    real(kind=real_kind):: min_max_dx,max_max_dx,min_min_dx,max_min_dx
    real(kind=real_kind):: avg_area,avg_max_eig,avg_min_eig,avg_max_dx,avg_min_dx

    real(kind=longdouble_kind):: fvm_corners(nc+1)  ! fvm cell corners on reference element
    real(kind=longdouble_kind):: fvm_points(nc)     ! fvm cell centers on reference element
    real(kind=longdouble_kind):: spelt_refnep(1:nep)
    real(kind=real_kind)      :: xtmp

    real(kind=real_kind):: lambda_max,lambda_vis,lambda
    real(kind=real_kind):: min_gw,stable_hv
    real(kind=real_kind):: min_hypervis,max_hypervis,avg_hypervis
    real(kind=real_kind):: max_unif_dx,max_eig_hypervis

    real(kind=real_kind):: dtnu   ! timestep*viscosity parameter (for CFL calc)
    real(kind=real_kind):: noreast,nw,se,sw,xx,yy
    real(kind=real_kind),allocatable:: zeta(:,:,:)
    integer                         :: rowind,colind
    integer                         :: nrot,npt

    real*8 XX0,WG1,WG2,WG3,WG4,WG5,WG6,WG7,X0p,X0n,X02,X0s
    integer     :: my_unit,my_number
    character*80:: FileName
     
    ! Test element grid
    !---------------------
    print *,' np=',np
    call test_gauss(np)
    call test_gausslobatto(np)

    ! Set Resolution dependent values
    !----------------------------------
    if(IO_SEopt%ne.gt.0) then
     ! FIXED resolution
     !---------------------
      IO_SEopt%hypervis_scaling= 0
    elseif(IO_SEopt%ne.eq.0) then
      ! VARIABLE mesh
      !-----------------------
      IO_SEopt%hypervis_scaling= 3.2
      SEmesh_file = trim(IO_SEopt%InputPath)//trim(IO_SEopt%MeshFile)
      if(SEmesh_file.eq."NULL") then
        call endrun("Init_BASEelem: For NE=0 a (.g) mesh file must be provided")
      endif
    else
      write(iulog,*) 'Init_BASEelem: ne=',IO_SEopt%ne
      call endrun("Init_BASEgrid: Negative Resolution?")
    endif

    ! Use an exodus mesh file if ne=0
    !-----------------------------------
    if(IO_SEopt%ne.eq.0) then
      write (iulog,*) 'Opening Mesh File:', trim(SEmesh_file)
      IO_SEopt%max_elements_attached_to_node = 8 ! 7                  ! variable resolution
      IO_SEopt%s_nv                          = 2*IO_SEopt%max_elements_attached_to_node
      IO_SEopt%max_corner_elem = IO_SEopt%max_elements_attached_to_node-3
      IO_SEopt%max_neigh_edges = 4 + 4*IO_SEopt%max_corner_elem
      call init_mesh_mod(IO_SEopt)
      call MeshOpen(SEmesh_file)
    else
      call init_mesh_mod(IO_SEopt)
    endif

    ! set cubed_sphere_map
    !-------------------------
    if(IO_SEopt%cubed_sphere_map < 0) then
      if(IO_SEopt%ne .eq. 0) then
        ! element_local for var-res grids
        !---------------------------------
        IO_SEopt%cubed_sphere_map=2
      else
        ! default is equi-angle gnomonic
        !---------------------------------
        IO_SEopt%cubed_sphere_map=0
      endif
    endif
    write (iulog,*) "Reference element projection: cubed_sphere_map=",IO_SEopt%cubed_sphere_map
    call init_cube_mod(IO_SEopt)

    ! Create the cube topology
    !-----------------------------
    if(trim(IO_SEopt%topology) == "cube") then
      print*,' Creating CUBE topology'
    
      ! Get the element counts
      !-------------------------
      if(MeshUseMeshFile) then
        IO_SEopt%nelem = MeshCubeElemCount()
        nelem_edge     = MeshCubeEdgeCount()
      else
        IO_SEopt%nelem = CubeElemCount()
        nelem_edge     = CubeEdgeCount()
      endif
      IO_SEopt%nelemd = IO_SEopt%nelem
      print *,'      nelem=',IO_SEopt%nelem
      print *,' nelem_edge=',nelem_edge

      ! Allocate GridVertex and GridEdge structures
      !---------------------------------------------
      call init_gridgraph_mod(IO_SEopt)
      allocate(GridVertex(IO_SEopt%nelem))
      allocate(GridEdge  (nelem_edge    ))
      do jj =1,IO_SEopt%nelem
        call allocate_gridvertex_nbrs(GridVertex(jj))
      end do
      print *,' Grid allocation done'
      print *,' Grid vertex nbrs allocation done'

      ! Set up Grids
      !---------------------
      if(MeshUseMeshFile) then
        call MeshCubeTopology(GridEdge,GridVertex)
      else
        call CubeTopology(GridEdge,GridVertex)
      endif
      print *,' Grid topology done'

    else
      print*,' Unknown topology=',trim(IO_SEopt%topology)
      stop
    endif

    ! Set up grid partitioning
    !---------------------------
    print *,' partmethod=',IO_SEopt%partmethod
    if(IO_SEopt%partmethod .eq. SFCURVE) then
      write(iulog,*)"partitioning graph using SF Curve..."
      call init_spacecurve_mod(IO_SEopt)
      call genspacepart(GridEdge,GridVertex)
    else
      write(iulog,*)"partitioning graph using Metis..."
!      call genmetispart(GridEdge,GridVertex)
      stop
    endif

    !=====================================================================
    !=====================================================================
 
    ! Allocate some work space
    !---------------------------
    allocate(edgeptrP      (nelem_edge))
    allocate(edgeptrP_ghost(nelem_edge))
    allocate(Global2Local(IO_SEopt%nelem))
    edgeptrP      (:) = 0
    edgeptrP_ghost(:) = 0
    Global2Local  (:) = -1

    ! Initialize Pointers for edge buffer storage
    !----------------------------------------------
    icount = 1
    do i0=1,nelem_edge
      if( (icount+1) .le. nelem_edge ) then
        ii = GridEdge(i0)%tail_face
        edgeptrP(icount+1) = edgeptrP(icount) + GridEdge(i0)%tail%nbrs_wgt(ii)
        edgeptrP_ghost(icount+1) = edgeptrP_ghost(icount)             &
                                 + GridEdge(i0)%tail%nbrs_wgt_ghost(ii)
        icount = icount + 1
      endif
    end do

    ! Allocate element datastructures
    !----------------------------------
    call init_element_mod(IO_SEopt)
    allocate(elem(IO_SEopt%nelem))
    call allocate_element_desc(elem)
    print *,' nelem=',IO_SEopt%nelem
    print *,' elem allocated'
    print *,' elem_desc allocated'

    ! Initialize edge pointers for each element
    !--------------------------------------------
    do il=1,IO_SEopt%nelem
      ig               = GridVertex(il)%number
      Global2Local(ig) = il
      elem(il)%desc%putmapP      (:) = -1
      elem(il)%desc%getmapP      (:) = -1
      elem(il)%desc%putmapP_ghost(:) = -1
      elem(il)%desc%getmapP_ghost(:) = -1
      elem(il)%desc%reverse      (:) = .FALSE.
    end do

    ! Set edge pointers for sending/receiving 
    ! edge data between neighbors
    !--------------------------------------------
    do ii=1,nelem_edge
      ! Send.
      !--------
      il   = Global2Local(GridEdge(ii)%tail%number)
      face =              GridEdge(ii)%tail_face
      if(face.ge.5) then
        direction = GridEdge(ii)%tail_dir
        loc = MOD(direction,IO_SEopt%max_corner_elem)
        direction = (direction-loc)/IO_SEopt%max_corner_elem
        loc = direction + (direction-5)*(IO_SEopt%max_corner_elem-1) + loc
      else
        loc = face
      endif
      if(il.gt.0) then
        elem(il)%desc%reverse      (loc) = GridEdge      (ii)%reverse
        elem(il)%desc%putmapP      (loc) = edgeptrP      (ii) + 0
        elem(il)%desc%putmapP_ghost(loc) = edgeptrP_ghost(ii) + 1
      endif
      ! Receive
      !---------
      il   = Global2Local(GridEdge(ii)%head%number)
      face =              GridEdge(ii)%head_face
      if(face.ge.5) then
        direction = GridEdge(ii)%head_dir
        loc = MOD(direction,IO_SEopt%max_corner_elem)
        direction = (direction-loc)/IO_SEopt%max_corner_elem
        loc = direction + (direction-5)*(IO_SEopt%max_corner_elem-1) + loc
      else
        loc = face
      endif
      if(il.gt.0) then
        elem(il)%desc%globalID     (loc) = GridEdge      (ii)%tail%number
        elem(il)%desc%getmapP      (loc) = edgeptrP      (ii)  + 0
        elem(il)%desc%getmapP_ghost(loc) = edgeptrP_ghost(ii) + 1
      endif
    end do ! ii=1,nelem_edge 

    ! Done with the work space
    !---------------------------
    deallocate(edgeptrP      )
    deallocate(edgeptrP_ghost)
    deallocate(Global2Local  )

    ! Finish setting edge and vertex 
    ! information for each element
    !----------------------------------------
    do ie=1,IO_SEopt%nelem
      ! Count the actual number of edges
      !-----------------------------------
      elem(ie)%desc%actual_neigh_edges = 0
      do ii=1,IO_SEopt%max_neigh_edges
        if(elem(ie)%desc%globalID(ii) > 0) then
          elem(ie)%desc%actual_neigh_edges = elem(ie)%desc%actual_neigh_edges + 1
        endif
      end do
      ! Sort edges
      !-----------------------------------
      do l1 =(   1),(IO_SEopt%max_neigh_edges-1)
      do l2 =(l1+1),(IO_SEopt%max_neigh_edges  )
        l1id = elem(ie)%desc%loc2buf(l1)
        l2id = elem(ie)%desc%loc2buf(l2)
        if(elem(ie)%desc%globalID(l2id) > elem(ie)%desc%globalID(l1id)) then
          l1id = elem(ie)%desc%loc2buf(l2)
          elem(ie)%desc%loc2buf(l2) = elem(ie)%desc%loc2buf(l1)
          elem(ie)%desc%loc2buf(l1) = l1id
        endif
      end do
      end do
      elem(ie)%vertex   = GridVertex(ie)
      elem(ie)%GlobalId = GridVertex(ie)%number
      elem(ie)%LocalId  = ie
    end do  !ie=1,nelem
    !=====================================================================
    !=====================================================================

    ! Initialize cube element coordinates
    !-------------------------------------
    gp = gausslobatto(np)
    if(trim(IO_SEopt%topology) == "cube") then
      print *,' Initializing cube element coordinates'
      if(MeshUseMeshFile) then
        call MeshSetCoordinates(elem)
      else
        do ie=1,IO_SEopt%nelem
          call set_corner_coordinates(elem(ie))
        end do
        call assign_node_numbers_to_elem(elem,GridVertex)
      endif
      do ie=1,IO_SEopt%nelem
        call cube_init_atomic(elem(ie),gp%points)
      end do
    endif

    ! Initialize mass matrix
    !-----------------------
    call init_mass_matrix_mod(IO_SEopt)
    call init_edge_mod(IO_SEopt)
    print *,' Running mass_matrix'
    call mass_matrix(elem)
    print *,' Running mass_matrix done'
    if(trim(IO_SEopt%topology) == "cube") then
      area = 0
      do ie=1,IO_SEopt%nelem
        area = area + sum(elem(ie)%mp(:,:)*elem(ie)%metdet(:,:))
      end do
      area = 4*DD_PI/area
      print *,' Element Total Area = ',area
      do ie=1,IO_SEopt%nelem
        call cube_init_atomic(elem(ie),gp%points,area)
        call rotation_init_atomic(elem(ie),IO_SEopt%rot_type)
      end do
    endif
    print *,' Re-running mass_matrix'
    call mass_matrix(elem)

    ! Determine global degree of freedom for each gridpoint
    !-------------------------------------------------------
    call init_dof_mod(IO_SEopt)
    call global_dof(elem)

    ! Create Unique Indices
    !------------------------
    do ie=1,IO_SEopt%nelem
      call CreateUniqueIndex(elem(ie)%GlobalId,elem(ie)%gdofP,elem(ie)%idxP)
    end do

    call SetElemOffset(elem,IO_SEopt%GlobalUniqueCols)
    do ie=1,IO_SEopt%nelem
      elem(ie)%idxV => elem(ie)%idxP
    end do

    ! Clean up: Done with GridEdge and GridVertex
    !---------------------------------------------
    deallocate(GridEdge)
    do ie =1,IO_SEopt%nelem
      call deallocate_gridvertex_nbrs(GridVertex(ie))
    end do
    deallocate(GridVertex)

    !======================
    ! Test Global Integral
    !======================
    allocate(Htest(np,np,IO_SEopt%nelem))
    Htest(:,:,:)= 1.0_real_kind
    I_sphere    = 0.0_real_kind
    do ie=1,IO_SEopt%nelem
    do jj=1,np
    do ii=1,np
      Darea = elem(ie)%mp(ii,jj)*elem(ie)%metdet(ii,jj)
      I_sphere = I_sphere + Darea*Htest(ii,jj,ie)
    end do
    end do
    end do
    I_sphere = I_sphere /(4.0_real_kind*DD_PI)
    deallocate(Htest)

    ! Calc some min/max/avg values
    !------------------------------
    min_area   = 1.0d99
    max_area   = 0.0
    max_ratio  = 0.0
    min_max_eig= 1.0d99
    max_max_eig= 0.0
    min_min_eig= 1.0d99
    max_min_eig= 0.0
    min_max_dx = 1.0d99
    max_max_dx = 0.0
    min_min_dx = 1.0d99
    max_min_dx = 0.0
    avg_area   = 0.0_real_kind
    avg_max_eig= 0.0_real_kind
    avg_min_eig= 0.0_real_kind
    avg_max_dx = 0.0_real_kind
    avg_min_dx = 0.0_real_kind
    do ie=1,IO_SEopt%nelem
      elem(ie)%area = sum(elem(ie)%spheremp(:,:))
      min_area      = min(min_area,elem(ie)%area)
      max_area      = max(max_area,elem(ie)%area)
      min_max_eig   = min(min_max_eig,elem(ie)%max_eig)
      max_max_eig   = max(max_max_eig,elem(ie)%max_eig)
      max_ratio     = max(max_ratio,elem(ie)%max_eig_ratio)
      min_min_eig   = min(min_min_eig,elem(ie)%min_eig)
      max_min_eig   = max(max_min_eig,elem(ie)%min_eig)
      min_max_dx    = min(min_max_dx,elem(ie)%dx_long)
      max_max_dx    = max(max_max_dx,elem(ie)%dx_long)
      min_min_dx    = min(min_min_dx,elem(ie)%dx_short)
      max_min_dx    = max(max_min_dx,elem(ie)%dx_short)
      avg_area   = avg_area   + elem(ie)%area
      avg_max_eig= avg_max_eig+ elem(ie)%max_eig
      avg_min_eig= avg_min_eig+ elem(ie)%min_eig
      avg_max_dx = avg_max_dx + elem(ie)%dx_long
      avg_min_dx = avg_min_dx + elem(ie)%dx_short
    end do
    avg_area   = avg_area   /dble(IO_SEopt%nelem)
    avg_max_eig= avg_max_eig/dble(IO_SEopt%nelem)
    avg_min_eig= avg_min_eig/dble(IO_SEopt%nelem)
    avg_max_dx = avg_max_dx /dble(IO_SEopt%nelem)
    avg_min_dx = avg_min_dx /dble(IO_SEopt%nelem)

    ! Physical units for area
    !-------------------------
    min_area = min_area*rearth*rearth/1000000_real_kind
    max_area = max_area*rearth*rearth/1000000_real_kind
    avg_area = avg_area*rearth*rearth/1000000_real_kind

    ! for an equation du/dt = i c u, leapfrog is stable for |c u dt| < 1
    ! Consider a gravity wave at the equator, c=340m/s  
    ! u = exp(i kmax x/ a ) with x = longitude,  and kmax =  pi a / dx, 
    ! u = exp(i pi x / dx ),   so du/dt = c du/dx becomes du/dt = i c pi/dx u
    ! stable for dt < dx/(c*pi)
    ! CAM 26 level AMIP simulation: max gravity wave speed 341.75 m/s
    !--------------------------------------------------------------------------
    write(iulog,*) ' '
    write(iulog,*) '  Running Global Integral Diagnostic...'
    write(iulog,*) '  Area of unit sphere is',I_sphere
    write(iulog,*) '  Should be 1.0 to round off...'
    write(iulog,'(a,f9.3)') '  Element area:  max/min',(max_area/min_area)
    if(.not.MeshUseMeshFile) then
      write(iulog,'(a,f6.3,f8.2)') '  Average equatorial node spacing (deg, km) = ', &
                                      dble(90)/dble(IO_SEopt%ne*(np-1)),             &
                                      DD_PI*rearth/(2000.0d0*dble(IO_SEopt%ne*(np-1)))
    endif
    write(iulog,'(a,2f9.3)' ) '  Min eigenvalue of Dinv (min, max): ', &
                                        min_min_eig, max_min_eig
    write(iulog,'(a,2f9.3)' ) '  Max eigenvalue of Dinv (min, max): ', &
                                        min_max_eig, max_max_eig
    write(iulog,'(a,1e11.3)') '  Max eigenvalue ratio (element distortion): ', &
                                        max_ratio
    write(iulog,'(a,3f8.2)' ) '  dx for CFL (based on Dinv eigenvalue):ave,min,max= ', &
                                        avg_min_dx, min_min_dx, max_min_dx
    write(iulog,'(a,3f8.2)' ) '  dx based on sqrt element area: ave,min,max = ', &
                  sqrt(avg_area)/(np-1),sqrt(min_area)/(np-1),sqrt(max_area)/(np-1)

    ! Initialize Derivative Structure
    !---------------------------------
    xtmp=nc
    do ii=1,(nc+1)
      fvm_corners(ii)= 2*(ii-1)/xtmp - 1  ! [-1,1] including end points
    end do
    do ii=1,nc
      fvm_points(ii)= (fvm_corners(ii)+fvm_corners(ii+1))/2
    end do
    xtmp=(nep-1)
    do ii=1,nep
      spelt_refnep(ii)= 2*(ii-1)/xtmp - 1
    end do
    call init_derivative_mod(IO_SEopt)
    call derivinit(deriv,fvm_corners,fvm_points,spelt_refnep)

    !======================
    ! Print CFL
    !======================
    ! Eigenvalues calculated by folks at UMich (Paul U & Jared W)
    !----------------------------------------------------------------
    select case (np)
      case (2)
            lambda_max = 0.5d0
      case (3)
            lambda_max = 1.5d0
            lambda_vis = 12.0d0
      case (4)
            lambda_max = 2.74d0
            lambda_vis = 30.0d0
      case (5)
            lambda_max = 4.18d0
            lambda_vis = 91.6742d0
      case (6)
            lambda_max = 5.86d0
            lambda_vis = 190.1176d0
      case (7)
            lambda_max = 7.79d0
            lambda_vis = 374.7788d0
      case (8)
            lambda_max = 10.0d0
            lambda_vis = 652.3015d0
      case DEFAULT
            lambda_max = 0.0d0
            lambda_vis = 0.0d0
    end select
    if(lambda_max.eq.0d0) then
      write(iulog,*) 'lambda_max not calculated for NP = ',np
      write(iulog,*) 'Estimate of gravity wave timestep will be incorrect'
    endif
    min_gw = minval(gp%weights)

    do ie=1,IO_SEopt%nelem
      elem(ie)%variable_hyperviscosity = 1.0
    end do

    !*******************************************************************
    ! dtnu NEEDS to be set based on timestep, splits, tstep_type, etc..
    ! for now just set it to 1 so the program wont fail. For meaningful
    ! output and variable_hyperviscosity values, this value shoud be set
    ! correctly.
    !*******************************************************************
    dtnu = 1.d0

    if(IO_SEopt%hypervis_power /= 0) then
      min_hypervis = 1d99
      max_hypervis = 0
      avg_hypervis = 0

      ! viscosity in namelist specified for smallest element:
      !---------------------------------------------------------------
      max_unif_dx = min_max_dx
      if(IO_SEopt%fine_ne > 0) then
        ! viscosity in namelist specified for regions with a resolution
        ! equivilant to a uniform grid with ne=fine_ne !!! og: is this in km?
        ! np=4? yes
        !-------------------------------------------------------------------------------
        max_unif_dx = (111.28*30)/dble(IO_SEopt%fine_ne)
      endif

      ! note: if L = eigenvalue of metinv, then associated length scale (km) is
      ! dx = 1.0d0/( sqrt(L)*0.5d0*dble(np-1)*rrearth*1000.0d0)
      !
      !       for viscosity *tensor*, we take at each point: 
      !            nu1 = nu*(dx1/max_unif_dx)**3.2      dx1 associated with
      !            eigenvalue 1
      !            nu2 = nu*(dx2/max_unif_dx)**3.2      dx2 associated with
      !            eigenvalue 2
      !       with this approach:
      !          - with this formula, no need to adjust for CFL violations
      !          - if nu comes from a 3.2 scaling that is stable for coarse and
      !          fine resolutions,
      !            this formulat will be stable.  
      !          - gives the correct answer in long skinny rectangles:
      !            large viscosity in the long direction, small viscosity in the
      !            short direction 
      !-------------------------------------------------------------------------------- 
      max_eig_hypervis = 0
      do ie=1,IO_SEopt%nelem
        ! variable viscosity based on map from ulatlon -> ucontra
        !---------------------------------------------------------
        ! dx_long
        !---------
        elem(ie)%variable_hyperviscosity = &
                 sqrt((elem(ie)%dx_long/max_unif_dx)**IO_SEopt%hypervis_power)
        ! dx_short
        !---------
!       elem(ie)%variable_hyperviscosity = &
!                sqrt((elem(ie)%dx_short/1.25/min_min_dx)**IO_SEopt%hypervis_power)
        ! geometric mean:
        !----------------
!       elem(ie)%variable_hyperviscosity = &
!                sqrt((sqrt(elem(ie)%dx_short*elem(ie)%dx_long)/min_min_dx/1.9)**IO_SEopt%hypervis_power)
        ! average
        !---------
!       elem(ie)%variable_hyperviscosity = &
!                sqrt(((elem(ie)%dx_short+elem(ie)%dx_long)/3.8/min_min_dx)**IO_SEopt%hypervis_power)

        ! sqrt(area), area = unit sphere.  other variables in km
        !---------------------------------------------------------
!        elem(ie)%variable_hyperviscosity = &
!                sqrt((sqrt(elem(ie)%area)*Rearth/2.4e3/max_unif_dx)**IO_SEopt%hypervis_power)

        elem(ie)%hv_courant = dtnu*(elem(ie)%variable_hyperviscosity(1,1)**2) &
                                  *(lambda_vis**2)*((rrearth*elem(ie)%max_eig)**4)
        ! Check to see if this is stable
        !---------------------------------
        if(elem(ie)%hv_courant.gt.IO_SEopt%max_hypervis_courant) then
          stable_hv = sqrt( IO_SEopt%max_hypervis_courant &
                           /(dtnu*(lambda_vis**2)*((rrearth*elem(ie)%max_eig)**4)))
          ! make sure that: elem(ie)%hv_courant <=  max_hypervis_courant 
          !-------------------------------------------------------------
          elem(ie)%variable_hyperviscosity = stable_hv
          elem(ie)%hv_courant = dtnu*(stable_hv**2) &
                                    *(lambda_vis**2)*((rrearth*elem(ie)%max_eig)**4)
        endif
        max_eig_hypervis = max(max_eig_hypervis,elem(ie)%hv_courant/dtnu)
        min_hypervis     = min(min_hypervis,elem(ie)%variable_hyperviscosity(1,1))
        max_hypervis     = max(max_hypervis,elem(ie)%variable_hyperviscosity(1,1))
        avg_hypervis     = avg_hypervis + elem(ie)%variable_hyperviscosity(1,1)
      end do
      avg_hypervis = avg_hypervis/dble(IO_SEopt%nelem)


      ! apply DSS (aka assembly procedure) to variable_hyperviscosity (makes
      ! continuous)
      !------------------------------------------------------------------------------------
      call initEdgeBuffer(edgebuf,1)
      allocate(zeta(np,np,IO_SEopt%nelem))
      do ie=1,IO_SEopt%nelem
        zeta(:,:,ie) = elem(ie)%variable_hyperviscosity(:,:)*elem(ie)%spheremp(:,:)
        call edgeVpack(edgebuf,zeta(1,1,ie),1,0,elem(ie)%desc)
      end do
      do ie=1,IO_SEopt%nelem
        call edgeVunpack(edgebuf,zeta(1,1,ie),1,0,elem(ie)%desc)
        elem(ie)%variable_hyperviscosity(:,:) = zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
      end do
      call FreeEdgeBuffer(edgebuf)

      ! replace hypervis w/ bilinear based on continuous corner values
      !----------------------------------------------------------------
      do ie=1,IO_SEopt%nelem
        noreast = elem(ie)%variable_hyperviscosity(np,np)
        nw      = elem(ie)%variable_hyperviscosity( 1,np)
        se      = elem(ie)%variable_hyperviscosity(np, 1)
        sw      = elem(ie)%variable_hyperviscosity( 1, 1)
        do ii=1,np
          xx = gp%points(ii)
          do jj=1,np
            yy = gp%points(jj)
            elem(ie)%variable_hyperviscosity(ii,jj)=( (1.0d0-xx)*(1.0d0-yy)*sw     &
                                                     +(1.0d0-xx)*(yy+1.0d0)*nw     &
                                                     +(xx+1.0d0)*(1.0d0-yy)*se     &
                                                     +(xx+1.0d0)*(yy+1.0d0)*noreast)/4.d0
          end do
        end do
      end do
      deallocate(zeta)
    elseif(IO_SEopt%hypervis_scaling /= 0) then
      ! tensorHV.  New eigenvalues are the eigenvalues of the tensor V
      ! formulas here must match what is in cube_mod.F90
      ! for tensorHV, we scale out the rearth dependency
      !------------------------------------------------------------------
      lambda = max_max_eig**2
      max_eig_hypervis = (lambda_vis**2)*(max_max_eig**4)*(lambda**(-IO_SEopt%hypervis_scaling/2))
    else
      ! constant coefficient formula:
      !-------------------------------
      max_eig_hypervis = (lambda_vis**2)*((rrearth*max_max_eig)**4)
    endif

    if(IO_SEopt%hypervis_scaling /= 0) then
      !this is a code for smoothing V for tensor HV
      !------------------------------------------------

      ! if nu=0, we have extra 4 DSS here!
      ! rowind, colind are from 1 to 2 cause they correspond to 2D tensor in
      ! lat/lon
      !-----------------------------------------------------------------------------

      ! IF DSSED V NEEDED
      !-------------------
      call initEdgeBuffer(edgebuf,1)
      allocate(zeta(np,np,IO_SEopt%nelem))
      do rowind=1,2
      do colind=1,2
        do ie=1,IO_SEopt%nelem
          zeta(:,:,ie) = elem(ie)%tensorVisc(rowind,colind,:,:)*elem(ie)%spheremp(:,:)
          call edgeVpack(edgebuf,zeta(1,1,ie),1,0,elem(ie)%desc)
        end do
        do ie=1,IO_SEopt%nelem
          call edgeVunpack(edgebuf,zeta(1,1,ie),1,0,elem(ie)%desc)
          elem(ie)%tensorVisc(rowind,colind,:,:) = zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
        end do
      end do 
      end do 
      deallocate(zeta)
      call FreeEdgeBuffer(edgebuf)

      ! IF BILINEAR MAP OF V NEEDED
      !------------------------------
      do rowind=1,2
      do colind=1,2
        ! replace hypervis w/ bilinear based on continuous corner values
        !----------------------------------------------------------------
        do ie=1,IO_SEopt%nelem
          noreast = elem(ie)%tensorVisc(rowind,colind,np,np)
          nw      = elem(ie)%tensorVisc(rowind,colind, 1,np)
          se      = elem(ie)%tensorVisc(rowind,colind,np, 1)
          sw      = elem(ie)%tensorVisc(rowind,colind, 1, 1)
          do ii=1,np
            xx = gp%points(ii)
            do jj=1,np
              yy = gp%points(jj)
              elem(ie)%tensorVisc(rowind,colind,ii,jj)=( (1.0d0-xx)*(1.0d0-yy)*sw     &
                                                        +(1.0d0-xx)*(yy+1.0d0)*nw     &
                                                        +(xx+1.0d0)*(1.0d0-yy)*se     &
                                                        +(xx+1.0d0)*(yy+1.0d0)*noreast)/4.d0
            end do
          end do
        end do
      end do !rowind
      end do !colind
    endif
    deallocate(gp%points)
    deallocate(gp%weights)

#if 0
!===================================================================================
! This output will be meanigful if nu, and a bunch of other values are set, 
! There is no effect on SE differential options (other than variable_hyperviscosity)
! DROP THIS SECTION FOR NOW
!===================================================================================
    write(iulog,'(a,f10.2)') 'CFL estimates in terms of S=time step stability region'
    write(iulog,'(a,f10.2)') '(i.e. advection w/leapfrog: S=1, viscosity w/forward Euler: S=2)'
    if(rk_stage_user>0) then
       write(iulog,'(a,f10.2,a)') 'SSP preservation (120m/s) RKSSP euler step dt  < S *', &
                                   min_gw/(120.0d0*max_max_eig*rrearth),'s'
    endif
    write(iulog,'(a,f10.2,a)') 'Stability: advective (120m/s)   dt_tracer < S *', &
                                1/(120.0d0*max_max_eig*lambda_max*rrearth),'s'
    write(iulog,'(a,f10.2,a)') 'Stability: gravity wave(342m/s)   dt_dyn  < S *', &
                                1/(342.0d0*max_max_eig*lambda_max*rrearth),'s'
    if(nu>0) then
      if(IO_SEopt%hypervis_order==1) then
        write(iulog,'(a,f10.2,a)') 'Stability: viscosity dt < S *',1/(((rrearth*max_max_eig)**2)*lambda_vis),'s'
      endif
      if(IO_SEopt%hypervis_order==2) then
        ! counrant number = dtnu*max_eig_hypervis  < S
        !  dt < S  1/nu*max_eig
        !------------------------------------------------------------------
        write(iulog,'(a,f10.2,a)') "Stability: nu_q   hyperviscosity dt < S *", 1/(nu_q*max_eig_hypervis),'s'
        write(iulog,'(a,f10.2,a)') "Stability: nu_vor hyperviscosity dt < S *", 1/(nu*max_eig_hypervis),'s'

        ! bug in nu_div implimentation:
        ! we apply nu_ration=(nu_div/nu) in laplace, so it is applied 2x
        ! making the effective nu_div = nu * (nu_div/nu)**2 
        ! should be fixed - but need all CAM defaults adjusted, 
        ! so we have to coordiante with CAM
        !------------------------------------------------------------------
        nu_div_actual = nu_div**2/nu
        write(iulog,'(a,f10.2,a)') "Stability: nu_div hyperviscosity dt < S *", 1/(nu_div_actual*max_eig_hypervis),'s'
      endif
    endif
    if(nu_top>0) then
      write(iulog,'(a,f10.2,a)') 'TOP3 viscosity CFL: dt < S*', &
                                1.0d0/(4*nu_top*((rrearth*max_max_eig)**2)*lambda_vis),'s'
    endif
    if(IO_SEopt%hypervis_power /= 0) then
      write(iulog,'(a,3e11.4)')'Hyperviscosity (dynamics): ave,min,max = ', &
                                nu*(/avg_hypervis**2,min_hypervis**2,max_hypervis**2/)
!     print*, 'fine_ne = ', IO_SEopt%fine_ne
!     print*, 'Using max_unif_dx = ', max_unif_dx
    endif
#endif

    ! Write out elem diagnostics
    !----------------------------
    if(.FALSE.) then
      do ie=1,IO_SEopt%nelemd
        my_unit   = 23
        my_number = elem(ie)%GlobalId
        FileName  = 'Elem_NNNNN_Diag.txt'
        write(FileName(6:10),'(i5.5)') my_number
        open(my_unit,file='./Element_Diag/'//trim(FileName),form='formatted')
          write(my_unit,*) ' Element Number=',my_number
          write(my_unit,*) '      LocalId=',elem(ie)%LocalId
          write(my_unit,*) '     GlobalId=',elem(ie)%GlobalId
          write(my_unit,*) '      FaceNum=',elem(ie)%FaceNum
          write(my_unit,*) '      corners=',elem(ie)%corners
          write(my_unit,*) '    corners3D=',elem(ie)%corners3D
          write(my_unit,*) '       u2qmap=',elem(ie)%u2qmap
          write(my_unit,*) '      spherep=',elem(ie)%spherep
          write(my_unit,*) '        cartp=',elem(ie)%cartp
          write(my_unit,*) ' D=',elem(ie)%D
          write(my_unit,*) ' Dinv=',elem(ie)%Dinv
          write(my_unit,*) ' '
          write(my_unit,*) 'elem()%idxP: '
          write(my_unit,*) '            idxP%is=',elem(ie)%idxP%is
          write(my_unit,*) '            idxP%ie=',elem(ie)%idxP%ie
          write(my_unit,*) '  idxP%NumUniquePts=',elem(ie)%idxP%NumUniquePts
          write(my_unit,*) 'idxP%UniquePtOffset=',elem(ie)%idxP%UniquePtOffset
          write(my_unit,*) '            idxP%ia=',elem(ie)%idxP%ia
          write(my_unit,*) '            idxP%ja=',elem(ie)%idxP%ja
          write(my_unit,*) ' '
          write(my_unit,*) 'elem()%vertex: '
          write(my_unit,*) 'vertex%face_number=',elem(ie)%vertex%face_number
          write(my_unit,*) '     vertex%number=',elem(ie)%vertex%number
          write(my_unit,*) ' vertex%SpaceCurve=',elem(ie)%vertex%SpaceCurve
          write(my_unit,*) '       vertex%nbrs=',elem(ie)%vertex%nbrs
          write(my_unit,*) '  vertex%nbrs_face=',elem(ie)%vertex%nbrs_face
        close(my_unit)
      end do
    endif

    ! End Routine
    !--------------
    return
  end subroutine Init_BASEelem
  !==================================================================


  !==================================================================
  subroutine SE_Get_Ncol(O_ncol)
    !
    ! SE_Get_Ncol: Return the number of unique SE gridpoints.
    !
    !==========================================================
    !
    ! Passed Variables
    !----------------
    integer,intent(out):: O_ncol

    if(SEinitialized) then
      O_ncol = SEncol
    else
      write(iulog,*) "WARNING: SEinitialized = ",SEinitialized
      O_ncol = 0
    endif

    ! End Routine
    !--------------
    return
  end subroutine SE_Get_Ncol
  !==================================================================


  !==================================================================
  subroutine SE_Get_Nelem(O_nelem)
    !
    ! SE_Get_Nelem: Return the number of SE elements.
    !
    !==========================================================
    !
    ! Passed Variables
    !----------------
    integer,intent(out):: O_nelem

    if(SEinitialized) then
      O_nelem = SEnelem
    else
      write(iulog,*) "WARNING: SEinitialized = ",SEinitialized
      O_nelem = 0
    endif

    ! End Routine
    !--------------
    return
  end subroutine SE_Get_Nelem
  !==================================================================


  !==================================================================
  subroutine SE_Get_Face(I_SEelem,O_FaceNumber)
    type(SEelem_t)      ,intent(in ):: I_SEelem(SEnelem)
    real(kind=real_kind),intent(out):: O_FaceNumber(np,np,SEnelem)

    integer ie
 
    do ie=1,SEnelem
     O_FaceNumber(:,:,ie) = real(I_SEelem(ie)%face_number,real_kind)
    end do

    ! End Routine
    !--------------
    return
  end subroutine SE_Get_Face
  !==================================================================


  !==================================================================
  subroutine SE_Get_Number(I_SEelem,O_Number)
    type(SEelem_t)      ,intent(in ):: I_SEelem(SEnelem)
    real(kind=real_kind),intent(out):: O_Number(np,np,SEnelem)

    integer ie
 
    do ie=1,SEnelem
     O_Number(:,:,ie) = real(I_SEelem(ie)%number,real_kind)
    end do

    ! End Routine
    !--------------
    return
  end subroutine SE_Get_Number
  !==================================================================

end module SE_Element_mod

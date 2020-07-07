module SE_Options
!============================================================================
!
! Purpose: Provide a SEoptions data structure containing parameters and
!          values which control SE grid and processing behavior. 
!
! Author: Patrick Callaghan
!
! Description: This module initializes the SEoptions datastructure with 
!              default values. The values are over written with values
!              contained in the provided namelist file. The module also
!              provides as global values constants which make the usage
!              option values readable when they are utilized.
!
! VERSION: 1.0
!
! Revisions:
!   Feb 2015 - Original Version
!
! Usage:
!   Init_SEoptions(O_SEopt,I_NameListFile)
!                 Routine to initialize the given SEoptions_t datastructure.
!                 If the optional namelist file is provided, it is read in 
!                 to set option values.
!
!   set_VarMesh_dimensions() : PRIVATE routine to adjust values that change 
!                 when a variable mesh is specified.
!
!============================================================================
  ! Useful modules
  !------------------
  use SE_Constants

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private
  save

  ! Routines
  !----------
  public :: set_VarMesh_dimensions
  public :: Init_SEoptions
  private:: find_group_name

  ! Spectral Element Option Constants/Types
  !----------------------------------------
  public :: west
  public :: east
  public :: south
  public :: north
  public :: swest
  public :: seast
  public :: nwest
  public :: neast
  public :: INTERNAL_EDGE
  public :: EXTERNAL_EDGE
  public :: RECURSIVE
  public :: KWAY
  public :: VOLUME
  public :: WRECURSIVE
  public :: SFCURVE

  private:: max_elements_attached_to_node
  private:: s_nv                         
  private:: max_corner_elem              
  private:: max_neigh_edges              
  private:: ne
  private:: MeshFile 
  private:: nelem                       
  private:: nelemd                      
  private:: npart
  private:: GlobalUniqueCols
  private:: hypervis_power       
  private:: hypervis_scaling    
  private:: fine_ne              
  private:: max_hypervis_courant 
  private:: hypervis_order       
  private:: partmethod
  private:: cubed_sphere_map  
  private:: topology 
  private:: rot_type 

  public :: SEoptions_t

  ! Control - indexing of element sides
  !-------------------------------------
  integer,parameter:: west  = 1
  integer,parameter:: east  = 2
  integer,parameter:: south = 3
  integer,parameter:: north = 4
  integer,parameter:: swest = 5
  integer,parameter:: seast = 6
  integer,parameter:: nwest = 7
  integer,parameter:: neast = 8

  integer,parameter:: INTERNAL_EDGE  = 0
  integer,parameter:: EXTERNAL_EDGE  = 1

  ! Type of partitioning methods 
  !-------------------------------
  integer,parameter:: RECURSIVE  = 0
  integer,parameter:: KWAY       = 1
  integer,parameter:: VOLUME     = 2
  integer,parameter:: WRECURSIVE = 3
  integer,parameter:: SFCURVE    = 4

  ! SET DEFAULT VALUES: Dimensions and constants for SPELT
  !----------------------------------------------------------
  integer:: max_elements_attached_to_node  = 4
  integer:: s_nv                           = 6
  integer:: max_corner_elem                = 1  !max_elements_attached_to_node-3
  integer:: max_neigh_edges                = 8  ! 4 + 4*max_corner_elem
  integer:: ne                             = 0  ! 30 ! 60 ! 120 !PFC 30
  character(len=MAX_STRING_LEN):: MeshFile = 'NULL' 
  integer:: nelem                               ! total number of elements
  integer:: nelemd                              ! number of elements per MPI task
  integer:: npart                          = 1
  integer:: GlobalUniqueCols

  ! SET DEFAULT VALUES:
  ! three types of hyper viscosity are supported right now:
  !  (1)  const hv:   nu * del^2 del^2
  !  (2) scalar hv:   nu(lat,lon) * del^2 del^2
  !  (3) tensor hv,   nu * ( \div * tensor * \grad ) * del^2
  !
  !  (1) default:  hypervis_power=0, hypervis_scaling=0
  !  (2) Original version for var-res grids. (M. Levy)
  !               scalar coefficient within each element
  !               hypervisc_scaling=0
  !               set hypervis_power>0 and set fine_ne, max_hypervis_courant
  !  (3) tensor HV var-res grids 
  !               tensor within each element:
  !               set hypervis_scaling > 0 (typical values would be 3.2 or 4.0)
  !               hypervis_power=0
  !               (\div * tensor * \grad) operator uses cartesian laplace
  !------------------------------------------------------------------------
  real(kind=real_kind):: hypervis_power       = 0    ! if not 0, use variable hyperviscosity 
                                                            !   based on element area
  real(kind=real_kind):: hypervis_scaling     = 3.2    ! 0    ! use tensor hyperviscosity
  real(kind=real_kind):: fine_ne              = 120    ! 0
  real(kind=real_kind):: max_hypervis_courant = 9991.9 ! 1d99 ! upper bound for Courant number
  real(kind=real_kind):: hypervis_order       = 2  

  ! SET DEFAULT VALUES:
  !  partition method
  !-----------------------------------
  integer             :: partmethod           = SFCURVE

  ! SET DEFAULT VALUES:
  ! cubed_sphere_map:
  !       -1 = chosen at run time
  !        0 = equi-angle Gnomonic (default)
  !        1 = equi-spaced Gnomonic (not yet coded)
  !        2 = element-local projection  (for var-res)
  !        3 = parametric (not yet coded)
  !---------------------------------------------
  integer            :: cubed_sphere_map  = -1    

  character(len=MAX_STRING_LEN):: topology = 'cube'          ! options: "cube" is supported
  character(len=MAX_STRING_LEN):: rot_type = 'contravariant' 

  ! SE_options_t
  !----------------
  type SEoptions_t
    integer                      :: max_elements_attached_to_node
    integer                      :: s_nv                         
    integer                      :: max_corner_elem
    integer                      :: max_neigh_edges
    integer                      :: ne
    integer                      :: nelem              
    integer                      :: nelemd             
    integer                      :: npart
    integer                      :: GlobalUniqueCols
    real(kind=real_kind)         :: hypervis_power  
    real(kind=real_kind)         :: hypervis_scaling
    real(kind=real_kind)         :: fine_ne         
    real(kind=real_kind)         :: max_hypervis_courant
    real(kind=real_kind)         :: hypervis_order      
    integer                      :: partmethod 
    integer                      :: cubed_sphere_map 
    character(len=MAX_STRING_LEN):: topology
    character(len=MAX_STRING_LEN):: rot_type
    character(len=MAX_STRING_LEN):: MeshFile
    character(len=MAX_FILE_LEN)  :: InputPath
    character(len=MAX_FILE_LEN)  :: OutputPath
    character(len=MAX_STRING_LEN):: SCRIP_filename
    character(len=MAX_STRING_LEN):: LATLON_filename
    character(len=MAX_STRING_LEN):: PHYS_filename
    character(len=MAX_STRING_LEN):: GRID_filename
    logical                      :: create_SCRIP_file
    logical                      :: create_LATLON_file
    logical                      :: create_PHYS_file
    logical                      :: create_GRID_file
  end type SEoptions_t

contains

  !==========================================================================
  subroutine set_VarMesh_dimensions()
    ! set_VarMesh_dimensions:
    !
    !  To accomodate Variable Resolution meshes, 
    !  reset dimension values when mesh_file is used.
    !================================================================

    ! new "params"
    !--------------
    max_elements_attached_to_node = 7  ! variable resolution
    s_nv                          = 2*max_elements_attached_to_node

    !recalculate these
    !------------------
    max_corner_elem = max_elements_attached_to_node-3
    max_neigh_edges = 4 + 4*max_corner_elem

    ! End Routine
    !-------------
    return
  end subroutine set_VarMesh_dimensions
  !==========================================================================


  !==========================================================================
  subroutine Init_SEoptions(O_SEopt,I_NameListFile)
    ! Init_SEoptions:
    !
    !  If a file is given, open it and read SE options, 
    !  otherwise read in values from 'SEoptions.nl'
    !================================================================
    use err_exit,only: iulog,endrun
    ! Passed Variables
    !-------------------
    type(SEoptions_t)         ,intent(out):: O_SEopt
    character(len=*) ,optional,intent(in ):: I_NameListFile
    ! Local Values
    !---------------
    character(len=MAX_STRING_LEN):: GridPath           = './' 
    character(len=MAX_STRING_LEN):: GridName           = 'NULL' 
    character(len=MAX_FILE_LEN)  :: InputPath          = './'
    character(len=MAX_FILE_LEN)  :: OutputPath         = './'
    character(len=MAX_STRING_LEN):: SCRIP_filename     = 'NULL' 
    character(len=MAX_STRING_LEN):: LATLON_filename    = 'NULL' 
    character(len=MAX_STRING_LEN):: PHYS_filename      = 'NULL' 
    character(len=MAX_STRING_LEN):: GRID_filename      = 'NULL' 
    logical                      :: create_SCRIP_file  = .true.
    logical                      :: create_LATLON_file = .true.
    logical                      :: create_PHYS_file   = .false.
    logical                      :: create_GRID_file   = .false.

    real(kind=real_kind):: rotate_grid
    real(kind=real_kind):: nu, nu_s, nu_p, nu_q, nu_div, nu_top
    integer             :: rk_stage_user
    integer             :: NLunit,ierr

    namelist /ctl_nl/ partmethod,           &  ! Mesh partitioning method (METIS)
                      npart,                &
                      topology,             &  ! Mesh topology
                      cubed_sphere_map,     &
                      ne,                   &  ! element resolution factor
                      MeshFile,             &  ! Name of mesh file
                      fine_ne,              &
                      InputPath,            &
                      OutputPath,           &
                      create_SCRIP_file,    &
                      SCRIP_filename,       &
                      create_LATLON_file,   &
                      LATLON_filename,      &
                      create_PHYS_file,     &
                      PHYS_filename,        &
                      create_GRID_file,     &
                      GRID_filename,        &
                      hypervis_scaling,     &
                      hypervis_power,       &
                      hypervis_order,       &
                      max_hypervis_courant, &
                      rot_type,             &
                      rotate_grid,          &
                      nu,                   &
                      nu_s,                 &
                      nu_p,                 &
                      nu_q,                 &
                      nu_div,               &
                      nu_top,               &
                      rk_stage_user

    namelist /input_nl/ GridPath, &
                        GridName
    
    ! Optionally read in values from the given namelist file
    !-------------------------------------------------------
    NLunit = 21
    if(present(I_NameListFile)) then
      ! Open namelist file and read values
      !-------------------------------------
      print *,' OPENING: namelist file = ',trim(I_NameListFile)
      open(NLunit,file=trim(I_NameListFile),status='old')
      call find_group_name(NLunit,'input_nl',status=ierr)
      if(ierr.eq.0) then
        read(NLunit,input_nl,iostat=ierr)
        if(ierr.ne.0) then
          call endrun("ERROR Reading input_nl namelist values")
        endif
      else
        call endrun("ERROR Finding input_nl namelist in file")
      endif
      close(NLunit)
    endif

    ! Construct file names
    !---------------------------
    InputPath       = trim(GridPath)
    OutputPath      = trim(GridPath)
    MeshFile        = trim(GridName)//'_EXODUS.nc'
    SCRIP_filename  = trim(GridName)//'_np4_SCRIP.nc'
    LATLON_filename = trim(GridName)//'_np4_LATLON.nc'

    ! Copy final values to the output data structure
    !------------------------------------------------
    O_SEopt%max_elements_attached_to_node = max_elements_attached_to_node
    O_SEopt%s_nv                          = s_nv                         
    O_SEopt%max_corner_elem               = max_corner_elem
    O_SEopt%max_neigh_edges               = max_neigh_edges
    O_SEopt%ne                            = ne
    O_SEopt%nelem                         = nelem              
    O_SEopt%nelemd                        = nelemd             
    O_SEopt%npart                         = npart
    O_SEopt%GlobalUniqueCols              = GlobalUniqueCols
    O_SEopt%hypervis_power                = hypervis_power  
    O_SEopt%hypervis_scaling              = hypervis_scaling
    O_SEopt%fine_ne                       = fine_ne         
    O_SEopt%max_hypervis_courant          = max_hypervis_courant
    O_SEopt%hypervis_order                = hypervis_order      
    O_SEopt%partmethod                    = partmethod 
    O_SEopt%cubed_sphere_map              = cubed_sphere_map 
    O_SEopt%topology                      = topology
    O_SEopt%rot_type                      = rot_type
    O_SEopt%MeshFile                      = MeshFile
    O_SEopt%InputPath                     = InputPath
    O_SEopt%OutputPath                    = OutputPath
    O_SEopt%SCRIP_filename                = SCRIP_filename
    O_SEopt%LATLON_filename               = LATLON_filename
    O_SEopt%PHYS_filename                 = PHYS_filename
    O_SEopt%GRID_filename                 = GRID_filename
    O_SEopt%create_SCRIP_file             = create_SCRIP_file
    O_SEopt%create_LATLON_file            = create_LATLON_file
    O_SEopt%create_PHYS_file              = create_PHYS_file
    O_SEopt%create_GRID_file              = create_GRID_file

    ! End Routine
    !-------------
    return
  end subroutine Init_SEoptions
  !==========================================================================


  !==========================================================================
  subroutine find_group_name(unit, group, status)
    ! Purpose:
    ! Search a file that contains namelist input for the specified namelist group name.
    ! Leave the file positioned so that the current record is the first record of the
    ! input for the specified group.
    !
    ! Method:
    ! Read the file line by line.  Each line is searched for an '&' which may only
    ! be preceded by blanks, immediately followed by the group name which is case
    ! insensitive.  If found then backspace the file so the current record is the
    ! one containing the group name and return success.  Otherwise return -1.
    !
    ! Author:  B. Eaton, August 2007
    !==================================================================
    !
    ! Passed Variables
    !-------------------
    integer,         intent(in) :: unit     ! fortran unit attached to file
    character(len=*),intent(in) :: group    ! namelist group name
    integer,         intent(out):: status   ! 0 for success, -1 if group name not found
    !
    ! Local Values
    !---------------
    integer           :: len_grp
    integer           :: ios    ! io status
    character(len=80) :: inrec  ! first 80 characters of input record
    character(len=80) :: inrec2 ! left adjusted input record

    len_grp = len_trim(group)

    ios = 0
    do while (ios <= 0)
      read(unit,'(a)',iostat=ios,end=100) inrec
      if(ios <= 0) then  
        ! ios < 0  indicates an end of record condition
        ! look for group name in this record
        !---------------------------------------------

        ! remove leading blanks then check for leading '&'
        !-------------------------------------------------
        inrec2 = adjustl(inrec)
        if(inrec2(1:1) == '&') then
          ! check for group name
          !-------------------------------
          if(trim(group) == inrec2(2:len_grp+1)) then
            ! found group name.  backspace to leave file position 
            ! at this record
            !-------------------------------------------------------
            backspace(unit)
            status = 0
            return
          endif
        endif
      endif
    end do

100 continue  ! end of file processing
    status = -1

    ! End Routine
    !---------------
    return
  end subroutine find_group_name
  !==========================================================================

end module SE_Options

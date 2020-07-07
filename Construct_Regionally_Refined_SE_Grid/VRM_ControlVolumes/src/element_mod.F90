module element_mod

  ! Useful modules
  !------------------
  use SE_Constants          ,only: real_kind,long_kind,int_kind,np,npsq,nlev,nlevp,qsize_d
  use coordinate_systems_mod,only: spherical_polar_t, cartesian2D_t, cartesian3D_t, distance
  use edge_mod,              only: edgedescriptor_t, rotation_t
  use gridgraph_mod,         only: gridvertex_t

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private

  public :: timelevels
  public :: StateComponents

  public :: elem_state_t
  public :: derived_state_t
  public :: elem_accum_t
  public :: index_t
  public :: element_t

  public :: GetColumnIdP
  public :: GetColumnIdV
  public :: element_coordinates
  public :: element_var_coordinates
  public :: element_var_coordinates3D
  public :: allocate_element_desc

  ! PFC: Hack to set global option values
  !--------------------------------------
  public :: init_element_mod
  private:: max_neigh_edges
  integer:: max_neigh_edges

  ! Global Parameters
  !-------------------
  integer               ,parameter:: timelevels      = 3
  integer(kind=int_kind),parameter:: StateComponents = 8 ! num prognistics variables 
                                                         ! (for prim_restart_mod.F90)

  ! Type Definitions
  !------------------
  !=====================================================================
  ! =========== PRIMITIVE-EQUATION DATA-STRUCTURES =====================
  !=====================================================================
  type elem_state_t
    ! prognostic variables for preqx solver
    !
    ! prognostics must match those in prim_restart_mod.F90
    ! vertically-lagrangian code advects dp3d instead of ps_v
    ! tracers Q, Qdp always use 2 level time scheme
    !----------------------------------------------------------------
    real(kind=real_kind):: v   (np,np,2,nlev,timelevels) ! velocity                         1
    real(kind=real_kind):: T   (np,np,nlev,timelevels)   ! temperature                      2
    real(kind=real_kind):: dp3d(np,np,nlev,timelevels)   ! delta p on levels                8
    real(kind=real_kind):: lnps(np,np,timelevels)        ! log surface pressure             3
    real(kind=real_kind):: ps_v(np,np,timelevels)        ! surface pressure                 4
    real(kind=real_kind):: phis(np,np)                   ! surface geopotential(prescribed) 5
    real(kind=real_kind):: Q   (np,np,nlev,qsize_d)      ! Tracer concentration             6
    real(kind=real_kind):: Qdp (np,np,nlev,qsize_d,2)    ! Tracer mass                      7
  end type elem_state_t

  type derived_state_t
    ! diagnostic variables for preqx solver
    !
    ! storage for subcycling tracers/dynamics
    ! if (compute_mean_flux==1) vn0=time_avg(U*dp) else vn0=U at tracer-time t
    !--------------------------------------------------------------------------
    real(kind=real_kind):: vn0(np,np,2,nlev)             ! velocity for SE tracer advection
    real(kind=real_kind):: vstar(np,np,2,nlev)           ! velocity on Lagrangian surfaces
    real(kind=real_kind):: dpdiss_biharmonic(np,np,nlev) ! mean dp dissipation tendency,if nu_p>0
    real(kind=real_kind):: dpdiss_ave(np,np,nlev)        ! mean dp used to compute psdiss_tens
    ! diagnostics for explicit timestep
    !-------------------------------------------------------
    real(kind=real_kind):: phi(np,np,nlev)               ! geopotential
    real(kind=real_kind):: omega_p(np,np,nlev)           ! vertical tendency (derived)       
    real(kind=real_kind):: eta_dot_dpdn(np,np,nlevp)     ! mean vertical flux from dynamics
    ! semi-implicit diagnostics: computed in explict-component, reused in Helmholtz-component.
    !-------------------------------------------------------
    real(kind=real_kind):: grad_lnps(np,np,2)            ! gradient of log surface pressure
    real(kind=real_kind):: zeta(np,np,nlev)              ! relative vorticity
    real(kind=real_kind):: div(np,np,nlev,timelevels)    ! divergence
    ! tracer advection fields used for consistency and limiters
    !-------------------------------------------------------
    real(kind=real_kind):: dp(np,np,nlev)                ! for dp_tracers at physics timestep
    real(kind=real_kind):: divdp(np,np,nlev)             ! divergence of dp
    real(kind=real_kind):: divdp_proj(np,np,nlev)        ! DSSed divdp
    ! forcing terms for CAM
    !-------------------------------------------------------
    real(kind=real_kind):: FQ(np,np,nlev,qsize_d, 1)     ! tracer forcing  
    real(kind=real_kind):: FM(np,np,2,nlev, 1)           ! momentum forcing
    real(kind=real_kind):: FT(np,np,nlev, 1)             ! temperature forcing
    real(kind=real_kind):: omega_prescribed(np,np,nlev)  ! prescribed vertical tendency
    ! forcing terms for both CAM and HOMME
    ! FQps for conserving dry mass in the presence of precipitation
    !------------------------------------------------------------
    real(kind=real_kind):: pecnd(np,np,nlev)      ! pressure perturbation from condensate
    real(kind=real_kind):: FQps(np,np,timelevels) ! forcing of FQ on ps_v 
  end type derived_state_t

  type elem_accum_t
    ! the "4" timelevels represents data computed at:
    !  1  t-.5   
    !  2  t+.5   after dynamics
    !  3  t+.5   after forcing
    !  4  t+.5   after Robert
    ! after calling TimeLevelUpdate, all times above decrease by 1.0
    !----------------------------------------------------------------
    real(kind=real_kind):: KEner(np,np,4)       
    real(kind=real_kind):: PEner(np,np,4)       
    real(kind=real_kind):: IEner(np,np,4)    
    real(kind=real_kind):: IEner_wet(np,np,4)    
    real(kind=real_kind):: Qvar(np,np,qsize_d,4)     ! Q variance at half time levels   
    real(kind=real_kind):: Qmass(np,np,qsize_d,4)    ! Q mass at half time levels
    real(kind=real_kind):: Q1mass(np,np,qsize_d)     ! Q mass at full time levels
  end type elem_accum_t

  !=====================================================================
  ! ============= DATA-STRUCTURES COMMON TO ALL SOLVERS ================
  !=====================================================================
  type index_t
    integer(kind=int_kind):: ia(npsq),ja(npsq)
    integer(kind=int_kind):: is,ie
    integer(kind=int_kind):: NumUniquePts
    integer(kind=int_kind):: UniquePtOffset
  end type index_t

  type element_t
    integer(kind=int_kind) :: LocalId
    integer(kind=int_kind) :: GlobalId
    ! Coordinate values of element points
    !-------------------------------------------
    type(spherical_polar_t):: spherep(np,np)    ! Spherical coords of GLL points
    ! Equ-angular gnomonic projection coordinates
    !----------------------------------------------
    type(cartesian2D_t)    :: cartp(np,np)      ! gnomonic coords of GLL points 
    type(cartesian2D_t)    :: corners(4)        ! gnomonic coords of element corners
    real(kind=real_kind)   :: u2qmap(4,2)       ! bilinear map from ref element to quad in 
    ! 3D cartesian coordinates
    !----------------------------------------------
    type (cartesian3D_t)   :: corners3D(4)  
    ! Element diagnostics
    !----------------------------------------------
    real(kind=real_kind)   :: area               ! Area of element
    real(kind=real_kind)   :: max_eig            ! max singular value of metinv
    real(kind=real_kind)   :: min_eig            ! min singular value of metinv
    real(kind=real_kind)   :: max_eig_ratio      ! max ratio of singular values
    real(kind=real_kind)   :: dx_short           ! short length scale
    real(kind=real_kind)   :: dx_long            ! long length scale
    real(kind=real_kind)   :: variable_hyperviscosity(np,np) ! hyperviscosity based on above
    real(kind=real_kind)   :: hv_courant                     ! hyperviscosity courant number
    real(kind=real_kind)   :: tensorVisc(2,2,np,np)          !og, matrix V for tensor viscosity
    ! Edge connectivity information
    !----------------------------------------------
    integer(kind=int_kind) :: node_numbers(4)
    integer(kind=int_kind) :: node_multiplicity(4) ! number of elements sharing corner node
    type(GridVertex_t)     :: vertex               ! element grid vertex information
    type(EdgeDescriptor_t) :: desc
!PFC    type(elem_state_t)     :: state
!PFC    type(derived_state_t)  :: derived
!PFC    type(elem_accum_t)     :: accum
    ! Metric terms 
    !----------------------------------------------
    real(kind=real_kind)   :: met   (2,2,np,np) ! metric tensor on velocity and pressure grid
    real(kind=real_kind)   :: metinv(2,2,np,np) ! metric tensor on velocity and pressure grid
    real(kind=real_kind)   :: metdet    (np,np) ! g = SQRT(det(g_ij)) on velocity and pressure grid
    real(kind=real_kind)   :: rmetdet   (np,np) ! 1/metdet on velocity pressure grid
    real(kind=real_kind)   :: D     (2,2,np,np) ! Map covariant field on cube to vector field on sphere
    real(kind=real_kind)   :: Dinv  (2,2,np,np) ! Map vector field on the sphere to covariant v on cube
    ! Convert vector fields from spherical to rectangular components
    ! The transpose of this operation is its pseudoinverse.
    !------------------------------------------------------------
    real(kind=real_kind)   :: vec_sphere2cart(np,np,3,2)
    ! Mass matrix terms for an element on a cube face
    !------------------------------------------------------------
    real(kind=real_kind)   :: mp (np,np)    ! mass matrix on v and p grid
    real(kind=real_kind)   :: rmp(np,np)    ! inverse mass matrix on v and p grid
    ! Mass matrix terms for an element on the sphere
    ! This mass matrix is used when solving the equations in weak form
    ! with the natural (surface area of the sphere) inner product
    !------------------------------------------------------------
!+++ARH
    real(kind=real_kind)   :: tmp_elemarea(np,np)
!---ARH
    real(kind=real_kind)   :: spheremp(np,np)       ! mass matrix on v and p grid
    real(kind=real_kind)   :: rspheremp(np,np)      ! inverse mass matrix on v and p grid
    integer(kind=long_kind):: gdofP(np,np)          ! global degree of freedom (P-grid)
    real(kind=real_kind)   :: fcor(np,np)           ! Coreolis term
    type(index_t)          :: idxP
    type(index_t),pointer  :: idxV
    integer                :: FaceNum
    ! force element_t to be a multiple of 8 bytes.  
    ! on BGP, code will crash (signal 7, or signal 15) if 8 byte alignment is off
    ! check core file for:
    ! core.63:Generated by interrupt..(Alignment Exception DEAR=0xa1ef671c 
    !                                                      ESR=0x01800000 CCR0=0x4800a002)
    !------------------------------------------------------------
    integer :: dummy
  end type element_t


contains
  !==================================================================
  subroutine init_element_mod(I_SEopt)
    ! PFC: Hack to set global option values
    !--------------------------------------
    use SE_Options,only: SEoptions_t
    ! Passed Variables
    !------------------
    type(SEoptions_t),intent(in):: I_SEopt

    max_neigh_edges = I_SEopt%max_neigh_edges

    ! End Routine
    !-------------
    return
  end subroutine init_element_mod
  !==================================================================


  !==================================================================
  function GetColumnIdP(elem,ii,jj) result(col_id)
    ! GetColumnIdP: Get unique identifier for a Physics column on the P-grid
    !=======================================================================
    !
    ! Passed Variables
    !------------------
    type(element_t),intent(in):: elem
    integer        ,intent(in):: ii,jj
    integer                   :: col_id

    col_id = elem%gdofP(ii,jj)

    ! End Function
    !-------------
    return
  end function GetColumnIdP
  !==================================================================


  !==================================================================
  function GetColumnIdV(elem,ii,jj) result(col_id)
    ! GetColumnIdV: Get unique identifier for a Physics column on the V-grid
    !=========================================================================
    !
    ! Passed Variables
    !-------------------
    type(element_t),intent(in):: elem
    integer        ,intent(in):: ii,jj
    integer                   :: col_id
   
    col_id = elem%gdofP(ii,jj)

    ! End Function
    !-------------
    return
  end function GetColumnIdV
  !==================================================================


  !==================================================================
  function element_coordinates(start,end,points) result(cart)
    ! element_coordinates: Initialize 2D rectilinear element colocation points
    !=========================================================================
    use SE_Constants,only: longdouble_kind
    !
    ! Passed Variables
    !-------------------
    type(cartesian2D_t)       ,intent(in):: start
    type(cartesian2D_t)       ,intent(in):: end
    real(kind=longdouble_kind),intent(in):: points(:)
    !
    ! Local Values
    !--------------
    type(cartesian2D_t)       :: cart(size(points),size(points))
    type(cartesian2D_t)       :: length, centroid
    real(kind=longdouble_kind):: yy

    integer ii,jj

    length%x   = 0.50D0*(end%x-start%x)
    length%y   = 0.50D0*(end%y-start%y)
    centroid%x = 0.50D0*(end%x+start%x) 
    centroid%y = 0.50D0*(end%y+start%y) 
    do jj=1,size(points)
      yy = centroid%y + length%y*points(jj)
      do ii=1,size(points)
        cart(ii,jj)%x = centroid%x + length%x*points(ii)
        cart(ii,jj)%y = yy
      end do
    end do

    ! End Function
    !-------------
    return
  end function element_coordinates
  !==================================================================


  !==================================================================
  function element_var_coordinates(c,points) result(cart)
    ! element_var_coordinates:
    !
    !=======================================================
    use SE_Constants, only : longdouble_kind
    !
    ! Passed Variables
    !-----------------
    type(cartesian2D_t)       , intent(in):: c(4)
    real(kind=longdouble_kind), intent(in):: points(:)
    type(cartesian2D_t)                   :: cart(size(points),size(points))
    !
    ! Local Values
    !----------------
    real(kind=longdouble_kind):: p(size(points))
    real(kind=longdouble_kind):: q(size(points))

    integer ii,jj

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do jj=1,size(points)
    do ii=1,size(points)
      cart(ii,jj)%x = p(ii)*p(jj)*c(1)%x &
                    + q(ii)*p(jj)*c(2)%x &
                    + q(ii)*q(jj)*c(3)%x &
                    + p(ii)*q(jj)*c(4)%x 
      cart(ii,jj)%y = p(ii)*p(jj)*c(1)%y &
                    + q(ii)*p(jj)*c(2)%y &
                    + q(ii)*q(jj)*c(3)%y &
                    + p(ii)*q(jj)*c(4)%y 
    end do
    end do

    ! End Function
    !-------------
    return
  end function element_var_coordinates
  !==================================================================


  !==================================================================
  function element_var_coordinates3D(c,points) result(cart)
    ! element_var_coordinates3D:
    !
    !=======================================================
    use SE_Constants, only : longdouble_kind
    !
    ! Passed Values
    !-----------------
    type(cartesian3D_t)       ,intent(in):: c(4)
    real(kind=longdouble_kind),intent(in):: points(:)
    type(cartesian3D_t)                  :: cart(size(points),size(points))
    !
    ! Local Values
    !--------------
    real(kind=longdouble_kind):: p(size(points))
    real(kind=longdouble_kind):: q(size(points)),rr

    integer ii,jj

    p(:) = (1.0D0-points(:))/2.0D0
    q(:) = (1.0D0+points(:))/2.0D0

    do jj=1,size(points)
    do ii=1,size(points)
      cart(ii,jj)%x = p(ii)*p(jj)*c(1)%x &
                    + q(ii)*p(jj)*c(2)%x &
                    + q(ii)*q(jj)*c(3)%x &
                    + p(ii)*q(jj)*c(4)%x 
      cart(ii,jj)%y = p(ii)*p(jj)*c(1)%y &
                    + q(ii)*p(jj)*c(2)%y &
                    + q(ii)*q(jj)*c(3)%y &
                    + p(ii)*q(jj)*c(4)%y 
      cart(ii,jj)%z = p(ii)*p(jj)*c(1)%z &
                    + q(ii)*p(jj)*c(2)%z &
                    + q(ii)*q(jj)*c(3)%z &
                    + p(ii)*q(jj)*c(4)%z 

      ! project back to sphere:
      !-------------------------
      rr            = distance(cart(ii,jj))
      cart(ii,jj)%x = cart(ii,jj)%x/rr
      cart(ii,jj)%y = cart(ii,jj)%y/rr
      cart(ii,jj)%z = cart(ii,jj)%z/rr
    end do
    end do

    ! End Function
    !-------------
    return
  end function element_var_coordinates3D
  !==================================================================


  !==================================================================
  subroutine allocate_element_desc(elem)
    ! allocate_element_desc:
    !
    !========================================================
    ! 
    ! Passed Variables
    !------------------
    type(element_t),intent(inout):: elem(:)
    !
    ! Local Values
    !-------------
    integer num,ii,jj

    num = size(elem)
    
    do jj=1,num
      allocate(elem(jj)%desc%putmapP      (max_neigh_edges))
      allocate(elem(jj)%desc%getmapP      (max_neigh_edges))
      allocate(elem(jj)%desc%putmapP_ghost(max_neigh_edges))
      allocate(elem(jj)%desc%getmapP_ghost(max_neigh_edges))
      allocate(elem(jj)%desc%reverse      (max_neigh_edges))
      allocate(elem(jj)%desc%globalID     (max_neigh_edges))
      allocate(elem(jj)%desc%loc2buf      (max_neigh_edges))
      do ii=1,max_neigh_edges
        elem(jj)%desc%loc2buf (ii)=ii
        elem(jj)%desc%globalID(ii)=-1
      end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine allocate_element_desc
  !==================================================================

end module element_mod

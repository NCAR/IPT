#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _BEGIN_FACE 1
#define _END_FACE   4
#undef _FACE_6
#undef _FACE_5

module cube_mod

  ! Useful modules
  !----------------
  use SE_Constants          ,only: real_kind, long_kind, longdouble_kind
  use SE_Constants          ,only: DD_PI, rearth
  use coordinate_systems_mod,only: spherical_polar_t,cartesian3D_t   ,cartesian2d_t    , &
                                   projectpoint     ,cubedsphere2cart,spherical_to_cart, &
                                   sphere_tri_area  ,dist_threshold  ,change_coordinates
  use err_exit              ,only: endrun

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private

  public :: nfaces
  public :: nInnerElemEdge
  public :: nCornerElemEdge
  public :: face_t
  public :: cube_face_coord_t
  public :: rotate_grid
  public :: cube_init_atomic
  private:: coordinates_atomic
  private:: elem_jacobians
  private:: metric_atomic
  public :: covariant_rot
  public :: contravariant_rot
  public :: dmap
  private:: dmap_equiangular
  public :: vmap
  private:: dmap_elementlocal
  private:: coreolis_init_atomic
  public :: rotation_init_atomic
  public :: set_corner_coordinates
  public :: assign_node_numbers_to_elem
  public :: convert_gbl_index
  public :: CubeTopology
!PFC   public :: cube_assemble
  public :: CubeEdgeCount
  public :: CubeElemCount
  public :: CubeSetupEdgeIndex
  public :: ref2sphere
  private:: ref2sphere_double
  private:: ref2sphere_longdouble
  private:: ref2sphere_equiangular_double
  private:: ref2sphere_equiangular_longdouble
  private:: ref2sphere_elementlocal_double
  private:: ref2sphere_elementlocal_longdouble
  private:: ref2sphere_elementlocal_q

  ! PFC: Hack to set global option values
  !--------------------------------------
  public :: init_cube_mod
  private:: hypervis_scaling
  private:: cubed_sphere_map
  private:: ne
  private:: partmethod
  real(kind=real_kind):: hypervis_scaling
  integer             :: cubed_sphere_map
  integer             :: ne
  integer             :: partmethod



  ! Parameter Values
  !------------------
  integer,parameter:: nfaces         = 6  ! number of faces on the cube
  integer,parameter:: nInnerElemEdge = 8  ! number of edges for an interior element
  integer,parameter:: nCornerElemEdge= 4  ! number of corner elements

  real(kind=real_kind),parameter:: cube_xstart = -0.25D0*DD_PI
  real(kind=real_kind),parameter:: cube_xend   =  0.25D0*DD_PI
  real(kind=real_kind),parameter:: cube_ystart = -0.25D0*DD_PI
  real(kind=real_kind),parameter:: cube_yend   =  0.25D0*DD_PI

  ! Type Definitions
  !--------------------
  type face_t
    type (spherical_polar_t):: sphere0       ! tangent point of face on sphere
    type (spherical_polar_t):: sw            ! sw corner of face on sphere
    type (spherical_polar_t):: se            ! se corner of face on sphere
    type (spherical_polar_t):: ne            ! ne corner of face on sphere
    type (spherical_polar_t):: nw            ! nw corner of face on sphere
    type (cartesian3D_t)    :: P0
    type (cartesian3D_t)    :: X0
    type (cartesian3D_t)    :: Y0
    integer                 :: number
    integer                 :: padding       ! pad the struct
  end type face_t

  type cube_face_coord_t
    real(real_kind)      :: x        ! x coordinate
    real(real_kind)      :: y        ! y coordinate
    type(face_t), pointer:: face     ! face
  end type cube_face_coord_t

  ! Global Data
  !---------------------
  real(kind=real_kind):: rotate_grid = 0  ! Rotate the North Pole:  used for JW baroclinic 
                                          ! test case Setting this only changes Coriolis.  
                                          ! User must also rotate initial condition

  ! Public interface to REFERECE element map
  !-----------------------------------------
#if HOMME_QUAD_PREC
  interface ref2sphere
    module procedure ref2sphere_double
    module procedure ref2sphere_longdouble
  end interface
#else
  ! both routines have identical arguments in this case, cant use interface
  interface ref2sphere
    module procedure ref2sphere_double
  end interface
#endif


contains
  !==================================================================
  subroutine init_cube_mod(I_SEopt)
    ! PFC: Hack to set global option values
    !--------------------------------------
    use SE_Options,only: SEoptions_t
    ! Passed Variables
    !------------------
    type(SEoptions_t),intent(in):: I_SEopt

  print *,' PFC: init_cube_mod: ne=',ne,I_SEopt%ne
    hypervis_scaling = I_SEopt%hypervis_scaling
    cubed_sphere_map = I_SEopt%cubed_sphere_map
    ne               = I_SEopt%ne
    partmethod       = I_SEopt%partmethod
  print *,' PFC: init_cube_mod: ne=',ne,I_SEopt%ne

    ! End Routine
    !-------------
    return
  end subroutine init_cube_mod
  !==================================================================


  !==================================================================
  subroutine cube_init_atomic(elem,gll_points,alpha_in)
    !  cube_init_atomic:  Initialize element descriptors for 
    !                     cube sphere case for each element ... 
    !============================================================
    use element_mod ,only: element_t
    use SE_Constants,only: np
    !
    ! Passed Variables
    !--------------------
    type(element_t)           ,intent(inout):: elem
    real(kind=longdouble_kind)              :: gll_points(np)
    real(kind=real_kind)      ,optional     :: alpha_in
    !
    ! Local Values
    !--------------
    real(kind=real_kind):: alpha=1

    ! Optionally change alpha
    !--------------------------
    if(present(alpha_in)) alpha=alpha_in
    
    ! Inialize values
    !-------------------
    elem%FaceNum           = elem%vertex%face_number
    elem%desc%use_rotation = 0

    call coordinates_atomic  (elem,gll_points)
    call metric_atomic       (elem,gll_points,alpha)
    call coreolis_init_atomic(elem)
!    call solver_weights_atomic(elem)

    ! End Routine
    !-------------
    return
  end subroutine cube_init_atomic
  !==================================================================


  !==================================================================
  subroutine coordinates_atomic(elem,gll_points)
    ! coordinates_atomic: Initialize element coordinates for 
    !                     cube-sphere case ... (atomic) 
    !===============================================================================
    use element_mod ,only: element_t, element_var_coordinates
    use SE_Constants,only: np
    !
    ! Passed Variables
    !-------------------
    type(element_t)           :: elem
    real(kind=longdouble_kind):: gll_points(np)
    !
    ! Local Values
    !--------------
    real(kind=real_kind):: area1,area2
    type(cartesian3d_t) :: quad(4)
    integer             :: face_no,ii,jj

    ! compute the corners in Cartesian coordinates
    !---------------------------------------------
    face_no = elem%vertex%face_number
    do ii=1,4
      elem%corners3D(ii)=cubedsphere2cart(elem%corners(ii),face_no)
    end do

    ! compute lat/lon coordinates of each GLL point
    !-----------------------------------------------
    do ii=1,np
    do jj=1,np
      elem%spherep(ii,jj)=ref2sphere(gll_points(ii),gll_points(jj),elem)
    end do
    end do

    ! compute the [-pi/2,pi/2] cubed sphere coordinates:
    !----------------------------------------------------------
    elem%cartp=element_var_coordinates(elem%corners,gll_points)

    ! Matrix describing vector conversion to cartesian Zonal direction
    !------------------------------------------------------------------
    elem%vec_sphere2cart(:,:,1,1) = -sin(elem%spherep(:,:)%lon)
    elem%vec_sphere2cart(:,:,2,1) =  cos(elem%spherep(:,:)%lon)
    elem%vec_sphere2cart(:,:,3,1) =  0.0_real_kind

    ! Meridional direction
    !----------------------
    elem%vec_sphere2cart(:,:,1,2) = -sin(elem%spherep(:,:)%lat)*cos(elem%spherep(:,:)%lon)
    elem%vec_sphere2cart(:,:,2,2) = -sin(elem%spherep(:,:)%lat)*sin(elem%spherep(:,:)%lon)
    elem%vec_sphere2cart(:,:,3,2) =  cos(elem%spherep(:,:)%lat)

    ! End Routine
    !-------------
    return
  end subroutine coordinates_atomic
  !==================================================================


  !==================================================================
  subroutine elem_jacobians(coords, unif2quadmap)
    ! elem_jacobians: Calculate Jacobian associated with mapping from arbitrary 
    !                 quadrilateral to [-1,1]^2 along with its inverse and determinant
    !=================================================================================
    use SE_Constants, only : np
    !
    ! Passed Variables
    !-------------------
    type(cartesian2D_t) ,intent(in) :: coords(np,np)
    real(kind=real_kind),intent(out):: unif2quadmap(4,2) ! unif2quadmap=bilinear map from 
                                                         ! [-1,1]^2 -> arbitrary quadrilateral
    
    unif2quadmap(1,1)=( coords(1,1)%x+coords(np,1)%x+coords(np,np)%x+coords(1,np)%x)/4.0d0
    unif2quadmap(1,2)=( coords(1,1)%y+coords(np,1)%y+coords(np,np)%y+coords(1,np)%y)/4.0d0
    unif2quadmap(2,1)=(-coords(1,1)%x+coords(np,1)%x+coords(np,np)%x-coords(1,np)%x)/4.0d0
    unif2quadmap(2,2)=(-coords(1,1)%y+coords(np,1)%y+coords(np,np)%y-coords(1,np)%y)/4.0d0
    unif2quadmap(3,1)=(-coords(1,1)%x-coords(np,1)%x+coords(np,np)%x+coords(1,np)%x)/4.0d0
    unif2quadmap(3,2)=(-coords(1,1)%y-coords(np,1)%y+coords(np,np)%y+coords(1,np)%y)/4.0d0
    unif2quadmap(4,1)=( coords(1,1)%x-coords(np,1)%x+coords(np,np)%x-coords(1,np)%x)/4.0d0
    unif2quadmap(4,2)=( coords(1,1)%y-coords(np,1)%y+coords(np,np)%y-coords(1,np)%y)/4.0d0

    ! End Routine
    !-------------
    return
  end subroutine elem_jacobians
  !==================================================================


  !==================================================================
  subroutine metric_atomic(elem,gll_points,alpha)
    ! metric_atomic: Initialize cube-sphere metric terms:
    !
    ! equal angular elements (atomic)
    ! initialize:  
    !         metdet, rmetdet  (analytic)    = detD, 1/detD
    !         met                (analytic)    D^t D     (symmetric)
    !         metdet             (analytic)    = detD
    !         metinv             (analytic)    Dinv Dinv^t  (symmetic)
    !         D     (from subroutine vmap)
    !         Dinv  (computed directly from D)
    ! 
    ! ucontra = Dinv * u  =  metinv * ucov   
    ! ucov    = D^t * u   =  met * ucontra
    !
    ! we also compute DE = D*E, where 
    ! E = eigenvectors of metinv as a basis      metinv = E LAMBDA E^t
    !   
    ! ueig = E^t ucov  = E^t D^t u =  (DE)^t u  
    !  
    !
    ! so if we want to tweak the mapping by a factor alpha (so he weights add up to 4pi, 
    ! for example) we take:
    !    NEW       OLD     
    !       D = sqrt(alpha) D  and then rederive all quantities.  
    !    detD = alpha detD
    !    
    ! where alpha = 4pi/SEMarea, SEMarea = global sum elem(ie)%mv(i,j)*elem(ie)%metdet(i,j)
    ! 
    !==============================================================================
    use element_mod ,only: element_t
    use SE_Constants,only: np
    use SE_Constants,only: rrearth
    !
    ! Passed Variables
    !-------------------
    type(element_t)           :: elem
    real(kind=longdouble_kind):: gll_points(np)
    real(kind=real_kind)      :: alpha
    !
    ! Local Values
    !-------------
    real(kind=real_kind) :: rr         ! distance from origin for point on cube tangent 
                                       ! to unit sphere
    real(kind=real_kind) :: const
    real(kind=real_kind) :: norm
    real(kind=real_kind) :: detD       ! determinant of vector field mapping matrix.  
    real(kind=real_kind) :: x1         ! 1st cube face coordinate
    real(kind=real_kind) :: x2         ! 2nd cube face coordinate
    real(kind=real_kind) :: tmpD(2,2)
    real(kind=real_kind) :: MM(2,2),EE(2,2),eig(2),DE(2,2),DEL(2,2),VV(2,2)
    real(kind=real_kind) :: nu1, nu2, lamStar1, lamStar2
    real (kind=real_kind):: l1, l2, sc           ! eigen values of met
    real (kind=real_kind):: roundoff_err = 1e-11 !!! OG: this is a temporal fix
    integer              :: imaxM(2)
    integer              :: ii,jj
    
    !===================================================================================
    ! Initialize differential mapping operator to and from vector fields on the sphere to 
    ! contravariant vector fields on the cube
    !    i.e. dM/dx^i in Sadourney (1972) and it's inverse
    !===================================================================================

    ! MNL: Calculate Jacobians of bilinear map from cubed-sphere to ref element
    !--------------------------------------------------------------------------
    if(cubed_sphere_map==0) then
      call elem_jacobians(elem%cartp, elem%u2qmap)
    endif

    elem%max_eig      = 0.0d0
    elem%min_eig      = 1d99
    elem%max_eig_ratio= 0d0
    do jj=1,np
    do ii=1,np
      x1=gll_points(ii)
      x2=gll_points(jj)
      call dmap(elem%D(:,:,ii,jj),elem,x1,x2)

      ! Numerical metric tensor based on analytic D: met = D^T times D
      ! (D maps between sphere and reference element)
      !-------------------------------------------------------------
      elem%met(1,1,ii,jj) = elem%D(1,1,ii,jj)*elem%D(1,1,ii,jj) + &
                            elem%D(2,1,ii,jj)*elem%D(2,1,ii,jj)
      elem%met(1,2,ii,jj) = elem%D(1,1,ii,jj)*elem%D(1,2,ii,jj) + &
                            elem%D(2,1,ii,jj)*elem%D(2,2,ii,jj)
      elem%met(2,1,ii,jj) = elem%D(1,1,ii,jj)*elem%D(1,2,ii,jj) + &
                            elem%D(2,1,ii,jj)*elem%D(2,2,ii,jj)
      elem%met(2,2,ii,jj) = elem%D(1,2,ii,jj)*elem%D(1,2,ii,jj) + &
                            elem%D(2,2,ii,jj)*elem%D(2,2,ii,jj)

      ! compute D^-1...
      ! compute determinant of D mapping matrix... if not zero compute inverse
      !------------------------------------------------------------------------
      detD = elem%D(1,1,ii,jj)*elem%D(2,2,ii,jj) - elem%D(1,2,ii,jj)*elem%D(2,1,ii,jj)      

      elem%Dinv(1,1,ii,jj) =  elem%D(2,2,ii,jj)/detD
      elem%Dinv(1,2,ii,jj) = -elem%D(1,2,ii,jj)/detD
      elem%Dinv(2,1,ii,jj) = -elem%D(2,1,ii,jj)/detD
      elem%Dinv(2,2,ii,jj) =  elem%D(1,1,ii,jj)/detD

      ! L2 norm = sqrt max eigenvalue of metinv
      !         = 1/sqrt(min eigenvalue of met)
      ! l1 and l2 are eigenvalues of met
      ! (should both be positive, l1 > l2)
      !---------------------------------------------
      l1 = (elem%met(1,1,ii,jj) + elem%met(2,2,ii,jj)                      &
            + sqrt(4.0d0*elem%met(1,2,ii,jj)*elem%met(2,1,ii,jj)           &
                   +    (elem%met(1,1,ii,jj)-elem%met(2,2,ii,jj))**2))/2.0d0
      l2 = (elem%met(1,1,ii,jj) + elem%met(2,2,ii,jj)                      &
            - sqrt(4.0d0*elem%met(1,2,ii,jj)*elem%met(2,1,ii,jj)           &
                   +    (elem%met(1,1,ii,jj)-elem%met(2,2,ii,jj))**2))/2.0d0

      ! Max L2 norm of Dinv is sqrt of max eigenvalue of metinv
      ! max eigenvalue of metinv is 1/min eigenvalue of met
      !---------------------------------------------------------
      norm         = 1.0d0/sqrt(min(abs(l1),abs(l2)))
      elem%max_eig = max(norm, elem%max_eig)

      ! Min L2 norm of Dinv is sqrt of min eigenvalue of metinv
      ! min eigenvalue of metinv is 1/max eigenvalue of met
      !---------------------------------------------------------
      norm         = 1.0d0/sqrt(max(abs(l1),abs(l2)))
      elem%min_eig = min(norm, elem%min_eig)

      ! Need inverse of met if not calculated analytically
      !----------------------------------------------------
      elem%metdet (ii,jj) = abs(detD)
      elem%rmetdet(ii,jj) = 1.0D0/abs(detD)

      elem%metinv(1,1,ii,jj) =  elem%met(2,2,ii,jj)/(detD*detD)
      elem%metinv(1,2,ii,jj) = -elem%met(1,2,ii,jj)/(detD*detD)
      elem%metinv(2,1,ii,jj) = -elem%met(2,1,ii,jj)/(detD*detD)
      elem%metinv(2,2,ii,jj) =  elem%met(1,1,ii,jj)/(detD*detD)

      ! matricies for tensor hyper-viscosity
      ! compute eigenvectors of metinv (probably same as computed above)
      !-----------------------------------------------------------------
      MM(:,:) = elem%metinv(:,:,ii,jj)

      eig(1) = (MM(1,1) + MM(2,2)                                   &
                + sqrt(4.0d0*MM(1,2)*MM(2,1) + (MM(1,1)-MM(2,2))**2))/2.0d0
      eig(2) = (MM(1,1) + MM(2,2)                                   &
                - sqrt(4.0d0*MM(1,2)*MM(2,1) + (MM(1,1)-MM(2,2))**2))/2.0d0
          
      ! use DE to store M - Lambda, to compute eigenvectors
      !----------------------------------------------------
      DE(:,:)=MM(:,:)
      DE(1,1)=DE(1,1)-eig(1)
      DE(2,2)=DE(2,2)-eig(1)

      imaxM = maxloc(abs(DE))
      if(maxval(abs(DE))==0) then
        EE(1,1)=1.0d0
        EE(2,1)=0.0d0
      elseif((imaxM(1)==1).and.(imaxM(2)==1)) then
        EE(2,1)=1.0d0
        EE(1,1)=(-DE(2,1)/DE(1,1))
      elseif((imaxM(1)==1).and.(imaxM(2)==2)) then
        EE(2,1)=1.0d0
        EE(1,1)=(-DE(2,2)/DE(1,2))
      elseif((imaxM(1)==2).and.(imaxM(2)==1)) then
        EE(1,1)=1.0d0
        EE(2,1)=(-DE(1,1)/DE(2,1))
      elseif((imaxM(1)==2).and.(imaxM(2)==2)) then
        EE(1,1)=1.0d0
        EE(2,1)=(-DE(1,2)/DE(2,2))
      else
        call endrun('Impossible error in cube_mod.F90::metric_atomic()')
      endif

      ! the other eigenvector is orthgonal:
      !------------------------------------
      EE(1,2)=-EE(2,1)
      EE(2,2)= EE(1,1)

      !normalize columns
      !--------------------
      EE(:,1)=EE(:,1)/sqrt(sum(EE(:,1)*EE(:,1)))
      EE(:,2)=EE(:,2)/sqrt(sum(EE(:,2)*EE(:,2)))

      !============================================================================
      ! OBTAINING TENSOR FOR HV:
      !
      ! Instead of the traditional scalar Laplace operator \grad \cdot \grad
      ! we introduce \grad \cdot V \grad
      ! where V = D E LAM LAM^* E^T D^T. 
      ! Recall (metric_tensor)^{-1}=(D^T D)^{-1} = E LAM E^T.
      ! Here, LAM=diag( 4/((np-1)dx)^2 , 4/((np-1)dy)^2 ) = diag( 4/(dx_elem)^2, 4/(dy_elem)^2 )
      ! Note that metric tensors and LAM correspondingly are quantities on a unit sphere.
      !
      ! This motivates us to use V = D E LAM LAM^* E^T D^T
      ! where LAM^* = diag( nu1, nu2 ) where nu1, nu2 are HV coefficients scaled 
      !                                like (dx)^{hv_scaling/2}, (dy)^{hv_scaling/2}.
      ! (Halves in powers come from the fact that HV consists of two Laplace iterations.)
      !
      ! Originally, we took LAM^* = diag(
      !  1/(eig(1)**(hypervis_scaling/4.0d0))*(rearth**(hypervis_scaling/2.0d0))
      !  1/(eig(2)**(hypervis_scaling/4.0d0))*(rearth**(hypervis_scaling/2.0d0)) ) = 
      !  = diag( lamStar1, lamStar2)
      !  \simeq ((np-1)*dx_sphere / 2 )^hv_scaling/2 = SQRT(OPERATOR_HV)
      ! because 1/eig(...) \simeq (dx_on_unit_sphere)^2 .
      ! Introducing the notation OPERATOR = lamStar^2 is useful for conversion formulas.
      !
      ! This leads to the following conversion formula: nu_const is nu used for 
      ! traditional HV on uniform grids
      ! nu_tensor = nu_const * OPERATOR_HV^{-1}, so
      ! nu_tensor = nu_const *((np-1)*dx_sphere / 2 )^{ - hv_scaling} or
      ! nu_tensor = nu_const *(2/( (np-1) * dx_sphere) )^{hv_scaling} .
      ! dx_sphere = 2\pi *rearth/(np-1)/4/NE
      ! [nu_tensor] = [meter]^{4-hp_scaling}/[sec]
      !
      ! (1) Later developments:
      ! Apply tensor V only at the second Laplace iteration. 
      ! Thus, LAM^* should be scaled as (dx)^{hv_scaling}, (dy)^{hv_scaling},
      ! see this code below:
      !          DEL(1:2,1) = (lamStar1**2) *eig(1)*DE(1:2,1)
      !          DEL(1:2,2) = (lamStar2**2) *eig(2)*DE(1:2,2)
      !
      ! (2) Later developments:
      ! Bringing [nu_tensor] to 1/[sec]:
      !	  lamStar1=1/(eig(1)**(hypervis_scaling/4.0d0)) *(rearth**2.0d0)
      !	  lamStar2=1/(eig(2)**(hypervis_scaling/4.0d0)) *(rearth**2.0d0)
      ! OPERATOR_HV = ( (np-1)*dx_unif_sphere / 2 )^{hv_scaling} * rearth^4
      ! Conversion formula:
      ! nu_tensor = nu_const * OPERATOR_HV^{-1}, so
      ! nu_tensor = nu_const *( 2*rearth /((np-1)*dx))^{hv_scaling} * rearth^{-4.0}.
      !
      ! For the baseline coefficient nu=1e15 for NE30, 
      ! nu_tensor=7e-8 (BUT RUN TWICE AS SMALL VALUE FOR NOW) for hv_scaling=3.2
      ! and 
      ! nu_tensor=1.3e-6 for hv_scaling=4.0.
      !============================================================================

      ! matrix D*E
      !----------------
      DE(1,1)=sum(elem%D(1,:,ii,jj)*EE(:,1))
      DE(1,2)=sum(elem%D(1,:,ii,jj)*EE(:,2))
      DE(2,1)=sum(elem%D(2,:,ii,jj)*EE(:,1))
      DE(2,2)=sum(elem%D(2,:,ii,jj)*EE(:,2))

      lamStar1=1/(eig(1)**(hypervis_scaling/4.0d0)) *(rearth**2.0d0)
      lamStar2=1/(eig(2)**(hypervis_scaling/4.0d0)) *(rearth**2.0d0)

      ! matrix (DE) * Lam^* * Lam , tensor HV when V is applied at each Laplace calculation
      !          DEL(1:2,1) = lamStar1*eig(1)*DE(1:2,1)
      !          DEL(1:2,2) = lamStar2*eig(2)*DE(1:2,2)
      !
      ! matrix (DE) * (Lam^*)^2 * Lam, tensor HV when V is applied only once, at the last 
      ! Laplace calculation will only work with hyperviscosity, not viscosity
      !------------------------------------------------------------------------------------
      DEL(1:2,1) = (lamStar1**2) *eig(1)*DE(1:2,1)
      DEL(1:2,2) = (lamStar2**2) *eig(2)*DE(1:2,2)

      ! matrix (DE) * Lam^* * Lam  *E^t *D^t or (DE) * (Lam^*)^2 * Lam  *E^t *D^t 
      !---------------------------------------------------------------------------
      VV(1,1)=sum(DEL(1,:)*DE(1,:))
      VV(1,2)=sum(DEL(1,:)*DE(2,:))
      VV(2,1)=sum(DEL(2,:)*DE(1,:))
      VV(2,2)=sum(DEL(2,:)*DE(2,:))

      elem%tensorVisc(:,:,ii,jj)=VV(:,:)

    end do ! ii=1,np
    end do ! jj=1,np

    elem%dx_short = 1.0d0/(elem%max_eig*0.5d0*dble(np-1)*rrearth*1000.0d0)
    elem%dx_long  = 1.0d0/(elem%min_eig*0.5d0*dble(np-1)*rrearth*1000.0d0)

    !===============================================
    ! Initialize equal angular metric tensor on each 
    ! on velocity grid for unit sphere.
    !
    ! Initialize gdet = SQRT(ABS(DET(gij)))
    !
    ! These quantities are the same on every face
    ! of the cube.
    !=================================================

    ! mt: better might be to compute all these quantities 
    !     directly from D for consistency?
    !
    ! MNL: done
    !------------------------------------------------------
    elem%D      = elem%D      *sqrt(alpha) 
    elem%Dinv   = elem%Dinv   /sqrt(alpha) 
    elem%metdet = elem%metdet *alpha
    elem%rmetdet= elem%rmetdet/alpha
    elem%met    = elem%met    *alpha
    elem%metinv = elem%metinv /alpha

    ! End Routine
    !-------------
    return
  end subroutine metric_atomic
  !==================================================================


#if 0
  !==================================================================
  subroutine solver_weights_atomic(elem)
    ! solver_weights: For nonstaggered GaussLobatto elements, compute weights for 
    !                 redundant points in cg solver.
    !============================================================================
    use element_mod ,only: element_t
    use SE_Constants,only: np
    !
    ! Passed Variables
    !--------------------
    type(element_t)     :: elem
    !
    ! Local Values
    !--------------
    real(kind=real_kind):: xx
    integer             :: ii, jj

    ! compute cube face coordinates of element
    !-----------------------------------------
    do ii=1,np
    do jj=1,np
      if(ii==1) then
        if(jj==1) then
          xx = 1.0_real_kind/elem%node_multiplicity(1)
        elseif(jj==np) then
          xx = 1.0_real_kind/elem%node_multiplicity(4)
        else
          xx = 0.5_real_kind
        endif
      elseif(ii==np) then
        if(jj==1) then
          xx = 1.0_real_kind/elem%node_multiplicity(2)
        elseif(jj==np) then
          xx = 1.0_real_kind/elem%node_multiplicity(3)
        else
          xx = 0.5_real_kind
        endif
      elseif((jj==1).or.(jj==np)) then
        xx = 0.5_real_kind
      else
        xx = 1.0_real_kind
      endif
      elem%solver_wts(ii,jj) = xx
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine solver_weights_atomic
  !==================================================================
#endif


#if 1
  !==================================================================
  function covariant_rot(Da,Db) result(RR)
    ! covariant_rot: 2 x 2 matrix multiply:  Db^T * Da^-T
    !                for edge rotations: maps face a to face b
    !===========================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind):: Da(2,2)
    real(kind=real_kind):: Db(2,2)
    real(kind=real_kind):: RR(2,2)
    !
    ! Local Values
    !-------------
    real(kind=real_kind):: detDa

    detDa = Da(2,2)*Da(1,1) - Da(1,2)*Da(2,1)

    RR(1,1)=(Da(2,2)*Db(1,1) - Da(1,2)*Db(2,1))/detDa
    RR(1,2)=(Da(1,1)*Db(2,1) - Da(2,1)*Db(1,1))/detDa
    RR(2,1)=(Da(2,2)*Db(1,2) - Da(1,2)*Db(2,2))/detDa
    RR(2,2)=(Da(1,1)*Db(2,2) - Da(2,1)*Db(1,2))/detDa

    ! End Function
    !-------------
    return
  end function covariant_rot
  !==================================================================
#else
  !==================================================================
  function covariant_rot(Da,Db) result(RR)
    ! covariant_rot:  2 x 2 matrix multiply:  Db * Da^-1
    !                 for edge rotations: maps face a to face b
    ! ===========================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind):: Da(2,2)
    real(kind=real_kind):: Db(2,2)
    real(kind=real_kind):: RR(2,2)
    !
    ! Local Values
    !-------------
    real (kind=real_kind):: detDa

    detDa = Da(2,2)*Da(1,1) - Da(1,2)*Da(2,1)

    RR(1,1)=(Da(2,2)*Db(1,1) - Da(2,1)*Db(1,2))/detDa
    RR(1,2)=(Da(1,1)*Db(1,2) - Da(1,2)*Db(1,1))/detDa
    RR(2,1)=(Da(2,2)*Db(2,1) - Da(2,1)*Db(2,2))/detDa
    RR(2,2)=(Da(1,1)*Db(2,2) - Da(1,2)*Db(2,1))/detDa

    ! End Function
    !-------------
    return
  end function covariant_rot
  !==================================================================
#endif


  !==================================================================
  function contravariant_rot(Da,Db) result(RR)
    ! contravariant_rot: 2 x 2 matrix multiply:  Db^-1 * Da
    !                    that maps a contravariant vector field from an edge 
    !                    of cube face a to a contiguous edge of cube face b.
    !================================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind):: Da(2,2)
    real(kind=real_kind):: Db(2,2)
    real(kind=real_kind):: RR(2,2)
    !
    ! Local Values
    !-------------
    real(kind=real_kind):: detDb

    detDb = Db(2,2)*Db(1,1) - Db(1,2)*Db(2,1)

    RR(1,1)=(Da(1,1)*Db(2,2) - Da(2,1)*Db(1,2))/detDb
    RR(1,2)=(Da(1,2)*Db(2,2) - Da(2,2)*Db(1,2))/detDb
    RR(2,1)=(Da(2,1)*Db(1,1) - Da(1,1)*Db(2,1))/detDb
    RR(2,2)=(Da(2,2)*Db(1,1) - Da(1,2)*Db(2,1))/detDb

    ! End Function
    !-------------
    return
  end function contravariant_rot
  !==================================================================


  !==================================================================
  subroutine dmap(DD,elem,aa,bb)
    ! dmap: Initialize mapping that tranforms contravariant vector fields 
    !       on the reference element onto vector fields on the sphere. 
    !====================================================================
    use element_mod,only: element_t
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(out):: DD(2,2)
    type(element_t)                 :: elem
    real(kind=real_kind),intent(in ):: aa,bb

    if(cubed_sphere_map==0) then
      call dmap_equiangular(DD,elem,aa,bb)
    elseif(cubed_sphere_map==1) then
      call endrun('equi-distance gnomonic map not yet implemented')
    elseif(cubed_sphere_map==2) then
      call dmap_elementlocal(DD,elem,aa,bb)
    else
      call endrun('bad value of cubed_sphere_map')
    endif

    ! End Routine
    !-------------
    return
  end subroutine dmap
  !==================================================================


  !==================================================================
  subroutine dmap_equiangular(DD,elem,aa,bb)
    ! dmap_equiangular: Equiangular Gnomonic Projection Composition of equiangular 
    !                   Gnomonic projection to cubed-sphere face, followd by bilinear 
    !                   map to reference element
    !================================================================================
    use element_mod ,only: element_t
    use SE_Constants,only: np
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(out):: DD(2,2)
    type(element_t)                 :: elem
    real(kind=real_kind),intent(in ):: aa,bb
    !
    ! Local Values
    !-------------
    real(kind=real_kind):: tmpD(2,2),Jp(2,2),x1,x2,pi,pj,qi,qj
    real(kind=real_kind):: unif2quadmap(4,2)

#if 0
    ! we shoud get rid of elem%u2qmap() and routine cube_mod.F90::elem_jacobian()
    ! and replace with this code below:
    ! but this produces roundoff level changes
    !----------------------------------------------------------------
    unif2quadmap(1,1)=( elem%cartp( 1, 1)%x+elem%cartp(np, 1)%x &
                       +elem%cartp(np,np)%x+elem%cartp( 1,np)%x )/4.0d0
    unif2quadmap(1,2)=( elem%cartp( 1, 1)%y+elem%cartp(np, 1)%y &
                       +elem%cartp(np,np)%y+elem%cartp( 1,np)%y )/4.0d0
    unif2quadmap(2,1)=(-elem%cartp( 1, 1)%x+elem%cartp(np, 1)%x &
                       +elem%cartp(np,np)%x-elem%cartp( 1,np)%x )/4.0d0
    unif2quadmap(2,2)=(-elem%cartp( 1, 1)%y+elem%cartp(np, 1)%y &
                       +elem%cartp(np,np)%y-elem%cartp( 1,np)%y )/4.0d0
    unif2quadmap(3,1)=(-elem%cartp( 1, 1)%x-elem%cartp(np, 1)%x &
                       +elem%cartp(np,np)%x+elem%cartp( 1,np)%x )/4.0d0
    unif2quadmap(3,2)=(-elem%cartp( 1, 1)%y-elem%cartp(np, 1)%y &
                       +elem%cartp(np,np)%y+elem%cartp( 1,np)%y )/4.0d0
    unif2quadmap(4,1)=( elem%cartp( 1, 1)%x-elem%cartp(np, 1)%x &
                       +elem%cartp(np,np)%x-elem%cartp( 1,np)%x )/4.0d0
    unif2quadmap(4,2)=( elem%cartp( 1, 1)%y-elem%cartp(np, 1)%y &
                       +elem%cartp(np,np)%y-elem%cartp( 1,np)%y )/4.0d0
    Jp(1,1) = unif2quadmap(2,1) + unif2quadmap(4,1)*bb
    Jp(1,2) = unif2quadmap(3,1) + unif2quadmap(4,1)*aa
    Jp(2,1) = unif2quadmap(2,2) + unif2quadmap(4,2)*bb
    Jp(2,2) = unif2quadmap(3,2) + unif2quadmap(4,2)*aa
#else
    ! input (a,b) shold be a point in the reference element [-1,1]
    ! compute Jp(a,b)
    !----------------------------------------------------------------
    Jp(1,1) = elem%u2qmap(2,1) + elem%u2qmap(4,1)*bb
    Jp(1,2) = elem%u2qmap(3,1) + elem%u2qmap(4,1)*aa
    Jp(2,1) = elem%u2qmap(2,2) + elem%u2qmap(4,2)*bb
    Jp(2,2) = elem%u2qmap(3,2) + elem%u2qmap(4,2)*aa
#endif

    ! map (a,b) to the [-pi/2,pi/2] equi angular cube face:  x1,x2
    ! a = gp%points(i)
    ! b = gp%points(j)
    !----------------------------------------------------------------
    pi = (1-aa)/2
    pj = (1-bb)/2
    qi = (1+aa)/2
    qj = (1+bb)/2
    x1 = pi*pj*elem%corners(1)%x &
        +qi*pj*elem%corners(2)%x &
        +qi*qj*elem%corners(3)%x &
        +pi*qj*elem%corners(4)%x 
    x2 = pi*pj*elem%corners(1)%y &
        +qi*pj*elem%corners(2)%y &
        +qi*qj*elem%corners(3)%y &
        +pi*qj*elem%corners(4)%y 
    
    call vmap(tmpD,x1,x2,elem%vertex%face_number)

    ! Include map from element -> ref element in D
    !-----------------------------------------------
    DD(1,1) = tmpD(1,1)*Jp(1,1) + tmpD(1,2)*Jp(2,1)
    DD(1,2) = tmpD(1,1)*Jp(1,2) + tmpD(1,2)*Jp(2,2)
    DD(2,1) = tmpD(2,1)*Jp(1,1) + tmpD(2,2)*Jp(2,1)
    DD(2,2) = tmpD(2,1)*Jp(1,2) + tmpD(2,2)*Jp(2,2)

    ! End Routine
    !-------------
    return
  end subroutine dmap_equiangular
  !==================================================================


  !==================================================================
  subroutine vmap(DD,x1,x2,face_no) 
    ! vmap: Initialize mapping that tranforms contravariant vector fields on the 
    !       cube onto vector fields on the sphere. This follows Taylor's D matrix 
    !
    !       | cos(theta)dlambda/dx1  cos(theta)dlambda/dx2 |
    !   D = |                                              |
    !       |     dtheta/dx1              dtheta/dx2       |
    !
    !============================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(inout):: DD(2,2)
    real(kind=real_kind),intent(in   ):: x1
    real(kind=real_kind),intent(in   ):: x2
    integer             ,intent(in   ):: face_no
    !
    ! Local Values
    !-------------
    real (kind=real_kind) :: poledist  ! SQRT(TAN(x1)**2 +TAN(x2)**2)
    real (kind=real_kind) :: rr        ! distance from cube point to center of sphere
    real (kind=real_kind) :: D11
    real (kind=real_kind) :: D12
    real (kind=real_kind) :: D21
    real (kind=real_kind) :: D22

    rr=sqrt(1.0D0 + tan(x1)**2 + tan(x2)**2)

    if((face_no >= 1).and.(face_no <= 4)) then
      D11 = 1.0D0/(rr*cos(x1))
      D12 = 0.0D0
      D21 = -tan(x1)*tan(x2)/(cos(x1)*rr*rr)        
      D22 = 1.0D0/(rr*rr*cos(x1)*cos(x2)*cos(x2))
      DD(1,1) =  D11
      DD(1,2) =  D12
      DD(2,1) =  D21
      DD(2,2) =  D22
    elseif(face_no==6) then
      poledist=sqrt(tan(x1)**2 + tan(x2)**2)
      if(poledist <= DIST_THRESHOLD) then
        ! we set the D transform to the identity matrix 
        ! which works ONLY for swtc1, phi starting at 
        ! 3*PI/2... assumes lon at pole == 0
        !----------------------------------------------
        DD(1,1) =  1.0D0
        DD(1,2) =  0.0D0
        DD(2,1) =  0.0D0
        DD(2,2) =  1.0D0
      else
        D11 = -tan(x2)/(poledist*cos(x1)*cos(x1)*rr)
        D12 =  tan(x1)/(poledist*cos(x2)*cos(x2)*rr)
        D21 = -tan(x1)/(poledist*cos(x1)*cos(x1)*rr*rr)
        D22 = -tan(x2)/(poledist*cos(x2)*cos(x2)*rr*rr)
        DD(1,1) =  D11
        DD(1,2) =  D12
        DD(2,1) =  D21
        DD(2,2) =  D22
      endif
    elseif(face_no==5) then
      poledist=sqrt( tan(x1)**2 + tan(x2)**2)
      if(poledist <= DIST_THRESHOLD) then
        ! we set the D transform to the identity matrix 
        ! which works ONLY for swtc1, phi starting at 
        ! 3*PI/2... assumes lon at pole == 0, i.e. very specific
        !-------------------------------------------------------
        DD(1,1) =  1.0D0
        DD(1,2) =  0.0D0
        DD(2,1) =  0.0D0
        DD(2,2) =  1.0D0
      else
        D11 =  tan(x2)/(poledist*cos(x1)*cos(x1)*rr)
        D12 = -tan(x1)/(poledist*cos(x2)*cos(x2)*rr)
        D21 =  tan(x1)/(poledist*cos(x1)*cos(x1)*rr*rr)
        D22 =  tan(x2)/(poledist*cos(x2)*cos(x2)*rr*rr)
        DD(1,1) =  D11
        DD(1,2) =  D12
        DD(2,1) =  D21
        DD(2,2) =  D22
      endif
    endif

    ! End Routine
    !-------------
    return
  end subroutine vmap
  !==================================================================


  !==================================================================
  subroutine dmap_elementlocal(DD,elem,aa,bb)
    ! dmap_elementlocal: Initialize mapping that tranforms contravariant vector 
    !                    fields on the reference element onto vector fields on
    !                    the sphere. 
    !                    For Gnomonic, followed by bilinear, this code uses the old 
    !                    vmap() for unstructured grids, this code uses the parametric 
    !                    map that maps quads on the sphere directly to the reference 
    !                    element
    !==============================================================================
    use element_mod,only: element_t
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(out):: DD(2,2)
    type(element_t)                 :: elem
    real(kind=real_kind),intent(in ):: aa,bb
    !
    ! Local Values
    !------------------
    type(spherical_polar_t):: sphere
    type(cartesian3d_t)    :: corners(4)   
    real(kind=real_kind)   :: cc(3,4), qq(4), xx(3), rr, lam, th, dwrk(4,2)
    real(kind=real_kind)   :: sinlam, sinth, coslam, costh
    real(kind=real_kind)   :: D1(2,3), D2(3,3), D3(3,2), D4(3,2)
    integer :: ii,jj

    sphere = ref2sphere(aa,bb,elem)
    corners = elem%corners3D

    cc(1,1)=corners(1)%x
    cc(2,1)=corners(1)%y
    cc(3,1)=corners(1)%z
    cc(1,2)=corners(2)%x
    cc(2,2)=corners(2)%y
    cc(3,2)=corners(2)%z
    cc(1,3)=corners(3)%x
    cc(2,3)=corners(3)%y
    cc(3,3)=corners(3)%z
    cc(1,4)=corners(4)%x
    cc(2,4)=corners(4)%y
    cc(3,4)=corners(4)%z

    qq(1)=(1-aa)*(1-bb)
    qq(2)=(1+aa)*(1-bb)
    qq(3)=(1+aa)*(1+bb)
    qq(4)=(1-aa)*(1+bb)
    qq=qq/4.0d0

    do ii=1,3
      xx(ii)=sum(cc(ii,:)*qq(:))
    end do

    rr=sqrt(xx(1)**2 + xx(2)**2 + xx(3)**2)

    lam   =sphere%lon
    th    =sphere%lat
    sinlam=sin(lam)
    sinth =sin( th)
    coslam=cos(lam)
    costh =cos( th)

    D1(1,1)=-sinlam
    D1(1,2)= coslam  
    D1(1,3)=0.0d0 
    D1(2,1)=0.0d0
    D1(2,2)=0.0d0
    D1(2,3)=1.0d0

    D2(1,1)=(sinlam**2)*(costh**2)+(sinth**2)
    D2(1,2)=-sinlam*coslam*(costh**2)
    D2(1,3)=-coslam*sinth*costh
    D2(2,1)=-sinlam*coslam*(costh**2)
    D2(2,2)=(coslam**2)*(costh**2)+sinth**2
    D2(2,3)=-sinlam*sinth*costh
    D2(3,1)=-coslam*sinth           
    D2(3,2)=-sinlam*sinth
    D2(3,3)= costh

    dwrk(1,1)=-1+bb
    dwrk(1,2)=-1+aa
    dwrk(2,1)= 1-bb
    dwrk(2,2)=-1-aa
    dwrk(3,1)= 1+bb
    dwrk(3,2)= 1+aa
    dwrk(4,1)=-1-bb
    dwrk(4,2)= 1-aa

    dwrk=dwrk/4.0d0

    do ii=1,3
    do jj=1,2
      D3(ii,jj)=sum(cc(ii,:)*dwrk(:,jj))
    end do
    end do

    do ii=1,3
    do jj=1,2
      D4(ii,jj)=sum(D2(ii,:)*D3(:,jj))
    end do
    end do   

    do ii=1,2
    do jj=1,2
      DD(ii,jj)=sum(D1(ii,:)*D4(:,jj))
    end do
    end do

    DD=DD/rr

    ! End Routine
    !-------------
    return
  end subroutine dmap_elementlocal
  !==================================================================


  !==================================================================
  subroutine coreolis_init_atomic(elem)
    ! coreolis_init_atomic: Initialize coreolis term ...
    !=======================================================
    use element_mod ,only: element_t
    use SE_Constants,only: np
    use SE_Constants,only: omega
    !
    ! Passed Variables
    !------------------
    type(element_t):: elem
    !
    ! Local Values
    !------------------
    real(kind=real_kind):: lat,lon,rangle
    integer             :: ii,jj

    rangle = rotate_grid*DD_PI/180.d0
    do jj=1,np
    do ii=1,np
      if(rotate_grid /= 0) then
        lat = elem%spherep(ii,jj)%lat
        lon = elem%spherep(ii,jj)%lon
        elem%fcor(ii,jj)= 2.d0*omega*(-cos(lon)*cos(lat)*sin(rangle)+sin(lat)*cos(rangle))
      else
        elem%fcor(ii,jj)= 2.0d0*omega*sin(elem%spherep(ii,jj)%lat)
      endif
    end do
    end do

    ! End Routine
    !-------------
    return
  end subroutine coreolis_init_atomic
  !==================================================================


  !==================================================================
  subroutine rotation_init_atomic(elem, rot_type)
    ! rotation_init_atomic: Initialize cube rotation terms resulting
    !                       from changing cube face coordinate systems
    !==================================================================
    use element_mod ,only: element_t
    use SE_Constants,only: np
    use SE_Options  ,only: north, south, east, west, neast, seast, swest, nwest
    !
    ! Passed Variables
    !------------------
    type(element_t) :: elem
    character(len=*)   rot_type
    !
    ! Local Values
    !------------------
    integer             :: myface_no        ! current element face number
    integer             :: nbrface_no       ! neighbor element face number
    integer             :: inbr
    integer             :: nrot,irot
    integer             :: ii,jj,kk
    integer             :: ir,jr
    integer             :: start, cnt
    real(kind=real_kind):: Dloc(2,2,np)
    real(kind=real_kind):: Drem(2,2,np)
    real(kind=real_kind):: x1,x2

    myface_no = elem%vertex%face_number
    nrot      = 0

    do inbr=1,8
      cnt  = elem%vertex%nbrs_ptr(inbr+1) - elem%vertex%nbrs_ptr(inbr) 
      start= elem%vertex%nbrs_ptr(inbr) 
      do kk=0,(cnt-1)
        nbrface_no = elem%vertex%nbrs_face(start+kk)
        if(myface_no /= nbrface_no) nrot=nrot+1
      end do
    end do

    if(associated(elem%desc%rot)) then
      if(size(elem%desc%rot) > 0) then
        !         deallocate(elem%desc%rot)
        nullify(elem%desc%rot)
      endif
    endif

    ! =====================================================
    ! If there are neighbors on other cube faces, allocate 
    ! an array of rotation matrix structs.
    ! =====================================================
    if(nrot > 0) then
      allocate(elem%desc%rot(nrot))
      elem%desc%use_rotation=1
      irot=0          

      do inbr=1,8
        cnt  = elem%vertex%nbrs_ptr(inbr+1) - elem%vertex%nbrs_ptr(inbr) 
        start= elem%vertex%nbrs_ptr(inbr) 

        do kk= 0, cnt-1
          nbrface_no = elem%vertex%nbrs_face(start+kk)

          ! The cube edge (myface_no,nbrface_no) and inbr defines 
          ! a unique rotation given by (D^-1) on myface_no x (D) on nbrface_no
          !--------------------------------------------------------------------
          if((myface_no /= nbrface_no).and.(elem%vertex%nbrs(start+kk) /= -1)) then           
            irot=irot+1

            if(inbr <= 4) then      
              allocate(elem%desc%rot(irot)%R(2,2,np))  ! edge
            else                     
              allocate(elem%desc%rot(irot)%R(2,2,1 ))  ! corner
            endif

            ! must compute Dloc on my face, Drem on neighbor face, 
            ! for each point on edge or corner.
                
            ! Equatorial belt east/west neighbors
            !=========================================
            if((nbrface_no <= 4).and.(myface_no <= 4)) then
              if(inbr == west) then
                do jj=1,np
                  x1 = elem%cartp( 1,jj)%x
                  x2 = elem%cartp( 1,jj)%y
                  call Vmap(Dloc(1,1,jj), x1,x2, myface_no)
                  call Vmap(Drem(1,1,jj),-x1,x2,nbrface_no)
                end do
              elseif(inbr == east) then
                do jj=1,np
                  x1 = elem%cartp(np,jj)%x
                  x2 = elem%cartp(np,jj)%y
                  call Vmap(Dloc(1,1,jj), x1,x2, myface_no)
                  call Vmap(Drem(1,1,jj),-x1,x2,nbrface_no)
                end do
              elseif(inbr == swest) then
                x1 = elem%cartp(1,1)%x
                x2 = elem%cartp(1,1)%y
                call Vmap(Dloc(1,1,1), x1,x2, myface_no)
                call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
              elseif(inbr == nwest) then
                x1 = elem%cartp(1,np)%x
                x2 = elem%cartp(1,np)%y
                call Vmap(Dloc(1,1,1), x1,x2, myface_no)
                call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
              elseif(inbr == seast) then
                x1 = elem%cartp(np,1)%x
                x2 = elem%cartp(np,1)%y
                call Vmap(Dloc(1,1,1), x1,x2, myface_no)
                call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
              elseif(inbr == neast) then
                x1 = elem%cartp(np,np)%x
                x2 = elem%cartp(np,np)%y
                call Vmap(Dloc(1,1,1), x1,x2, myface_no)
                call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
              endif
            endif
                
            ! Northern Neighbors of Equatorial Belt
            !=========================================
            if((myface_no <= 4).and.(nbrface_no == 6)) then
              if(inbr == north) then
                do ii=1,np
                  ir=np+1-ii
                  x1=elem%cartp(ii,np)%x
                  x2=elem%cartp(ii,np)%y
                  if(myface_no == 1) then
                    call Vmap(Dloc(1,1,ii), x1, x2, myface_no)
                    call Vmap(Drem(1,1,ii), x1,-x2,nbrface_no)
                  endif
                  if(myface_no == 2) then
                    call Vmap(Dloc(1,1,ii), x1,x2, myface_no)
                    call Vmap(Drem(1,1,ii), x2,x1,nbrface_no)
                  endif
                  if(myface_no == 3) then
                    call Vmap(Dloc(1,1,ir), x1,x2, myface_no)
                    call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                  endif
                  if(myface_no == 4) then
                    call Vmap(Dloc(1,1,ir), x1,x2,myface_no)
                    call Vmap(Drem(1,1,ir),-x2,-x1,nbrface_no)
                  endif
                end do
              elseif(inbr == nwest) then
                x1 = elem%cartp(1,np)%x
                x2 = elem%cartp(1,np)%y
                call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                if( myface_no == 1) call Vmap(Drem(1,1,1), x1,-x2,nbrface_no)
                if( myface_no == 2) call Vmap(Drem(1,1,1), x2, x1,nbrface_no)
                if( myface_no == 3) call Vmap(Drem(1,1,1),-x1, x2,nbrface_no)
                if( myface_no == 4) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
              elseif(inbr == neast) then
                x1 = elem%cartp(np,np)%x
                x2 = elem%cartp(np,np)%y
                call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                if( myface_no == 1) call Vmap(Drem(1,1,1), x1,-x2,nbrface_no)
                if( myface_no == 2) call Vmap(Drem(1,1,1), x2, x1,nbrface_no)
                if( myface_no == 3) call Vmap(Drem(1,1,1),-x1, x2,nbrface_no)
                if( myface_no == 4) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
              endif
            endif
                
            ! Southern Neighbors of Equatorial Belt
            !=========================================
            if((myface_no <= 4).and.(nbrface_no == 5)) then
              if(inbr == south) then
                do ii=1,np
                  ir=np+1-ii
                  x1 = elem%cartp(ii,1)%x
                  x2 = elem%cartp(ii,1)%y
                  if(myface_no == 1) then
                    call Vmap(Dloc(1,1,ii), x1, x2, myface_no)
                    call Vmap(Drem(1,1,ii), x1,-x2,nbrface_no)
                  endif
                  if(myface_no == 2) then
                    call Vmap(Dloc(1,1,ir), x1, x2, myface_no)
                    call Vmap(Drem(1,1,ir),-x2,-x1,nbrface_no)
                  endif
                  if(myface_no == 3) then
                    call Vmap(Dloc(1,1,ir), x1,x2, myface_no)
                    call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                  endif
                  if(myface_no == 4) then
                    call Vmap(Dloc(1,1,ii), x1,x2, myface_no)
                    call Vmap(Drem(1,1,ii), x2,x1,nbrface_no)
                  endif
                end do
              elseif(inbr == swest) then
                x1 = elem%cartp(1,1)%x
                x2 = elem%cartp(1,1)%y
                call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                if( myface_no == 1) call Vmap(Drem(1,1,1), x1,-x2,nbrface_no)
                if( myface_no == 2) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                if( myface_no == 3) call Vmap(Drem(1,1,1),-x1, x2,nbrface_no)
                if( myface_no == 4) call Vmap(Drem(1,1,1), x2, x1,nbrface_no)
              elseif(inbr == seast) then
                x1 = elem%cartp(np,1)%x
                x2 = elem%cartp(np,1)%y
                call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                if( myface_no == 1) call Vmap(Drem(1,1,1), x1,-x2,nbrface_no)
                if( myface_no == 2) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                if( myface_no == 3) call Vmap(Drem(1,1,1),-x1, x2,nbrface_no)
                if( myface_no == 4) call Vmap(Drem(1,1,1), x2, x1,nbrface_no)
              endif
            endif
                
            ! Neighbors of Northern Capping Face Number 6
            !=========================================
            if(myface_no == 6) then
              if(nbrface_no == 1) then
                if(inbr == south) then
                  do ii=1,np
                    x1 = elem%cartp(ii,1)%x
                    x2 = elem%cartp(ii,1)%y
                    call Vmap(Dloc(1,1,ii),x1, x2, myface_no)
                    call Vmap(Drem(1,1,ii),x1,-x2,nbrface_no)
                  end do
                elseif(inbr == swest) then
                  x1 = elem%cartp(1,1)%x
                  x2 = elem%cartp(1,1)%y
                  call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                  call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                elseif(inbr == seast) then
                  x1 = elem%cartp(np,1)%x
                  x2 = elem%cartp(np,1)%y
                  call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                  call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                endif
              elseif(nbrface_no == 2) then
                if(inbr == east) then
                  do jj=1,np
                    x1 = elem%cartp(np,jj)%x
                    x2 = elem%cartp(np,jj)%y
                    call Vmap(Dloc(1,1,jj),x1,x2, myface_no)
                    call Vmap(Drem(1,1,jj),x2,x1,nbrface_no)
                  end do
                elseif(inbr == seast) then
                  x1 = elem%cartp(np,1)%x
                  x2 = elem%cartp(np,1)%y
                  call Vmap(Dloc(1,1,1),x1,x2, myface_no)
                  call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                elseif(inbr == neast) then
                  x1 = elem%cartp(np,np)%x
                  x2 = elem%cartp(np,np)%y
                  call Vmap(Dloc(1,1,1),x1,x2, myface_no)
                  call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                endif
              elseif(nbrface_no == 3) then
                if(inbr == north) then
                  do ii=1,np
                    ir =np+1-ii
                    x1 = elem%cartp(ii,np)%x
                    x2 = elem%cartp(ii,np)%y
                    call Vmap(Dloc(1,1,ir), x1,x2, myface_no)
                    call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                  end do
                elseif(inbr == nwest) then
                  x1 = elem%cartp(1,np)%x
                  x2 = elem%cartp(1,np)%y
                  call Vmap(Dloc(1,1,1), x1,x2, myface_no)
                  call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                elseif(inbr == neast) then
                  x1 = elem%cartp(np,np)%x
                  x2 = elem%cartp(np,np)%y
                  call Vmap(Dloc(1,1,1), x1,x2, myface_no)
                  call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                endif
              elseif(nbrface_no == 4) then
                if(inbr == west) then
                  do jj=1,np
                    jr=np+1-jj
                    x1 = elem%cartp(1,jj)%x
                    x2 = elem%cartp(1,jj)%y
                    call Vmap(Dloc(1,1,jr), x1, x2, myface_no)
                    call Vmap(Drem(1,1,jr),-x2,-x1,nbrface_no)
                  end do
                elseif(inbr == swest) then
                  x1 = elem%cartp(1,1)%x
                  x2 = elem%cartp(1,1)%y
                  call Vmap(Dloc(1,1,1), x1, x2, myface_no)
                  call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                elseif(inbr == nwest) then
                  x1 = elem%cartp(1,np)%x
                  x2 = elem%cartp(1,np)%y
                  call Vmap(Dloc(1,1,1), x1, x2, myface_no)
                  call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                endif
              endif
            endif
                
            ! Neighbors of South Capping Face Number 5
            !=========================================
            if(myface_no == 5) then
              if(nbrface_no == 1) then
                if(inbr == north) then
                  do ii=1,np
                    x1 = elem%cartp(ii,np)%x
                    x2 = elem%cartp(ii,np)%y
                    call Vmap(Dloc(1,1,ii),x1, x2, myface_no)
                    call Vmap(Drem(1,1,ii),x1,-x2,nbrface_no)
                  end do
                elseif(inbr == nwest) then
                  x1 = elem%cartp(1,np)%x
                  x2 = elem%cartp(1,np)%y
                  call Vmap(Dloc(:,:,1),x1, x2, myface_no)
                  call Vmap(Drem(:,:,1),x1,-x2,nbrface_no)
                elseif(inbr == neast) then
                  x1 = elem%cartp(np,np)%x
                  x2 = elem%cartp(np,np)%y
                  call Vmap(Dloc(1,1,1),x1, x2, myface_no)
                  call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                endif
              elseif(nbrface_no == 2) then
                if(inbr == east) then
                  do jj=1,np
                    jr=np+1-jj
                    x1 = elem%cartp(np,jj)%x
                    x2 = elem%cartp(np,jj)%y
                    call Vmap(Dloc(1,1,jr), x1, x2, myface_no)
                    call Vmap(Drem(1,1,jr),-x2,-x1,nbrface_no)
                  end do
                elseif(inbr == seast) then
                  x1 = elem%cartp(np,1)%x
                  x2 = elem%cartp(np,1)%y
                  call Vmap(Dloc(1,1,1), x1, x2, myface_no)
                  call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                elseif(inbr == neast) then
                  x1 = elem%cartp(np,np)%x
                  x2 = elem%cartp(np,np)%y
                  call Vmap(Dloc(1,1,1), x1, x2, myface_no)
                  call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                endif
              elseif(nbrface_no == 3) then
                if(inbr == south) then
                  do ii=1,np
                    ir=np+1-ii
                    x1 = elem%cartp(ii,1)%x
                    x2 = elem%cartp(ii,1)%y
                    call Vmap(Dloc(1,1,ir), x1,x2, myface_no)
                    call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                  end do
                elseif(inbr == swest) then
                  x1 = elem%cartp(1,1)%x
                  x2 = elem%cartp(1,1)%y
                  call Vmap(Dloc(1,1,1), x1,x2, myface_no)
                  call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                elseif(inbr == seast) then
                  x1 = elem%cartp(np,1)%x
                  x2 = elem%cartp(np,1)%y
                  call Vmap(Dloc(1,1,1), x1,x2, myface_no)
                  call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                endif
              elseif(nbrface_no == 4) then
                if(inbr == west) then
                  do jj=1,np
                    x1 = elem%cartp(1,jj)%x
                    x2 = elem%cartp(1,jj)%y
                    call Vmap(Dloc(1,1,jj),x1,x2, myface_no)
                    call Vmap(Drem(1,1,jj),x2,x1,nbrface_no)
                  end do
                elseif(inbr == swest) then
                  x1 = elem%cartp(1,1)%x
                  x2 = elem%cartp(1,1)%y
                  call Vmap(Dloc(1,1,1),x1,x2, myface_no)
                  call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                elseif(inbr == nwest) then
                  x1 = elem%cartp(1,np)%x
                  x2 = elem%cartp(1,np)%y
                  call Vmap(Dloc(1,1,1),x1,x2, myface_no)
                  call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                endif
              endif
            endif
                
            elem%desc%rot(irot)%nbr = inbr
            if(rot_type == "covariant") then
              do ii=1,size(elem%desc%rot(irot)%R(:,:,:),3)
                elem%desc%rot(irot)%R(:,:,ii)=covariant_rot(Dloc(:,:,ii),Drem(:,:,ii))
              end do
            elseif(rot_type == "contravariant") then
              do ii=1,size(elem%desc%rot(irot)%R(:,:,:),3)
                elem%desc%rot(irot)%R(:,:,ii)=contravariant_rot(Dloc(:,:,ii),Drem(:,:,ii))
              end do
            endif
                
          endif ! end of a unique rotation
        end do !k loop over neighbors in that direction
      end do !inbr loop
    endif !nrot > 0
    
    ! End Routine
    !-------------
    return
  end subroutine rotation_init_atomic
  !==================================================================
  

  !==================================================================
  subroutine set_corner_coordinates(elem)
    ! set_corner_coordinates:
    !
    !==========================================================
    use element_mod,only: element_t 
    !
    ! Passed Variables
    !-------------------
    type (element_t) :: elem 
    !
    ! Local Values
    !---------------
    integer              :: ii,ie,je,face_no,nn
    real (kind=real_kind):: dx,dy,startx,starty

    if(ne==0) call endrun('Error in set_corner_coordinates: ne is zero')

    ! compute cube face coordinates of element
    !------------------------------------------
    call convert_gbl_index(elem%vertex%number,ie,je,face_no)
    elem%vertex%face_number = face_no 

    dx = (cube_xend-cube_xstart)/ne
    dy = (cube_yend-cube_ystart)/ne

    startx = cube_xstart + ie*dx
    starty = cube_ystart + je*dy

    elem%corners(1)%x = startx
    elem%corners(1)%y = starty
    elem%corners(2)%x = startx+dx
    elem%corners(2)%y = starty
    elem%corners(3)%x = startx+dx
    elem%corners(3)%y = starty+dy
    elem%corners(4)%x = startx   
    elem%corners(4)%y = starty+dy

    do ii=1,4
      elem%node_multiplicity(ii) = 4
    end do  
    ie = ie + 1
    je = je + 1
    if((ie == 1).and.(je == 1)) then 
      elem%node_multiplicity(1) = 3
    elseif((ie == ne).and.(je ==  1)) then 
      elem%node_multiplicity(2) = 3
    elseif((ie == ne).and.(je == ne)) then
      elem%node_multiplicity(3) = 3
    elseif((ie ==  1).and.(je == ne)) then
      elem%node_multiplicity(4) = 3
    endif  

    ! End Routine
    !-------------
    return
  end subroutine set_corner_coordinates
  !==================================================================


  !==================================================================
  subroutine assign_node_numbers_to_elem(elements, GridVertex)
    ! assign_node_numbers_to_elem:
    !
    !============================================================
    use SE_Options   ,only: north, south, east, west, neast, seast, swest, nwest
    use element_mod  ,only: element_t
    use gridgraph_mod,only: GridVertex_t
    !
    ! Passed Variables
    !--------------------
    type(element_t)   ,intent(inout):: elements(:)
    type(GridVertex_t),intent(in   ):: GridVertex(:)
    !
    ! Local Values
    !----------------
    type (GridVertex_t):: vertex
    integer            :: connectivity(6*ne*ne, 4)
    integer            :: nn(4), en(4)
    integer            :: el, ii, jj, direction
    integer            :: current_node_num, tot_ne
    integer            :: start, cnt


    current_node_num = 0
    tot_ne           = 6*ne*ne

    if(ne==0) call endrun('Error in assign_node_numbers_to_elem: ne is zero')
    if(tot_ne /= size(GridVertex)) then
      call endrun('Error in assign_node_numbers_to_elem: GridVertex not correct length')
    endif

    connectivity(:,:) = 0

    do el = 1,tot_ne  
      vertex= GridVertex(el)
      en    = 0 
      do direction = 1,8
        cnt  = vertex%nbrs_ptr(direction+1) - vertex%nbrs_ptr(direction) 
        start= vertex%nbrs_ptr(direction) 

        do ii=0, cnt-1
          jj = vertex%nbrs(start+ii)
          if(jj /= -1) then
            nn(:) = connectivity(jj,:)
            select case (direction)
              case (north)
                   if (nn(1)/=0) en(4) = nn(1)
                   if (nn(2)/=0) en(3) = nn(2)
              case (south)
                   if (nn(4)/=0) en(1) = nn(4)
                   if (nn(3)/=0) en(2) = nn(3)
              case (east)
                   if (nn(1)/=0) en(2) = nn(1)
                   if (nn(4)/=0) en(3) = nn(4)
              case (west)
                   if (nn(2)/=0) en(1) = nn(2)
                   if (nn(3)/=0) en(4) = nn(3)
              case (neast)
                   if (nn(1)/=0) en(3) = nn(1)
              case (seast)
                   if (nn(4)/=0) en(2) = nn(4)
              case (swest)
                   if (nn(3)/=0) en(1) = nn(3)
              case (nwest)
                   if (nn(2)/=0) en(4) = nn(2)
            end select
          endif
        end do
      end do !direction

      do ii=1,4
        if(en(ii) == 0) then
          current_node_num = current_node_num + 1
          en(ii) = current_node_num
        end if
      end do
      connectivity(el,:) = en
    end do ! el = 1,tot_ne  

    if(current_node_num /= (6*ne*ne+2)) then
      call endrun('Error in assignment of node numbers: Failed Euler test')
    endif
    do el = 1,size(elements)
      elements(el)%node_numbers = connectivity(elements(el)%vertex%number, :)
    end do

    ! End Routine
    !-------------
    return
  end subroutine assign_node_numbers_to_elem
  !==================================================================


  !==================================================================
  subroutine convert_gbl_index(number,ie,je,face_no)
    ! convert_gbl_index: Convert global element index to cube index
    !==============================================================
    !
    ! Passed Variables
    !--------------------
    integer,intent(in ):: number
    integer,intent(out):: ie,je,face_no

    if(ne==0) call endrun('Error in cube_mod:convert_gbl_index: ne is zero')

    !  inverse of the function:      number = 1 + ie + ne*je + ne*ne*(face_no-1)
    !--------------------------------------------------------------------------
    face_no=((number-1)/(ne*ne))+1
    ie     =modulo(number-1,ne)
    je     =(number-1)/ne - (face_no-1)*ne

    ! End Routine
    !-------------
    return
  end subroutine convert_gbl_index
  !==================================================================
   

  !==================================================================
  subroutine CubeTopology(GridEdge, GridVertex)
    ! CubeTopology:
    !
    ! Since GridVertex fields must be allocated before calling this, it
    ! must be intent(inout).
    !og: is 'target' here necessary?
    !GridEdge : changed its 'out' attribute to 'inout'
    !==============================================================
    use SE_Constants  ,only: np
    use SE_Options    ,only: north,south,east,west,neast,seast,swest,nwest
    use SE_Options    ,only: RECURSIVE, SFCURVE
    use gridgraph_mod ,only: GridEdge_t, GridVertex_t, initgridedge, PrintGridEdge, &
                             allocate_gridvertex_nbrs, deallocate_gridvertex_nbrs 
    use spacecurve_mod,only: IsFactorable, genspacecurve
    !
    ! Passed Variables
    !-------------------
    type(GridEdge_t)  ,intent(inout),target:: GridEdge(:)
    type(GridVertex_t),intent(inout),target:: GridVertex(:)
    !
    ! Local Values
    !------------------
    integer           ,allocatable:: Mesh(:,:)
    integer           ,allocatable:: Mesh2(:,:),Mesh2_map(:,:,:),sfcij(:,:)
    type(GridVertex_t),allocatable:: GridElem(:,:,:)
    logical           ,allocatable:: nbrs_used(:,:,:,:)

    integer:: ii,jj,kk,ll,number,irev,ne2,i2,j2,sfc_index
    integer:: EdgeWgtP,CornerWgt
    integer:: ielem, nedge
    integer:: offset, ierr, loc
    

    if(ne==0) call endrun('Error in CubeTopology: ne is zero')

    ! Allocate and init arrays
    !------------------------
    allocate(GridElem(ne,ne,nfaces),stat=ierr)
    do kk = 1, nfaces
    do jj = 1, ne
    do ii = 1, ne
      call allocate_gridvertex_nbrs(GridElem(ii,jj,kk))
    end do
    end do
    end do
    if(ierr/=0) then
       call endrun('error in allocation of GridElem structure')
    end if

    allocate(nbrs_used(ne,ne,nfaces,8))
    nbrs_used(:,:,:,:) = .false.

    number    = 1
    EdgeWgtP  = np
    CornerWgt = 1
    do kk=1,nfaces
    do jj=1,ne
    do ii=1,ne
      ! Number elements
      !---------------------
      GridElem(ii,jj,kk)%nbrs          (:)=0
      GridElem(ii,jj,kk)%nbrs_wgt      (:)=0
      GridElem(ii,jj,kk)%nbrs_ptr      (:)=0
      GridElem(ii,jj,kk)%nbrs_wgt_ghost(:)=1  ! always this value
      GridElem(ii,jj,kk)%SpaceCurve       =0
      GridElem(ii,jj,kk)%number           =number 
      number=number+1
    end do
    end do
    end do

!   print *,'CubeTopology: Ne,IsFactorable,IsLoadBalanced:',ne,IsFactorable(ne),IsLoadBalanced(nelem,npart)

    allocate(Mesh(ne,ne))
    if(IsFactorable(ne)) then
      call GenspaceCurve(Mesh)
!      call PrintCurve(Mesh) 
    else
      ! find the smallest ne2 which is a power of 2 and ne2>ne
      !-------------------------------------------------------
      ne2=2**ceiling( log(real(ne))/log(2d0) )
      if(ne2<ne) call endrun('Fatel SFC error')

      allocate(Mesh2    (ne2,ne2  ))
      allocate(Mesh2_map(ne2,ne2,2))
      allocate(sfcij  (0:ne2*ne2,2))

      call GenspaceCurve(Mesh2)  ! SFC partition for ne2

      ! associate every element on the neXne mesh (Mesh)
      ! with its closest element on the ne2Xne2 mesh (Mesh2)
      ! Store this as a map from Mesh2 -> Mesh in Mesh2_map.
      ! elements in Mesh2 which are not mapped get assigned a value of 0
      !----------------------------------------------------------------------
      Mesh2_map=0
      do jj=1,ne
      do ii=1,ne
        ! map this element to an (i2,j2) element
        ! [ (i-.5)/ne , (j-.5)/ne ]  = [ (i2-.5)/ne2 , (j2-.5)/ne2 ]
        !------------------------------------------------------------
        i2=nint(((ii-.5)/ne)*ne2 + .5 )
        j2=nint(((jj-.5)/ne)*ne2 + .5 )
        if(i2<1  ) i2=1
        if(i2>ne2) i2=ne2
        if(j2<1  ) j2=1
        if(j2>ne2) j2=ne2
        Mesh2_map(i2,j2,1)=ii
        Mesh2_map(i2,j2,2)=jj
      end do
      end do

      ! create a reverse index array for Mesh2
      ! k = Mesh2(i,j) 
      ! (i,j) = (sfcij(k,1),sfci(k,2)) 
      !----------------------------------------
      do jj=1,ne2
      do ii=1,ne2
        kk=Mesh2(ii,jj)
        sfcij(kk,1)=ii
        sfcij(kk,2)=jj
      end do
      end do

      ! generate a SFC for Mesh with the same ordering as the 
      ! elements in Mesh2 which map to Mesh.
      !-------------------------------------------------------
      sfc_index=0
      do kk=0,ne2*ne2-1
        i2=sfcij(kk,1)
        j2=sfcij(kk,2)
        ii=Mesh2_map(i2,j2,1)
        jj=Mesh2_map(i2,j2,2)
        if(ii/=0) then
          ! (i2,j2) element maps to (i,j) element
          !---------------------------------------
          Mesh(ii,jj)=sfc_index
          sfc_index=sfc_index+1
        endif
      end do
#if 0
      print *,'SFC Mapping to non powers of 2,3 used.  Mesh:'  
      do jj=1,ne
        write(*,'(99i3)') (Mesh(ii,jj),ii=1,ne)
      end do
      call PrintCurve(Mesh2) 
#endif
      deallocate(Mesh2    )
      deallocate(Mesh2_map)
      deallocate(sfcij    )
    endif

    !  Setup the space-filling curve for face 1
    ! -------------------------------------------
    offset=0
    do jj=1,ne
    do ii=1,ne
      GridElem(ii,jj,1)%SpaceCurve = offset + Mesh(ii,ne-jj+1)
    end do
    end do

    !  Setup the space-filling curve for face 2
    ! -------------------------------------------
    offset = offset + ne*ne
    do jj=1,ne
    do ii=1,ne
      GridElem(ii,jj,2)%SpaceCurve = offset + Mesh(ii,ne-jj+1)
    end do
    end do

    !  Setup the space-filling curve for face 6
    ! -------------------------------------------
    offset = offset + ne*ne
    do jj=1,ne
    do ii=1,ne
      GridElem(ii,jj,6)%SpaceCurve = offset + Mesh(ne-ii+1,ne-jj+1)
    end do
    end do

    !  Setup the space-filling curve for face 4
    ! -------------------------------------------
    offset = offset + ne*ne
    do jj=1,ne
    do ii=1,ne
      GridElem(ii,jj,4)%SpaceCurve = offset + Mesh(ne-jj+1,ii)
    end do
    end do

    !  Setup the space-filling curve for face 5
    ! -------------------------------------------
    offset = offset + ne*ne
    do jj=1,ne
    do ii=1,ne
      GridElem(ii,jj,5)%SpaceCurve = offset + Mesh(ii,jj)
    end do
    end do

    !  Setup the space-filling curve for face 3
    ! -------------------------------------------
    offset = offset + ne*ne
    do jj=1,ne
    do ii=1,ne
      GridElem(ii,jj,3)%SpaceCurve = offset + Mesh(ii,jj)
    end do
    end do

    !==================
    ! face interiors
    !==================
    do kk=1,6
      ! setup  SOUTH, WEST, SW neighbors
      !---------------------------------
      do jj=2,ne
      do ii=2,ne
        nbrs_used(ii,jj,kk, west) = .true.
        nbrs_used(ii,jj,kk,south) = .true.
        nbrs_used(ii,jj,kk,swest) = .true.

        GridElem(ii,jj,kk)%nbrs     ( west) = GridElem(ii-1,jj  ,kk)%number
        GridElem(ii,jj,kk)%nbrs_face( west) = kk
        GridElem(ii,jj,kk)%nbrs_wgt ( west) = EdgeWgtP
        GridElem(ii,jj,kk)%nbrs     (south) = GridElem(ii  ,jj-1,kk)%number
        GridElem(ii,jj,kk)%nbrs_face(south) = kk
        GridElem(ii,jj,kk)%nbrs_wgt (south) = EdgeWgtP
        GridElem(ii,jj,kk)%nbrs     (swest) = GridElem(ii-1,jj-1,kk)%number
        GridElem(ii,jj,kk)%nbrs_face(swest) = kk
        GridElem(ii,jj,kk)%nbrs_wgt (swest) = CornerWgt
       end do
       end do

       !  setup EAST, NORTH, NE neighbors
       !---------------------------------
       do jj=1,ne-1
       do ii=1,ne-1
         nbrs_used(ii,jj,kk, east) = .true.
         nbrs_used(ii,jj,kk,north) = .true.
         nbrs_used(ii,jj,kk,neast) = .true.
             
         GridElem(ii,jj,kk)%nbrs     ( east) = GridElem(ii+1,jj  ,kk)%number
         GridElem(ii,jj,kk)%nbrs_face( east) = kk
         GridElem(ii,jj,kk)%nbrs_wgt ( east) = EdgeWgtP
         GridElem(ii,jj,kk)%nbrs     (north) = GridElem(ii  ,jj+1,kk)%number
         GridElem(ii,jj,kk)%nbrs_face(north) = kk
         GridElem(ii,jj,kk)%nbrs_wgt (north) = EdgeWgtP
         GridElem(ii,jj,kk)%nbrs     (neast) = GridElem(ii+1,jj+1,kk)%number
         GridElem(ii,jj,kk)%nbrs_face(neast) = kk
         GridElem(ii,jj,kk)%nbrs_wgt (neast) = CornerWgt
       end do
       end do

       ! Setup the remaining SOUTH, EAST, and SE neighbors
       !---------------------------------------------------
       do jj=2,ne
       do ii=1,ne-1
         nbrs_used(ii,jj,kk,south) = .true.
         nbrs_used(ii,jj,kk, east) = .true.
         nbrs_used(ii,jj,kk,seast) = .true.
             
         GridElem(ii,jj,kk)%nbrs     (south) = GridElem(ii  ,jj-1,kk)%number
         GridElem(ii,jj,kk)%nbrs_face(south) = kk
         GridElem(ii,jj,kk)%nbrs_wgt (south) = EdgeWgtP
         GridElem(ii,jj,kk)%nbrs     ( east) = GridElem(ii+1,jj  ,kk)%number
         GridElem(ii,jj,kk)%nbrs_face( east) = kk
         GridElem(ii,jj,kk)%nbrs_wgt ( east) = EdgeWgtP
         GridElem(ii,jj,kk)%nbrs     (seast) = GridElem(ii+1,jj-1,kk)%number
         GridElem(ii,jj,kk)%nbrs_face(seast) = kk
         GridElem(ii,jj,kk)%nbrs_wgt (seast) = CornerWgt
       end do
       end do

       ! Setup the remaining NORTH, WEST, and NW neighbors
       !---------------------------------------------------
       do jj=1,ne-1
       do ii=2,ne
         nbrs_used(ii,jj,kk,north) = .true.
         nbrs_used(ii,jj,kk, west) = .true.
         nbrs_used(ii,jj,kk,nwest) = .true.
             
         GridElem(ii,jj,kk)%nbrs     (north) = GridElem(ii  ,jj+1,kk)%number
         GridElem(ii,jj,kk)%nbrs_face(north) = kk
         GridElem(ii,jj,kk)%nbrs_wgt (north) = EdgeWgtP
         GridElem(ii,jj,kk)%nbrs     ( west) = GridElem(ii-1,jj  ,kk)%number
         GridElem(ii,jj,kk)%nbrs_face( west) = kk
         GridElem(ii,jj,kk)%nbrs_wgt ( west) = EdgeWgtP
         GridElem(ii,jj,kk)%nbrs     (nwest) = GridElem(ii-1,jj+1,kk)%number
         GridElem(ii,jj,kk)%nbrs_face(nwest) = kk
         GridElem(ii,jj,kk)%nbrs_wgt (nwest) = CornerWgt
       end do
       end do
    end do

    !======================
    ! west/east "belt" edges
    !======================
    do kk=1,4
    do jj=1,ne
      nbrs_used( 1,jj,kk,west) = .true.
      nbrs_used(ne,jj,kk,east) = .true.
          
      GridElem( 1,jj,kk)%nbrs     (west) = GridElem(ne,jj,modulo(kk+2,4)+1)%number
      GridElem( 1,jj,kk)%nbrs_face(west) = modulo(kk+2,4)+1
      GridElem( 1,jj,kk)%nbrs_wgt (west) = EdgeWgtP
      GridElem(ne,jj,kk)%nbrs     (east) = GridElem( 1,jj,modulo(kk  ,4)+1)%number
      GridElem(ne,jj,kk)%nbrs_face(east) = modulo(kk  ,4)+1
      GridElem(ne,jj,kk)%nbrs_wgt (east) = EdgeWgtP

      !  Special rules for corner 'edges'
      !----------------------------------
      if(jj /= 1) then
        nbrs_used( 1,jj,kk,swest) = .true.
        nbrs_used(ne,jj,kk,seast) = .true.
             
        GridElem( 1,jj,kk)%nbrs     (swest) = GridElem(ne,jj-1,modulo(kk+2,4)+1)%number
        GridElem( 1,jj,kk)%nbrs_face(swest) = modulo(kk+2,4)+1
        GridElem( 1,jj,kk)%nbrs_wgt (swest) = CornerWgt
        GridElem(ne,jj,kk)%nbrs     (seast) = GridElem( 1,jj-1,modulo(kk  ,4)+1)%number
        GridElem(ne,jj,kk)%nbrs_face(seast) = modulo(kk  ,4)+1
        GridElem(ne,jj,kk)%nbrs_wgt (seast) = CornerWgt
      endif
      if(jj /= ne) then
        nbrs_used( 1,jj,kk,nwest) = .true.
        nbrs_used(ne,jj,kk,neast) = .true.
             
        GridElem( 1,jj,kk)%nbrs     (nwest) = GridElem(ne,jj+1,modulo(kk+2,4)+1)%number
        GridElem( 1,jj,kk)%nbrs_face(nwest) = modulo(kk+2,4)+1
        GridElem( 1,jj,kk)%nbrs_wgt (nwest) = CornerWgt
        GridElem(ne,jj,kk)%nbrs     (neast) = GridElem(1 ,jj+1,modulo(kk  ,4)+1)%number
        GridElem(ne,jj,kk)%nbrs_face(neast) = modulo(kk  ,4)+1
        GridElem(ne,jj,kk)%nbrs_wgt (neast) = CornerWgt
      endif
    end do
    end do

    !==================================
    ! south edge of 1 / north edge of 5
    !==================================
    do ii=1,ne
      nbrs_used(ii, 1,1,south) = .true.
      nbrs_used(ii,ne,5,north) = .true.
              
      GridElem(ii, 1,1)%nbrs     (south) = GridElem(ii,ne,5)%number
      GridElem(ii, 1,1)%nbrs_face(south) = 5
      GridElem(ii, 1,1)%nbrs_wgt (south) = EdgeWgtP
      GridElem(ii,ne,5)%nbrs     (north) = GridElem(ii,1 ,1)%number
      GridElem(ii,ne,5)%nbrs_face(north) = 1
      GridElem(ii,ne,5)%nbrs_wgt (north) = EdgeWgtP

      !  Special rules for corner 'edges'
      !------------------------------------
      if(ii /= 1) then
        nbrs_used(ii, 1,1,swest) = .true.
        nbrs_used(ii,ne,5,nwest) = .true.
          
        GridElem(ii, 1,1)%nbrs     (swest) = GridElem(ii-1,ne,5)%number
        GridElem(ii, 1,1)%nbrs_face(swest) = 5
        GridElem(ii, 1,1)%nbrs_wgt (swest) = CornerWgt
        GridElem(ii,ne,5)%nbrs     (nwest) = GridElem(ii-1,1 ,1)%number
        GridElem(ii,ne,5)%nbrs_face(nwest) = 1
        GridElem(ii,ne,5)%nbrs_wgt (nwest) = CornerWgt
      endif
      if(ii /= ne) then
        nbrs_used(ii, 1,1,seast) = .true.
        nbrs_used(ii,ne,5,neast) = .true.
          
        GridElem(ii, 1,1)%nbrs     (seast) = GridElem(ii+1,ne,5)%number
        GridElem(ii, 1,1)%nbrs_face(seast) = 5
        GridElem(ii, 1,1)%nbrs_wgt (seast) = CornerWgt
        GridElem(ii,ne,5)%nbrs     (neast) = GridElem(ii+1,1 ,1)%number
        GridElem(ii,ne,5)%nbrs_face(neast) = 1
        GridElem(ii,ne,5)%nbrs_wgt (neast) = CornerWgt
      endif
    end do

    !==================================
    ! south edge of 2 / east edge of 5
    !==================================
    do ii=1,ne
      irev=ne+1-ii
      nbrs_used(ii, 1,2,south) = .true.
      nbrs_used(ne,ii,5, east) = .true.
       
      GridElem(ii, 1,2)%nbrs     (south) = GridElem(ne,irev,5)%number
      GridElem(ii, 1,2)%nbrs_face(south) = 5
      GridElem(ii, 1,2)%nbrs_wgt (south) = EdgeWgtP
      GridElem(ne,ii,5)%nbrs     ( east) = GridElem(irev,1 ,2)%number
      GridElem(ne,ii,5)%nbrs_face( east) = 2
      GridElem(ne,ii,5)%nbrs_wgt ( east) = EdgeWgtP

      !  Special rules for corner 'edges'
      !-------------------------------------
      if(ii /= 1) then
        nbrs_used(ii, 1,2,swest) = .true.
        nbrs_used(ne,ii,5,seast) = .true.
          
        GridElem(ii, 1,2)%nbrs     (swest) = GridElem(ne,irev+1,5)%number
        GridElem(ii, 1,2)%nbrs_face(swest) = 5
        GridElem(ii, 1,2)%nbrs_wgt (swest) = CornerWgt
        GridElem(ne,ii,5)%nbrs     (seast) = GridElem(irev+1,1 ,2)%number
        GridElem(ne,ii,5)%nbrs_face(seast) = 2
        GridElem(ne,ii,5)%nbrs_wgt (seast) = CornerWgt
      endif
      if(ii /= ne) then
        nbrs_used(ii, 1,2,seast) = .true.
        nbrs_used(ne,ii,5,neast) = .true.
          
        GridElem(ii, 1,2)%nbrs     (seast) = GridElem(ne,irev-1,5)%number
        GridElem(ii, 1,2)%nbrs_face(seast) = 5
        GridElem(ii, 1,2)%nbrs_wgt (seast) = CornerWgt
        GridElem(ne,ii,5)%nbrs     (neast) = GridElem(irev-1,1 ,2)%number
        GridElem(ne,ii,5)%nbrs_face(neast) = 2
        GridElem(ne,ii,5)%nbrs_wgt (neast) = CornerWgt
      endif
    enddo

    !==================================
    ! south edge of 3 / south edge of 5
    !==================================
    do ii=1,ne
      irev=ne+1-ii
      nbrs_used(ii,1,3,south) = .true.
      nbrs_used(ii,1,5,south) = .true.
       
      GridElem(ii,1,3)%nbrs     (south) = GridElem(irev,1,5)%number
      GridElem(ii,1,3)%nbrs_face(south) = 5
      GridElem(ii,1,3)%nbrs_wgt (south) = EdgeWgtP
      GridElem(ii,1,5)%nbrs     (south) = GridElem(irev,1,3)%number
      GridElem(ii,1,5)%nbrs_face(south) = 3
      GridElem(ii,1,5)%nbrs_wgt (south) = EdgeWgtP

      !  Special rules for corner 'edges'
      !-----------------------------------
      if(ii /= 1) then
        nbrs_used(ii,1,3,swest) = .true.
        nbrs_used(ii,1,5,swest) = .true.
          
        GridElem(ii,1,3)%nbrs     (swest) = GridElem(irev+1,1,5)%number
        GridElem(ii,1,3)%nbrs_face(swest) = 5
        GridElem(ii,1,3)%nbrs_wgt (swest) = CornerWgt
        GridElem(ii,1,5)%nbrs     (swest) = GridElem(irev+1,1,3)%number
        GridElem(ii,1,5)%nbrs_face(swest) = 3
        GridElem(ii,1,5)%nbrs_wgt (swest) = CornerWgt
      endif
      if(ii /= ne) then
        nbrs_used(ii,1,3,seast) = .true.
        nbrs_used(ii,1,5,seast) = .true.
          
        GridElem(ii,1,3)%nbrs     (seast) = GridElem(irev-1,1,5)%number
        GridElem(ii,1,3)%nbrs_face(seast) = 5
        GridElem(ii,1,3)%nbrs_wgt (seast) = CornerWgt
        GridElem(ii,1,5)%nbrs     (seast) = GridElem(irev-1,1,3)%number
        GridElem(ii,1,5)%nbrs_face(seast) = 3
        GridElem(ii,1,5)%nbrs_wgt (seast) = CornerWgt
      endif
    end do

    !==================================
    ! south edge of 4 / west edge of 5
    !==================================
    do ii=1,ne
      irev=ne+1-ii
      nbrs_used(ii, 1,4,south) = .true.
      nbrs_used( 1,ii,5, west) = .true.
       
      GridElem(ii, 1,4)%nbrs     (south) = GridElem(1,ii,5)%number
      GridElem(ii, 1,4)%nbrs_face(south) = 5
      GridElem(ii, 1,4)%nbrs_wgt (south) = EdgeWgtP
      GridElem( 1,ii,5)%nbrs     ( west) = GridElem(ii,1,4)%number
      GridElem( 1,ii,5)%nbrs_face( west) = 4
      GridElem( 1,ii,5)%nbrs_wgt ( west) = EdgeWgtP

      !  Special rules for corner 'edges'
      !------------------------------------
      if(ii /= 1) then
        nbrs_used(ii, 1,4,swest) = .true.
        nbrs_used( 1,ii,5,swest) = .true.
          
        GridElem(ii, 1,4)%nbrs     (swest) = GridElem(1,ii-1,5)%number
        GridElem(ii, 1,4)%nbrs_face(swest) = 5
        GridElem(ii, 1,4)%nbrs_wgt (swest) = CornerWgt
        GridElem( 1,ii,5)%nbrs     (swest) = GridElem(ii-1,1,4)%number
        GridElem( 1,ii,5)%nbrs_face(swest) = 4
        GridElem( 1,ii,5)%nbrs_wgt (swest) = CornerWgt
      endif
      if(ii /= ne) then
        nbrs_used(ii, 1,4,seast) = .true.
        nbrs_used( 1,ii,5,nwest) = .true.
          
        GridElem(ii, 1,4)%nbrs     (seast) = GridElem(1,ii+1,5)%number
        GridElem(ii, 1,4)%nbrs_face(seast) = 5
        GridElem(ii, 1,4)%nbrs_wgt (seast) = CornerWgt
        GridElem( 1,ii,5)%nbrs     (nwest) = GridElem(ii+1,1,4)%number
        GridElem( 1,ii,5)%nbrs_face(nwest) = 4
        GridElem( 1,ii,5)%nbrs_wgt (nwest) = CornerWgt
      endif
    end do

    !==================================
    ! north edge of 1 / south edge of 6
    !==================================
    do ii=1,ne
      nbrs_used(ii,ne,1,north) = .true.
      nbrs_used(ii, 1,6,south) = .true.
       
      GridElem(ii,ne,1)%nbrs     (north) = GridElem(ii, 1,6)%number
      GridElem(ii,ne,1)%nbrs_face(north) = 6
      GridElem(ii,ne,1)%nbrs_wgt (north) = EdgeWgtP
      GridElem(ii, 1,6)%nbrs     (south) = GridElem(ii,ne,1)%number
      GridElem(ii, 1,6)%nbrs_face(south) = 1
      GridElem(ii, 1,6)%nbrs_wgt (south) = EdgeWgtP

      !  Special rules for corner 'edges'
      !---------------------------------
      if(ii /= 1) then
        nbrs_used(ii,ne,1,nwest) = .true.
        nbrs_used(ii, 1,6,swest) = .true.
          
        GridElem(ii,ne,1)%nbrs     (nwest) = GridElem(ii-1, 1,6)%number
        GridElem(ii,ne,1)%nbrs_face(nwest) = 6
        GridElem(ii,ne,1)%nbrs_wgt (nwest) = CornerWgt
        GridElem(ii, 1,6)%nbrs     (swest) = GridElem(ii-1,ne,1)%number
        GridElem(ii, 1,6)%nbrs_face(swest) = 1
        GridElem(ii, 1,6)%nbrs_wgt (swest) = CornerWgt
      endif
      if(ii /= ne) then
        nbrs_used(ii,ne,1,neast) = .true.
        nbrs_used(ii, 1,6,seast) = .true.
          
        GridElem(ii,ne,1)%nbrs     (neast) = GridElem(ii+1, 1,6)%number
        GridElem(ii,ne,1)%nbrs_face(neast) = 6
        GridElem(ii,ne,1)%nbrs_wgt (neast) = CornerWgt
        GridElem(ii, 1,6)%nbrs     (seast) = GridElem(ii+1,ne,1)%number
        GridElem(ii, 1,6)%nbrs_face(seast) = 1
        GridElem(ii, 1,6)%nbrs_wgt (seast) = CornerWgt
      endif
    end do

    !==================================
    ! north edge of 2 / east edge of 6
    !==================================
    do ii=1,ne
      nbrs_used(ii,ne,2,north) = .true.
      nbrs_used(ne,ii,6, east) = .true.
       
      GridElem(ii,ne,2)%nbrs     (north) = GridElem(ne,ii,6)%number
      GridElem(ii,ne,2)%nbrs_face(north) = 6
      GridElem(ii,ne,2)%nbrs_wgt (north) = EdgeWgtP
      GridElem(ne,ii,6)%nbrs     ( east) = GridElem(ii,ne,2)%number
      GridElem(ne,ii,6)%nbrs_face( east) = 2
      GridElem(ne,ii,6)%nbrs_wgt ( east) = EdgeWgtP

      !  Special rules for corner 'edges'
      !-----------------------------------
      if(ii /= 1) then
        nbrs_used(ii,ne,2,nwest) = .true.
        nbrs_used(ne,ii,6,seast) = .true.
          
        GridElem(ii,ne,2)%nbrs     (nwest) = GridElem(ne,ii-1,6)%number
        GridElem(ii,ne,2)%nbrs_face(nwest) = 6
        GridElem(ii,ne,2)%nbrs_wgt (nwest) = CornerWgt
        GridElem(ne,ii,6)%nbrs     (seast) = GridElem(ii-1,ne,2)%number
        GridElem(ne,ii,6)%nbrs_face(seast) = 2
        GridElem(ne,ii,6)%nbrs_wgt (seast) = CornerWgt
      endif
      if(ii /= ne) then
        nbrs_used(ii,ne,2,neast) = .true.
        nbrs_used(ne,ii,6,neast) = .true.
          
        GridElem(ii,ne,2)%nbrs     (neast) = GridElem(ne,ii+1,6)%number
        GridElem(ii,ne,2)%nbrs_face(neast) = 6
        GridElem(ii,ne,2)%nbrs_wgt (neast) = CornerWgt
        GridElem(ne,ii,6)%nbrs     (neast) = GridElem(ii+1,ne,2)%number
        GridElem(ne,ii,6)%nbrs_face(neast) = 2
        GridElem(ne,ii,6)%nbrs_wgt (neast) = CornerWgt
      endif
    end do

    !===================================
    ! north edge of 3 / north edge of 6
    !===================================
    do ii=1,ne
      irev=ne+1-ii
      nbrs_used(ii,ne,3,north) = .true.
      nbrs_used(ii,ne,6,north) = .true.
       
      GridElem(ii,ne,3)%nbrs     (north) = GridElem(irev,ne,6)%number
      GridElem(ii,ne,3)%nbrs_face(north) = 6
      GridElem(ii,ne,3)%nbrs_wgt (north) = EdgeWgtP
      GridElem(ii,ne,6)%nbrs     (north) = GridElem(irev,ne,3)%number
      GridElem(ii,ne,6)%nbrs_face(north) = 3
      GridElem(ii,ne,6)%nbrs_wgt (north) = EdgeWgtP

      !  Special rules for corner 'edges'
      !------------------------------------
      if(ii /= 1) then
        nbrs_used(ii,ne,3,nwest) = .true.
        nbrs_used(ii,ne,6,nwest) = .true.
          
        GridElem(ii,ne,3)%nbrs     (nwest) = GridElem(irev+1,ne,6)%number
        GridElem(ii,ne,3)%nbrs_face(nwest) = 6
        GridElem(ii,ne,3)%nbrs_wgt (nwest) = CornerWgt
        GridElem(ii,ne,6)%nbrs     (nwest) = GridElem(irev+1,ne,3)%number
        GridElem(ii,ne,6)%nbrs_face(nwest) = 3
        GridElem(ii,ne,6)%nbrs_wgt (nwest) = CornerWgt
      endif
      if(ii /= ne) then
        nbrs_used(ii,ne,3,neast) = .true.
        nbrs_used(ii,ne,6,neast) = .true.
          
        GridElem(ii,ne,3)%nbrs     (neast) = GridElem(irev-1,ne,6)%number
        GridElem(ii,ne,3)%nbrs_face(neast) = 6
        GridElem(ii,ne,3)%nbrs_wgt (neast) = CornerWgt
        GridElem(ii,ne,6)%nbrs     (neast) = GridElem(irev-1,ne,3)%number
        GridElem(ii,ne,6)%nbrs_face(neast) = 3
        GridElem(ii,ne,6)%nbrs_wgt (neast) = CornerWgt
      endif
    end do

    !===================================
    ! north edge of 4 / west edge of 6
    !===================================
    do ii=1,ne
      irev=ne+1-ii
      nbrs_used(ii,ne,4,north) = .true.
      nbrs_used( 1,ii,6, west) = .true.
       
      GridElem(ii,ne,4)%nbrs     (north) = GridElem(1,irev,6)%number
      GridElem(ii,ne,4)%nbrs_face(north) = 6
      GridElem(ii,ne,4)%nbrs_wgt (north) = EdgeWgtP
      GridElem( 1,ii,6)%nbrs     ( west) = GridElem(irev,ne,4)%number
      GridElem( 1,ii,6)%nbrs_face( west) = 4
      GridElem( 1,ii,6)%nbrs_wgt ( west) = EdgeWgtP

      !  Special rules for corner 'edges'
      !----------------------------------
      if(ii /= 1) then
        nbrs_used(ii,ne,4,nwest) = .true.
        nbrs_used( 1,ii,6,swest) = .true.
          
        GridElem(ii,ne,4)%nbrs     (nwest) = GridElem(1,irev+1,6)%number
        GridElem(ii,ne,4)%nbrs_face(nwest) = 6
        GridElem(ii,ne,4)%nbrs_wgt (nwest) = CornerWgt
        GridElem( 1,ii,6)%nbrs     (swest) = GridElem(irev+1,ne,4)%number
        GridElem( 1,ii,6)%nbrs_face(swest) = 4
        GridElem( 1,ii,6)%nbrs_wgt (swest) = CornerWgt
      endif
      if(ii /= ne) then
        nbrs_used(ii,ne,4,neast) = .true.
        nbrs_used( 1,ii,6,nwest) = .true.
          
        GridElem(ii,ne,4)%nbrs     (neast) = GridElem(1,irev-1,6)%number
        GridElem(ii,ne,4)%nbrs_face(neast) = 6
        GridElem(ii,ne,4)%nbrs_wgt (neast) = CornerWgt
        GridElem( 1,ii,6)%nbrs     (nwest) = GridElem(irev-1,ne,4)%number
        GridElem( 1,ii,6)%nbrs_face(nwest) = 4
        GridElem( 1,ii,6)%nbrs_wgt (nwest) = CornerWgt
      endif
    end do
    
    ielem = 1                       ! Element counter
    do kk=1,6
    do jj=1,ne
    do ii=1,ne
      GridVertex(ielem)%nbrs_ptr(1) = 1
      do ll=1,8
        loc =  GridVertex(ielem)%nbrs_ptr(ll)
        if(nbrs_used(ii,jj,kk,ll)) then
          GridVertex(ielem)%nbrs          (loc) = GridElem(ii,jj,kk)%nbrs          (ll)
          GridVertex(ielem)%nbrs_face     (loc) = GridElem(ii,jj,kk)%nbrs_face     (ll)
          GridVertex(ielem)%nbrs_wgt      (loc) = GridElem(ii,jj,kk)%nbrs_wgt      (ll)
          GridVertex(ielem)%nbrs_wgt_ghost(loc) = GridElem(ii,jj,kk)%nbrs_wgt_ghost(ll)
          GridVertex(ielem)%nbrs_ptr     (ll+1) = GridVertex(ielem)%nbrs_ptr(ll)+1
        else
          GridVertex(ielem)%nbrs_ptr     (ll+1) = GridVertex(ielem)%nbrs_ptr(ll)
        endif
      end do
      GridVertex(ielem)%number           = GridElem(ii,jj,kk)%number
      GridVertex(ielem)%processor_number = 0
      GridVertex(ielem)%SpaceCurve       = GridElem(ii,jj,kk)%SpaceCurve
      ielem=ielem+1
    end do
    end do
    end do

    do kk = 1, nfaces
    do jj = 1, ne
    do ii = 1, ne
      call deallocate_gridvertex_nbrs(GridElem(ii,jj,kk))
    end do
    end do
    end do

    deallocate(Mesh     )
    deallocate(GridElem )
    deallocate(nbrs_used)

#if 0
    if(OutputFiles) then
      close(7)
      close(8)
    endif
#endif

    !=======================================
    ! Generate cube graph...
    !=======================================
#if 0
    if(OutputFiles) then
      write(9,*) nelem,2*nelem      ! METIS requires this first line
    endif
#endif

    !============================================
    !  Setup the Grid edges (topology independent)
    !============================================
    call initgridedge(GridEdge,GridVertex)

    !============================================
    !  Setup the Grid edge Indirect addresses
    !          (topology dependent)
    !============================================
    nedge = size(GridEdge)
    do ii=1,nedge
      call CubeSetupEdgeIndex(GridEdge(ii))
    end do

    ! End Routine
    !-------------
    return
  end subroutine CubeTopology
  !==================================================================


!PFC   !==================================================================
!PFC   function cube_assemble(gbl,fld,elem,par,nelemd,nelem,ielem) result(ierr)
!PFC     ! cube_assemble: Assemble the cube field element by element this routine 
!PFC     !                is assumed to be single threaded...
!PFC     !=======================================================================
!PFC     use element_mod ,only: element_t
!PFC #ifdef _MPI
!PFC     use parallel_mod,only: parallel_t, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_STATUS_SIZE, MPI_REAL8,MPI_TAG
!PFC #else
!PFC     use parallel_mod,only: parallel_t
!PFC #endif
!PFC     !
!PFC     ! Passed Variables
!PFC     !------------------
!PFC     real(kind=real_kind):: gbl(:,:,:,:)    ! global output field 
!PFC     real(kind=real_kind):: fld(:,:,:)      ! local model field  
!PFC     type(element_t)     :: elem            ! element to assemble 
!PFC     type(parallel_t)    :: par             ! parallel structure 
!PFC     integer             :: nelemd          ! number of elements on the node
!PFC     integer             :: nelem           ! number of elements on the node
!PFC     integer             :: ielem           ! local element ctr 
!PFC     integer             :: ierr            ! returned error code
!PFC     !
!PFC     ! Local Values
!PFC     !---------------
!PFC     integer ie,je,face_no
!PFC     integer ibase,jbase
!PFC     integer elem_number
!PFC     integer ne1,ne2       ! element dimensions
!PFC     integer n1,n2         ! gbl face dimensions
!PFC     integer nface         ! number of faces (must be 6)
!PFC     integer nlyr          ! number of layers
!PFC #if defined(_MPI)
!PFC     integer ectr          ! global element counter
!PFC     integer tag
!PFC     integer count         ! w/o "::", triggers PGI 3.1 F90 bug 
!PFC     integer pe
!PFC     integer status(MPI_STATUS_SIZE)
!PFC     integer mpi_err
!PFC #endif      
!PFC 
!PFC     integer :: ii,jj,kk
!PFC 
!PFC     call endrun('Because convert_gbl_index is not used cube_assemble is broken. ')
!PFC     ne1   = size(fld,1)
!PFC     ne2   = size(fld,2)
!PFC     nlyr  = size(fld,3)
!PFC     n1    = size(gbl,1)
!PFC     n2    = size(gbl,2)
!PFC     nface = size(gbl,3)
!PFC 
!PFC     !=========================
!PFC     ! Enforce certain rules...
!PFC     !=========================
!PFC     ierr=0
!PFC     if(modulo(n1,ne1) /= 0) then
!PFC       ierr=-1
!PFC       return
!PFC     endif
!PFC     if(modulo(n2,ne2) /= 0) then 
!PFC       ierr=-2
!PFC       return
!PFC     endif
!PFC     if(nface /= 6) then
!PFC       ierr=-3
!PFC       return
!PFC     endif
!PFC 
!PFC     !=========================================================
!PFC     ! Perform global assembly procedure element by element ...
!PFC     !=========================================================
!PFC     if(par%rank == par%root) then
!PFC       if(ielem <= nelemd) then
!PFC         elem_number = elem%vertex%number
!PFC 
!PFC         call convert_gbl_index(elem_number,ie,je,face_no)
!PFC         if (face_no /= elem%vertex%face_number) call endrun('Error in getting face number')
!PFC 
!PFC         ibase=ie*ne1
!PFC         jbase=je*ne2
!PFC         do kk=1,nlyr
!PFC         do jj=1,ne2
!PFC         do ii=1,ne1
!PFC           gbl(ii+ibase,jj+jbase,face_no,kk)=fld(ii,jj,kk)
!PFC         end do
!PFC         end do
!PFC         end do
!PFC       endif
!PFC #if defined(_MPI)
!PFC       if(ielem==nelemd) then
!PFC         ectr=nelemd
!PFC         do while(ectr < nelem)
!PFC           pe    = MPI_ANY_SOURCE
!PFC           tag   = MPI_ANY_TAG
!PFC           count = ne1*ne2*nlyr
!PFC           call MPI_RECV(fld(1,1,1),count,MPI_REAL8,pe,tag,par%comm,status,mpi_err) 
!PFC           elem_number = status(MPI_TAG)
!PFC           ! call convert_gbl_index(elem_number,ie,je,face_no)
!PFC           call endrun('Because convert_gbl_index is not used for neghbors, the _MPI version needs to be fixed')
!PFC 
!PFC           ibase=ie*ne1
!PFC           jbase=je*ne2
!PFC           do kk=1,nlyr
!PFC           do jj=1,ne2
!PFC           do ii=1,ne1
!PFC             gbl(ii+ibase,jj+jbase,face_no,kk)=fld(ii,jj,kk)
!PFC           end do
!PFC           end do
!PFC           end do
!PFC           ectr=ectr+1
!PFC         end do ! while(ectr < nelem)
!PFC       endif
!PFC     else
!PFC       pe    = par%root
!PFC       tag   = elem%vertex%number
!PFC       count = ne1*ne2*nlyr
!PFC       call MPI_SEND(fld(1,1,1),count,MPI_REAL8,pe,tag,par%comm,mpi_err)
!PFC #endif
!PFC     endif
!PFC 
!PFC     ! End Function
!PFC     !-------------
!PFC     return
!PFC   end function cube_assemble
!PFC   !==================================================================


  !==================================================================
  function CubeEdgeCount()  result(nedge)
    ! CubeEdgeCount: Determine the number of Grid Edges
    !
    !===================================================================
    ! 
    ! Passed Variables
    !------------------
    integer:: nedge

    if(ne==0) call endrun('Error in CubeEdgeCount: ne is zero')
    nedge = nfaces*(ne*ne*nInnerElemEdge - nCornerElemEdge)

    ! End Function
    !-------------
    return
  end function CubeEdgeCount
  !==================================================================


  !==================================================================
  function CubeElemCount()  result(nelem)
    ! CubeElemCount: Determine the number of Grid Elem
    !
    !===================================================================
    !
    ! Passed Variables
    !------------------
    integer:: nelem

    if(ne==0) call endrun('Error in CubeElemCount: ne is zero')

    nelem = nfaces*ne*ne

    ! End Function
    !-------------
    return
  end function CubeElemCount
  !==================================================================


  !==================================================================
  subroutine CubeSetupEdgeIndex(Edge)
    ! CubeSetupEdgeIndex:
    !
    !===============================================================
    use gridgraph_mod,only: gridedge_t
    use SE_Constants ,only: np
    use SE_Options   ,only: north, south, east, west, neast, seast, swest, nwest
    !
    ! Passed Variables
    !------------------
    type(GridEdge_t),target:: Edge
    !
    ! Local Values
    !--------------
    integer            :: np0,sFace,dFace
    logical            :: reverse
    integer,allocatable:: forwardV(:), forwardP(:)
    integer,allocatable:: backwardV(:), backwardP(:)
    integer            :: i,ii

    ! map to correct location - for now all on same nbr side 
    ! have same wgt, so take the first one
    !----------------------------------------------------------
    ii    = Edge%tail_face
    ii    = Edge%tail%nbrs_ptr(ii)
    np0   = Edge%tail%nbrs_wgt(ii)

    sFace = Edge%tail_face
    dFace = Edge%head_face

    ! Do not reverse the indices
    !----------------------------
    reverse=.FALSE.

    ! Under special conditions use index reversal
    !-----------------------------------------------
    if(((sFace == south).and.(dFace ==  east)).or. &
       ((sFace == east ).and.(dFace == south)).or. &
       ((sFace == north).and.(dFace ==  west)).or. &
       ((sFace == west ).and.(dFace == north)).or. &
       ((sFace == south).and.(dFace == south)).or. &
       ((sFace == north).and.(dFace == north)).or. &
       ((sFace == east ).and.(dFace ==  east)).or. &
       ((sFace == west ).and.(dFace ==  west))     ) then
      reverse     =.TRUE.
      Edge%reverse=.TRUE.
    endif

    ! End Routine
    !-------------
    return
  end subroutine CubeSetupEdgeIndex
  !==================================================================


  !==================================================================
  function ref2sphere_double(aa,bb,elem) result(sphere)         
    ! ref2sphere_double:
    ! 
    !  HOMME mapping from sphere (or other manifold) to reference element
    !  one should be able to add any mapping here.  For each new map, 
    !  an associated dmap() routine (which computes the map derivative matrix) 
    !  must also be written
    !  Note that for conservation, the parameterization of element edges must be
    !  identical for adjacent elements.  (this is violated with HOMME's default
    !  equi-angular cubed-sphere mapping for non-cubed sphere grids, hence the 
    !  need for a new map)
    !=========================================================================
    use element_mod,only: element_t
    !
    ! Passed Variables
    !-----------------
    real(kind=real_kind)   :: aa,bb
    type(element_t)        :: elem
    type(spherical_polar_t):: sphere

    if(cubed_sphere_map == 0) then
      sphere = ref2sphere_equiangular_double(aa,bb,elem%corners,elem%facenum)
    elseif(cubed_sphere_map == 1) then
!      sphere = ref2sphere_gnomonic_double(aa,bb,corners,face_no)
    elseif(cubed_sphere_map == 2) then
      sphere = ref2sphere_elementlocal_double(aa,bb,elem)
    else
      call endrun('ref2sphere_double(): bad value of cubed_sphere_map')
    endif

    ! End Function
    !-------------
    return
  end function ref2sphere_double
  !==================================================================


  !==================================================================
  function ref2sphere_longdouble(aa,bb,elem) result(sphere)         
    ! ref2sphere_longdouble:
    !
    !================================================================
    use element_mod, only : element_t
    ! 
    ! Passed Variables
    !------------------
    real(kind=longdouble_kind):: aa,bb
    type(element_t)           :: elem
    type(spherical_polar_t)   :: sphere

    if(cubed_sphere_map == 0) then
      sphere = ref2sphere_equiangular_longdouble(aa,bb,elem%corners,elem%facenum)
    elseif(cubed_sphere_map == 1) then
!      sphere = ref2sphere_gnomonic_longdouble(aa,bb,corners,face_no)
    elseif(cubed_sphere_map == 2) then
      sphere = ref2sphere_elementlocal_longdouble(aa,bb,elem)
    else
       call endrun('ref2sphere_double(): bad value of cubed_sphere_map')
    endif

    ! End Function
    !-------------
    return
  end function ref2sphere_longdouble
  !==================================================================


  !==================================================================
  function ref2sphere_equiangular_double(aa,bb,corners,face_no) result(sphere)         
    ! ref2sphere_equiangular_double:
    !
    ! map a point in the referece element to the sphere
    !=======================================================================
    !
    ! Passed Variables
    !--------------------
    real(kind=real_kind)   :: aa,bb
    type(cartesian2d_t)    :: corners(4)
    integer,intent(in)     :: face_no
    type(spherical_polar_t):: sphere
    !
    ! Local Values
    !---------------------
    real(kind=real_kind):: pi,pj,qi,qj
    type(cartesian2d_t) :: cart   

    ! map (a,b) to the [-pi/2,pi/2] equi angular cube face:  x1,x2
    ! a = gp%points(i)
    ! b = gp%points(j)
    !----------------------------------------------------------------
    pi = (1-aa)/2
    pj = (1-bb)/2
    qi = (1+aa)/2
    qj = (1+bb)/2
    cart%x = pi*pj*corners(1)%x &
           + qi*pj*corners(2)%x &
           + qi*qj*corners(3)%x &
           + pi*qj*corners(4)%x 
    cart%y = pi*pj*corners(1)%y &
           + qi*pj*corners(2)%y &
           + qi*qj*corners(3)%y &
           + pi*qj*corners(4)%y 

    ! map from [pi/2,pi/2] equ angular cube face to sphere:   
    !------------------------------------------------------
    sphere=projectpoint(cart,face_no)

    ! End Function
    !-------------
    return
  end function ref2sphere_equiangular_double
  !==================================================================


  !==================================================================
  function ref2sphere_equiangular_longdouble(aa,bb,corners,face_no) result(sphere)         
    ! ref2sphere_equiangular_longdouble:
    !
    ! map a point in the referece element to the sphere
    !==========================================================================
    !
    ! Passed Variables
    !-------------------
    real(kind=longdouble_kind):: aa,bb
    type(cartesian2d_t)       :: corners(4)
    integer,intent(in)        :: face_no
    type(spherical_polar_t)   :: sphere
    !
    ! Local Values
    !-------------
    real(kind=real_kind):: pi,pj,qi,qj
    type(cartesian2d_t) :: cart   

    ! map (a,b) to the [-pi/2,pi/2] equi angular cube face:  x1,x2
    ! a = gp%points(i)
    ! b = gp%points(j)
    !--------------------------------------------------------------
    pi = (1-aa)/2
    pj = (1-bb)/2
    qi = (1+aa)/2
    qj = (1+bb)/2
    cart%x = pi*pj*corners(1)%x &
           + qi*pj*corners(2)%x &
           + qi*qj*corners(3)%x &
           + pi*qj*corners(4)%x 
    cart%y = pi*pj*corners(1)%y &
           + qi*pj*corners(2)%y &
           + qi*qj*corners(3)%y &
           + pi*qj*corners(4)%y 

    ! map from [pi/2,pi/2] equ angular cube face to sphere:   
    !-----------------------------------------------------
    sphere=projectpoint(cart,face_no)

    ! End Function
    !-------------
    return
  end function ref2sphere_equiangular_longdouble
  !==================================================================


  !==================================================================
  function ref2sphere_elementlocal_double(aa,bb,elem) result(sphere)
    ! ref2sphere_elementlocal_double:
    !
    ! ELEMENT LOCAL MAP (DOES NOT USE CUBE FACES)
    ! unlike gnomonic equiangular map, this map will map all straight lines to
    ! great circle arcs
    !
    ! map a point in the referece element to the quad on the sphere by a
    ! general map, without using faces the map works this way: first, fix
    ! a coordinate (say, X). Map 4 corners of the ref element (corners are
    ! (-1,-1),(-1,1),(1,1), and (1,-1)) into 4 X-components of the quad in
    ! physical space via a bilinear map. Do so for Y and Z components as
    ! well. It produces a map: Ref element (\xi, \eta) ---> A quad in XYZ
    ! (ess, a piece of a twisted plane) with vertices of our target quad.  though
    ! the quad lies in a plane and not on the sphere manifold, its
    ! vertices belong to the sphere (by initial conditions). The last step
    ! is to utilize a map (X,Y,X) --> (X,Y,Z)/SQRT(X**2+Y**2+Z**2) to
    ! project the quad to the unit sphere.
    !============================================================================
    use element_mod,only: element_t
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind)   :: aa,bb
    type(element_t)        :: elem
    type(spherical_polar_t):: sphere
    !
    ! Local Values
    !---------------
    real(kind=real_kind):: qq(4)

    qq(1)=(1-aa)*(1-bb)
    qq(2)=(1+aa)*(1-bb)
    qq(3)=(1+aa)*(1+bb)
    qq(4)=(1-aa)*(1+bb)
    qq=qq/4.0d0

    sphere=ref2sphere_elementlocal_q(qq,elem%corners3D)

    ! End Function
    !-------------
    return
  end function ref2sphere_elementlocal_double
  !==================================================================


  !==================================================================
  function ref2sphere_elementlocal_longdouble(aa,bb,elem) result(sphere)
    ! ref2sphere_elementlocal_longdouble:
    !
    !======================================================================
    use element_mod,only: element_t
    !
    ! Passed Variables
    !-------------------
    real(kind=longdouble_kind):: aa,bb
    type(element_t)           :: elem
    type(spherical_polar_t)   :: sphere
    !
    ! Local Values
    !---------------
    real(kind=real_kind):: qq(4)

    qq(1)=(1-aa)*(1-bb)
    qq(2)=(1+aa)*(1-bb)
    qq(3)=(1+aa)*(1+bb)
    qq(4)=(1-aa)*(1+bb)
    qq=qq/4.0d0

    sphere=ref2sphere_elementlocal_q(qq,elem%corners3D)

    ! End Function
    !-------------
    return
  end function ref2sphere_elementlocal_longdouble
  !==================================================================


  !==================================================================
  function ref2sphere_elementlocal_q(qq,corners) result(sphere)
    ! ref2sphere_elementlocal_q:
    !
    !================================================================
    !
    ! Passed Variables
    !---------------------
    real(kind=real_kind)   :: qq(4)
    type(cartesian3d_t)    :: corners(4)
    type(spherical_polar_t):: sphere
    !
    ! Local Values
    !-----------------
    type (cartesian3d_t):: cart   
    real(kind=real_kind):: cc(3,4), xx(3), rr

    integer ii

    ! 3D corners fo the quad
    !------------------------
    cc(1,1)=corners(1)%x
    cc(2,1)=corners(1)%y
    cc(3,1)=corners(1)%z
    cc(1,2)=corners(2)%x
    cc(2,2)=corners(2)%y
    cc(3,2)=corners(2)%z
    cc(1,3)=corners(3)%x
    cc(2,3)=corners(3)%y
    cc(3,3)=corners(3)%z
    cc(1,4)=corners(4)%x
    cc(2,4)=corners(4)%y
    cc(3,4)=corners(4)%z

    ! physical point on a plane (sliced), not yet on the sphere
    !-----------------------------------------------------------
    do ii=1,3
      xx(ii)=sum(cc(ii,:)*qq(:))
    end do

    ! distance from the plane point to the origin
    !--------------------------------------------
    rr=sqrt(xx(1)**2 + xx(2)**2 + xx(3)**2)

    ! projecting the plane point to the sphere
    !------------------------------------------
    cart%x=xx(1)/rr
    cart%y=xx(2)/rr
    cart%z=xx(3)/rr

    ! XYZ coords of the point to lon/lat
    !--------------------------------------
    sphere=change_coordinates(cart)

    ! End Function
    !-------------
    return
  end function ref2sphere_elementlocal_q
  !==================================================================

end module cube_mod

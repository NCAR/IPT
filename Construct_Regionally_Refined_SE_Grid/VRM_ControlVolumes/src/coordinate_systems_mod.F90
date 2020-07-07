module coordinate_systems_mod
  !==============================================================
  ! WARNING:  When using this class be sure that you know if the
  ! cubic coordinates are on the unit cube or the [-\pi/4,\pi/4] cube
  ! and if the spherical longitude is in [0,2\pi] or [-\pi,\pi]
  !==============================================================

  ! Useful modules
  !----------------
  use SE_Constants, only : real_kind, longdouble_kind

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private

  public :: DIST_THRESHOLD
  private:: one
  private:: two

  public :: cartesian2D_t
  public :: cartesian3D_t
  public :: spherical_polar_t

  private:: copy_cart2d
  private:: eq_cart2d
  public :: distance
  private:: distance_cart2D
  private:: distance_cart2D_v
  private:: distance_cart3D
  private:: distance_cart3D_v
  public :: change_coordinates
  public :: spherical_to_cart   !CE
!  private:: spherical_to_cart
  private:: spherical_to_cart_v
  private:: cart_to_spherical
  private:: aray_to_spherical
  private:: cart_to_spherical_v
  private:: unit_face_based_cube_to_unit_sphere
  public :: cart2spherical  !CE
  public :: projectpoint        ! should be called cubedsphere2spherical
  public :: cubedsphere2cart
  public :: sphere2cubedsphere
  public :: cart2cubedsphere
  public :: cube_face_number_from_cart
  public :: cube_face_number_from_sphere
  public :: cart2cubedspherexy
  public :: sphere_tri_area
  public :: surfareaxy

  ! Parameter Values
  !------------------
  real(kind=real_kind),parameter:: DIST_THRESHOLD= 1.0D-9
  real(kind=real_kind),parameter:: one=1.0D0
  real(kind=real_kind),parameter:: two=2.0D0

  ! Type Definitions
  !--------------------
  type cartesian2D_t
    real(real_kind):: x             ! x coordinate
    real(real_kind):: y             ! y coordinate
  end type cartesian2D_t

  type cartesian3D_t
    real(real_kind):: x             ! x coordinate
    real(real_kind):: y             ! y coordinate
    real(real_kind):: z             ! z coordinate
  end type cartesian3D_t

  type spherical_polar_t
    real(real_kind):: r             ! radius
    real(real_kind):: lon           ! longitude
    real(real_kind):: lat           ! latitude
  end type spherical_polar_t

  ! Public Interfaces
  !---------------------------------------
  interface assignment ( = )
     module procedure copy_cart2d
  end interface

  interface operator( == )
     module procedure eq_cart2d
  end interface

  interface distance
     module procedure distance_cart2D
     module procedure distance_cart2D_v
     module procedure distance_cart3D
     module procedure distance_cart3D_v
  end interface

  interface change_coordinates
     module procedure spherical_to_cart_v
     module procedure spherical_to_cart
     module procedure cart_to_spherical_v
     module procedure cart_to_spherical
     module procedure aray_to_spherical
  end interface

contains

  !==================================================================
  subroutine copy_cart2d(cart2,cart1)
    ! copy_cart2d: Overload assignment operator for cartesian2D_t
    !=================================================================
    !
    ! Passed Variables
    !------------------
    type(cartesian2D_t),intent(out):: cart2
    type(cartesian2D_t),intent(in ):: cart1

    cart2%x=cart1%x
    cart2%y=cart1%y

    ! End Routine
    !-------------
    return
  end subroutine copy_cart2d
  !==================================================================


  !==================================================================
  pure function eq_cart2d(cart2,cart1) result(is_same)
    ! eq_cart2d: Overload == operator for cartesian2D_t
    !=================================================
    !
    ! Passed Variables
    !------------------
    type(cartesian2D_t),intent(in):: cart2
    type(cartesian2D_t),intent(in):: cart1
    logical                       :: is_same    

    if(distance(cart1,cart2) < DIST_THRESHOLD) then
      is_same=.true.
    else
      is_same=.false.
    endif

    ! End Function
    !-------------
    return
  end function eq_cart2d
  !==================================================================


  !==================================================================
  pure function distance_cart2D(cart1,cart2) result(dist)
    ! distance_cart2D  : scalar version
    !                    computes distance between cartesian 2D coordinates
    !======================================================================
    !
    ! Passed Variables
    !------------------
    type(cartesian2D_t),intent(in)          :: cart1
    type(cartesian2D_t),intent(in),optional :: cart2
    real(real_kind)                         :: dist

    if(present(cart2)) then
      dist = SQRT((cart1%x-cart2%x)**2 + (cart1%y-cart2%y)**2)
    else
      dist = SQRT(cart1%x*cart1%x + cart1%y*cart1%y)
    endif

    ! End Function
    !-------------
    return
  end function distance_cart2D
  !==================================================================


  !==================================================================
  pure function distance_cart2D_v(cart1,cart2) result(dist)
    ! distance_cart2D_v: vector version
    !                    computes distance between cartesian 2D coordinates
    !======================================================================
    !
    ! Passed Variables
    !------------------
    type(cartesian2D_t),intent(in)         :: cart1(:)
    type(cartesian2D_t),intent(in),optional:: cart2(:)
    real(real_kind)                        :: dist(SIZE(cart1))
    !
    ! Local Values
    !--------------
    integer ii

    if(present(cart2)) then
      forall(ii=1:SIZE(cart1)) dist(ii) = distance_cart2D(cart1(ii),cart2(ii))
    else
      forall(ii=1:SIZE(cart1)) dist(ii) = distance_cart2D(cart1(ii))
    end if

    ! End Function
    !-------------
    return
  end function distance_cart2D_v
  !==================================================================


  !==================================================================
  pure function distance_cart3D(cart1,cart2) result(dist)
    ! distance_cart3D  : scalar version
    !
    !=========================================================
    !
    ! Passed Variables
    !------------------
    type(cartesian3D_t),intent(in)         :: cart1
    type(cartesian3D_t),intent(in),optional:: cart2
    real(real_kind)                        :: dist

    if(present(cart2)) then
      dist = SQRT((cart1%x-cart2%x)**2 + (cart1%y-cart2%y)**2 + (cart1%z-cart2%z)**2)
    else
      dist = SQRT(cart1%x*cart1%x + cart1%y*cart1%y + cart1%z*cart1%z)
    endif

    ! End Function
    !-------------
    return
  end function distance_cart3D
  !==================================================================


  !==================================================================
  pure function distance_cart3D_v(cart1,cart2) result(dist)
    ! distance_cart3D_v: vector version
    !
    !=========================================================
    !
    ! Passed Variables
    !------------------
    type(cartesian3D_t),intent(in)         :: cart1(:)
    type(cartesian3D_t),intent(in),optional:: cart2(:)
    real(real_kind)                        :: dist(SIZE(cart1))
    !
    ! Local Values
    !----------------
    integer ii

    if(present(cart2)) then
      forall(ii=1:SIZE(cart1)) dist(ii) = distance_cart3D(cart1(ii),cart2(ii))
    else
      forall(ii=1:SIZE(cart1)) dist(ii) = distance_cart3D(cart1(ii))
    end if

    ! End Function
    !-------------
    return
  end function distance_cart3D_v
  !==================================================================


  !==================================================================
  pure function spherical_to_cart(sphere) result (cart)
    ! spherical_to_cart: converts spherical polar {lon,lat} to 3D cartesian {x,y,z}
    !                    on unit sphere.  Note: spherical longitude is [0,2\pi]
    !=============================================================================
    !
    ! Passed Variables
    !--------------------
    type(spherical_polar_t),intent(in):: sphere
    type(cartesian3D_t)               :: cart

    cart%x=sphere%r*COS(sphere%lat)*COS(sphere%lon)
    cart%y=sphere%r*COS(sphere%lat)*SIN(sphere%lon)
    cart%z=sphere%r*SIN(sphere%lat)

    ! End Function
    !-------------
    return
  end function spherical_to_cart
  !==================================================================


  !==================================================================
  pure function spherical_to_cart_v(sphere) result (cart)
    ! spherical_to_cart_v: converts spherical polar {lon,lat} to 3D cartesian {x,y,z}
    !                      on unit sphere.  Note: spherical longitude is [0,2\pi]
    !===============================================================================
    !
    ! Passed Variables
    !---------------------
    type(spherical_polar_t),intent(in):: sphere(:)
    type(cartesian3D_t)               :: cart(SIZE(sphere))
    !
    ! Local Values
    !---------------
    integer ii

    forall(ii=1:SIZE(sphere)) cart(ii) = spherical_to_cart(sphere(ii))

    ! End Function
    !-------------
    return
  end function spherical_to_cart_v
  !==================================================================


  !==================================================================
  pure function cart_to_spherical(cart) result (sphere)
    ! cart_to_spherical: converts 3D cartesian {x,y,z} to spherical polar {lon,lat} 
    !                    on unit sphere. Note: spherical longitude is [0,2\pi]
    !                    scalar version
    !==========================================================================
    use SE_Constants,only: DD_PI
    !
    ! Passed Variables
    !------------------
    type(cartesian3D_t),intent(in):: cart         
    type(spherical_polar_t)       :: sphere

    sphere%r   = distance(cart)
    sphere%lat = ASIN(cart%z/sphere%r)
    sphere%lon = 0

    !----------------------------------------------------------------
    ! enforce three facts:
    ! 1) lon at poles is defined to be zero
    ! 2) Grid points must be separated by about .01 Meter (on earth)
    !    from pole to be considered "not the pole".
    ! 3) range of lon is { 0<= lon < 2*pi }
    !----------------------------------------------------------------

    ! if point is away from the POLE.  distance(cart) = distance from center of earth,
    ! so this was a bug:
    !  if (distance(cart) >= DIST_THRESHOLD) then 
    !------------------------------------------------------------------------
    if(abs(abs(sphere%lat)-DD_PI/2)  >= DIST_THRESHOLD) then
      sphere%lon=ATAN2(cart%y,cart%x)
      if(sphere%lon < 0) then
        sphere%lon=sphere%lon + 2*DD_PI
      endif
    endif

    ! End Function
    !-------------
    return
  end function cart_to_spherical
  !==================================================================


  !==================================================================
  pure function aray_to_spherical(coordinates) result (sphere)
    ! aray_to_spherical:
    !
    !=============================================================
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind),intent(in):: coordinates(3)
    type(spherical_polar_t)        :: sphere
    !
    ! Local Values
    !--------------
    type(cartesian3D_t):: cart

    cart%x = coordinates(1)
    cart%y = coordinates(2)
    cart%z = coordinates(3)
    sphere = cart_to_spherical(cart) 

    ! End Function
    !-------------
    return
  end function aray_to_spherical
  !==================================================================


  !==================================================================
  pure function cart_to_spherical_v(cart) result (sphere)
    ! cart_to_spherical: converts 3D cartesian {x,y,z} to spherical polar {lon,lat} 
    !                    on unit sphere. Note: spherical longitude is [0,2\pi]
    !                    vector version
    !==========================================================================
    use SE_Constants,only: DD_PI
    !
    ! Passed Variables
    !-------------------
    type(cartesian3D_t),intent(in):: cart(:)
    type(spherical_polar_t)       :: sphere(SIZE(cart))
    !
    ! Local Values
    !---------------
    integer ii

    forall(ii=1:SIZE(cart)) sphere(ii) = cart_to_spherical(cart(ii))

    ! End Function
    !-------------
    return
  end function cart_to_spherical_v
  !==================================================================


  !==================================================================
  function unit_face_based_cube_to_unit_sphere(cart, face_no) result(sphere)
    ! unit_face_based_cube_to_unit_sphere:
    !
    ! Note: Output spherical longitude is [-pi,pi]
    !
    ! Project from a UNIT cube to a UNIT sphere.  ie, the lenght of the cube edge is 2.
    ! Face 1 of the cube touches the sphere at longitude, latitude (0,0). The negative 
    ! x axis is negative longitude (ie. going west is negative), the positive x axis
    ! is increasing longitude.  Face 1 maps the Face 1 to the lat,lon on the sphere:
    !    [-1,1] x [-1,1] => [-\pi/4,\pi/4] x [-\pi/4, \pi/4]
    !
    ! Face 2 continues with increasing longitude (ie to the east of Face 1).  
    ! The left edge of Face 2 (negative x) is the right edge of Face 1 (positive x)
    ! The latitude is the same as Face 1, but the longitude increases:
    !    [-1,1] x [-1,1] => [\pi/4, 3\pi/4] x [-\pi/4, \pi/4]

    ! Face 3 continues with increasing longitude (ie to the east of Face 2).  
    ! Face 3 is like Face 1, but the x coordinates are reversed, ie. decreasing x 
    ! is increasing longitude:
    !    [-1,1] x [-1,1]  =    [-1,0] x [-1,1] U  [0,1] x [-1,1] =>
    !            [3\pi/4,\pi] x [-\pi, -3\pi/4]
    !
    ! Face 4 finally connects Face 3 to Face 1.  Like Face 2, but wtih opposite x
    !    [-1,1] x [-1,1] => [-3\pi/4, -\pi/4] x [-\pi/4, \pi/4]
    !
    ! Face 5 is along the bottom edges of Faces 1,2,3,and 4 so the latitude goes from
    ! -\pi/4 to -\pi/2.  The tricky part is lining up the longitude.  The zero longitude
    ! must line up with the center of Face 1. ATAN2(x,1) = 0 => x = 0.  
    ! So the (0,1) point on Face 5 is the zero longitude on the sphere.  The top edge of 
    ! Face 5 is the bottom edge of Face 1. 
    ! ATAN(x,0) = \pi/2 => x = 1, so the right edge of Face 5 is the bottom of Face 2.
    ! Continueing, the bottom edge of 5 is the bottom of 3.  Left of 5 is bottom of 4.
    !
    ! Face 6 is along the top edges of Faces 1,2,3 and 4 so the latitude goes from  
    ! \pi/4 to \pi/2.   The zero longitude must line up with the center of Face 1.  
    ! This is just like Face 5, but the y axis is reversed.  So the bottom edge of Face 6
    ! is the top edge of Face 1.  The right edge of Face 6 is the top of Face 2.  The
    ! top of 6 the top of 3 and the left of 6 the top of 4.
    !==================================================================================
    use SE_Constants,only: DD_PI
    use err_exit   ,only: endrun
    !
    ! Passed Variables
    !--------------------
    type (cartesian2d_t),intent(in):: cart   ! On face_no of a unit cube
    integer             ,intent(in):: face_no 
    type (spherical_polar_t)       :: sphere
    !
    ! Local Values
    !--------------
    real(kind=real_kind):: rr !, l_inf

    ! MNL: removing check that points are on the unit cube because we allow
    ! spherical grids to map beyond the extent of the cube (though we probably
    ! should still have an upper bound for how far past the edge the element lies)
    !    l_inf = MAX(ABS(cart%x), ABS(cart%y)) 
    !    if (1.01 < l_inf) then
    !      call endrun('unit_face_based_cube_to_unit_sphere: Input not on unit cube.')
    !    end if
    !-------------------------------------------------------------------------------
    sphere%r=one
    rr = SQRT( one + (cart%x)**2 + (cart%y)**2)
    select case (face_no)
      case (1) 
             sphere%lat=ASIN((cart%y)/rr)
             sphere%lon=ATAN2(cart%x,one)
      case (2) 
             sphere%lat=ASIN((cart%y)/rr)
             sphere%lon=ATAN2(one,-cart%x)
      case (3) 
             sphere%lat=ASIN((cart%y)/rr)
             sphere%lon=ATAN2(-cart%x,-one)
      case (4) 
             sphere%lat=ASIN((cart%y)/rr)
             sphere%lon=ATAN2(-one,cart%x)
      case (5) 
             if((ABS(cart%y) > DIST_THRESHOLD).or.(ABS(cart%x) > DIST_THRESHOLD)) then
               sphere%lon = ATAN2(cart%x, cart%y )
             else
               sphere%lon = 0.0D0     ! longitude is meaningless at south pole set to 0.0
             endif
             sphere%lat=ASIN(-one/rr)
      case (6) 
             if((ABS(cart%y) > DIST_THRESHOLD).or.(ABS(cart%x) > DIST_THRESHOLD)) then
               sphere%lon = ATAN2(cart%x, -cart%y)
             else
               sphere%lon = 0.0D0     ! longitude is meaningless at north pole set to 0.0
             endif
             sphere%lat=ASIN(one/rr)
      case default
             call endrun('unit_face_based_cube_to_unit_sphere: Face number not 1 to 6.')
    end select

    if(sphere%lon < 0.0D0) then
      sphere%lon = sphere%lon + two*DD_PI
    endif

    ! End Function
    !-------------
    return
  end function unit_face_based_cube_to_unit_sphere
  !==================================================================


  !==================================================================
  function cart2spherical(xx,yy, face_no) result(sphere)
    ! cart2spherical:
    !
    ! IMPORTANT: INPUT ARE the REAL cartesian from the cube sphere
    ! Note: Output spherical longitude is [-pi,pi]
    !
    ! Project from a UNIT cube to a UNIT sphere.  ie, the lenght of the cube edge is 2.
    ! Face 1 of the cube touches the sphere at longitude, latitude (0,0). The negative 
    ! x axis is negative longitude (ie. going west is negative), the positive x axis
    ! is increasing longitude.  Face 1 maps the Face 1 to the lat,lon on the sphere:
    !    [-1,1] x [-1,1] => [-\pi/4,\pi/4] x [-\pi/4, \pi/4]
    !
    ! Face 2 continues with increasing longitude (ie to the east of Face 1).  
    ! The left edge of Face 2 (negative x) is the right edge of Face 1 (positive x)
    ! The latitude is the same as Face 1, but the longitude increases:
    !    [-1,1] x [-1,1] => [\pi/4, 3\pi/4] x [-\pi/4, \pi/4]
    !
    ! Face 3 continues with increasing longitude (ie to the east of Face 2).  
    ! Face 3 is like Face 1, but the x coordinates are reversed, ie. decreasing x 
    ! is increasing longitude:
    !    [-1,1] x [-1,1]  =    [-1,0] x [-1,1] U  [0,1] x [-1,1] =>
    !            [3\pi/4,\pi] x [-\pi, -3\pi/4]
    !
    ! Face 4 finally connects Face 3 to Face 1.  Like Face 2, but wtih opposite x
    !    [-1,1] x [-1,1] => [-3\pi/4, -\pi/4] x [-\pi/4, \pi/4]
    !
    ! Face 5 is along the bottom edges of Faces 1,2,3,and 4 so the latitude goes from
    ! -\pi/4 to -\pi/2.  The tricky part is lining up the longitude.  The zero longitude
    ! must line up with the center of Face 1. ATAN2(x,1) = 0 => x = 0.  
    ! So the (0,1) point on Face 5 is the zero longitude on the sphere.  The top edge of 
    ! Face 5 is the bottom edge of Face 1. 
    ! ATAN(x,0) = \pi/2 => x = 1, so the right edge of Face 5 is the bottom of Face 2.
    ! Continueing, the bottom edge of 5 is the bottom of 3.  Left of 5 is bottom of 4.
    !
    ! Face 6 is along the top edges of Faces 1,2,3 and 4 so the latitude goes from  
    ! \pi/4 to \pi/2.   The zero longitude must line up with the center of Face 1.  
    ! This is just like Face 5, but the y axis is reversed.  So the bottom edge of Face 6
    ! is the top edge of Face 1.  The right edge of Face 6 is the top of Face 2.  The
    ! top of 6 the top of 3 and the left of 6 the top of 4.
    !===================================================================================
    use SE_Constants,only: DD_PI
    use err_exit   ,only: endrun
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind),intent(in):: xx,yy   ! On face_no of a unit cube
    integer             ,intent(in):: face_no 
    type(spherical_polar_t)        :: sphere
    !
    ! Local Values
    !----------------
    real(kind=real_kind):: rr !, l_inf

    ! MNL: removing check that points are on the unit cube because we allow
    ! spherical grids to map beyond the extent of the cube (though we probably
    ! should still have an upper bound for how far past the edge the element lies)
    !    l_inf = MAX(ABS(cart%x), ABS(cart%y)) 
    !    if (1.01 < l_inf) then
    !      call endrun('unit_face_based_cube_to_unit_sphere: Input not on unit cube.')
    !    end if
    !------------------------------------------------------------------------------------

    sphere%r=one
    rr = SQRT( one + xx**2 + yy**2)
    select case (face_no)
      case (1) 
             sphere%lat=ASIN(yy/rr)
             sphere%lon=ATAN2(xx,one)
      case (2) 
             sphere%lat=ASIN(yy/rr)
             sphere%lon=ATAN2(one,-xx)
      case (3) 
             sphere%lat=ASIN(yy/rr)
             sphere%lon=ATAN2(-xx,-one)
      case (4) 
             sphere%lat=ASIN(yy/rr)
             sphere%lon=ATAN2(-one,xx)
      case (5) 
             if((ABS(yy) > DIST_THRESHOLD).or.(ABS(xx) > DIST_THRESHOLD)) then
               sphere%lon=ATAN2(xx,yy)
             else
               sphere%lon= 0.0D0     ! longitude is meaningless at south pole set to 0.0
             end if
             sphere%lat=ASIN(-one/rr)
      case (6) 
             if((ABS(yy) > DIST_THRESHOLD).or.(ABS(xx) > DIST_THRESHOLD)) then
               sphere%lon = ATAN2(xx, -yy)
             else
               sphere%lon= 0.0D0     ! longitude is meaningless at north pole set to 0.0
             end if
             sphere%lat=ASIN(one/rr)
      case default
             call endrun('unit_face_based_cube_to_unit_sphere: Face number not 1 to 6.')
    end select

    if(sphere%lon < 0.0D0) then
      sphere%lon = sphere%lon + two*DD_PI
    endif

    ! End Function
    !-------------
    return
  end function cart2spherical
  !==================================================================


  !==================================================================
  function projectpoint(cartin, face_no) result(sphere)         
    ! projectpoint:
    !
    ! Note: Output spherical longitude is [-pi,pi]
    !
    ! Projection from a [-pi/4, \pi/4] sized cube.  
    ! This will be checked because unit_face_based_cube_to_unit_sphere checks the ranges.
    ! See unit_face_based_cube_to_unit_sphere for documentation.
    !=======================================================================
    !
    ! Passed Variables
    !-----------------
    type(cartesian2d_t),intent(in):: cartin   
    integer            ,intent(in):: face_no
    type(spherical_polar_t)       :: sphere
    !
    ! Local Values
    !----------------
    type(cartesian2d_t):: cart   

    !ASC  This is X and Y and not xhi eta ...
    !--------------------------------------
    cart%x = TAN(cartin%x)
    cart%y = TAN(cartin%y)

    sphere = unit_face_based_cube_to_unit_sphere(cart, face_no)

    ! End Function
    !-------------
    return
  end function projectpoint
  !==================================================================


  !==================================================================
  function cubedsphere2cart(cartin, face_no) result(cart)
    ! cubedsphere2cart:
    !
    ! takes a 2D point on a face of the cube of size [-\pi/4, \pi/4] and projects it 
    ! onto a 3D point on a cube of size [-1,1] in R^3
    !=====================================================================
    !
    ! Passed Variables
    !------------------
    type(cartesian2d_t),intent(in):: cartin   ! assumed to be cartesian coordinates of cube
    integer            ,intent(in):: face_no
    type(cartesian3D_t)           :: cart

    cart = spherical_to_cart(projectpoint(cartin, face_no))

    ! End Function
    !-------------
    return
  end function cubedsphere2cart
  !==================================================================


  !==================================================================
  pure function sphere2cubedsphere(sphere, face_no) result(cart)
    ! sphere2cubedsphere: onto a cube of size [-\pi/2,\pi/2] in R^3 the spherical 
    !                     longitude can be either in [0,2\pi] or [-\pi,\pi]
    !==============================================================================
    use SE_Constants,only: DD_PI
    !
    ! Passed Variables
    !---------------------
    type(spherical_polar_t),intent(in):: sphere
    integer                ,intent(in):: face_no
    type(cartesian2d_t)               :: cart
    !
    ! Local Values
    !---------------
    real(kind=real_kind):: xp,yp
    real(kind=real_kind):: lat,lon
    real(kind=real_kind):: pi,twopi, pi2, pi3, pi4

    lat  = sphere%lat
    lon  = sphere%lon
    pi   = DD_PI 
    twopi= 2.0D0 * pi
    pi2  = pi * 0.5D0 
    pi3  = pi * 1.5D0               
    pi4  = pi * 0.25D0

    select case (face_no)
      case  (1) 
              xp = lon
              if(pi < lon) xp=lon - twopi !if lon in [0,2\pi]
              yp = atan(tan(lat)/cos(xp))
      case  (2) 
              xp = lon - pi2
              yp = atan(tan(lat)/cos(xp))
      case  (3) 
              xp = lon - pi
              if(lon < 0) xp=lon + pi  !if lon in [0,2\pi]
              yp = atan(tan(lat)/cos(xp))
      case  (4) 
              xp = lon - pi3
              if(lon < 0) xp=lon + pi2  !if lon in [0,2\pi]
              yp = atan(tan(lat)/cos(xp))
      case  (5) 
              xp = atan(-sin(lon)/tan(lat))
              yp = atan(-cos(lon)/tan(lat))
      case  (6) 
              xp = atan( sin(lon)/tan(lat))
              yp = atan(-cos(lon)/tan(lat))
    end select

    ! coordinates on the cube:
    !---------------------------
    cart%x = xp
    cart%y = yp
    
    ! End Function
    !-------------
    return
  end function sphere2cubedsphere 
  !==================================================================


  !==================================================================
  pure function cart2cubedsphere(cart3D, face_no) result(cart)
    ! cart2cubedsphere:
    !
    ! Go from an arbitrary sized cube in 3D 
    ! to a [-\pi/4,\pi/4] sized cube with (face,2d) coordinates.  
    !
    !                        Z
    !                        |
    !                        |
    !                        |
    !                        |
    !                        ---------------Y
    !                       /
    !                      /
    !                     /
    !                    /
    !                   X
    !
    ! NOTE: Face 1 =>  X positive constant face of cube
    !       Face 2 =>  Y positive constant face of cube
    !       Face 3 =>  X negative constant face of cube
    !       Face 4 =>  Y negative constant face of cube
    !       Face 5 =>  Z negative constant face of cube
    !       Face 6 =>  Z positive constant face of cube
    !================================================================
    !
    ! Passed Variables
    !---------------------
    type(cartesian3D_t),intent(in):: cart3d
    integer            ,intent(in):: face_no
    type(cartesian2d_t)           :: cart   
    !
    ! Local Values
    !---------------
    real(kind=real_kind):: xx,yy

    select case (face_no) 
      case (1)
             xx =  cart3D%y/cart3D%x
             yy =  cart3D%z/cart3D%x
      case (2)
             xx = -cart3D%x/cart3D%y
             yy =  cart3D%z/cart3D%y
      case (3)
             xx =  cart3D%y/cart3D%x
             yy = -cart3D%z/cart3D%x
      case (4)
             xx = -cart3D%x/cart3D%y
             yy = -cart3D%z/cart3D%y
      case (5)
             xx = -cart3D%y/cart3D%z
             yy = -cart3D%x/cart3D%z
      case (6)
             xx =  cart3D%y/cart3D%z
             yy = -cart3D%x/cart3D%z
    end select

    cart%x = ATAN(xx)
    cart%y = ATAN(yy)

    ! End Function
    !-------------
    return
  end function cart2cubedsphere
  !==================================================================


  !==================================================================
  pure function cube_face_number_from_cart(cart) result(face_no)
    ! cube_face_number_from_cart:
    !
    ! This function divides three dimentional space up into 
    ! six sectors.  These sectors are then considered as the
    ! faces of the cube.  It should work for any (x,y,z) coordinate
    ! if on a sphere or on a cube.
    !=================================================================
    !
    ! Passed Variables
    !------------------
    type(cartesian3D_t),intent(in):: cart  
    integer                       :: face_no
    !
    ! Local Values
    !----------------
    real(real_kind):: xx,yy,zz

    xx=cart%x
    yy=cart%y
    zz=cart%z

    ! Divide the X-Y plane into for quadrants of 
    ! [-\pi/2,\pi/2], [\pi/2,3\pi/2], .....
    ! based on the lines X=Y and X=-Y.  This divides
    ! 3D space up into four sections.  Doing the same
    ! for the XZ and YZ planes divides space into six
    ! sections.  Can also be thought of as conic sections
    ! in the L_infinity norm.  
    !------------------------------------------------------

    if((yy < xx).and.(yy > -xx)) then      ! x>0, Face 1,5 or 6
      if(zz > xx) then
        face_no=6  ! north pole
      elseif(zz < -xx) then
        face_no=5  ! south pole
      else
        face_no=1
      endif
    elseif((yy > xx).and.(yy < -xx)) then  ! x<0
      if(zz > -xx) then
        face_no=6  ! north pole
      elseif(zz < xx) then
        face_no=5  ! south pole
      else 
        face_no=3
      endif
    elseif((yy > xx).and.(yy > -xx)) then  ! y>0
      if(zz > yy) then
        face_no=6  ! north pole
      elseif(zz < -yy) then
        face_no=5  ! south pole
      else 
        face_no=2
      endif
    elseif((yy < xx).and.(yy < -xx)) then  ! y<0
      if(zz > -yy) then
        face_no=6  ! north pole
      elseif(zz < yy) then
        face_no=5  ! south pole
      else 
        face_no=4
      endif
    else
      ! abs(y) = abs(x).  point is on cube edge, or on face 5 or 6:
      !-------------------------------------------------------------
      if(abs(xx) < zz) then
        face_no=6  ! north pole
      elseif(zz < -abs(xx)) then
        face_no=5  ! south pole
      elseif((0 < xx).and.(0 < yy)) then
        face_no=1
      elseif((xx < 0).and.(0 < yy)) then
        face_no=2
      elseif((xx < 0).and.(yy < 0)) then
        face_no=3
      else
        face_no=4
      endif
    endif
   
    ! End Function
    !-------------
    return
   end function cube_face_number_from_cart
  !==================================================================


  !==================================================================
  pure function cube_face_number_from_sphere(sphere) result(face_no)
    ! cube_face_number_from_sphere:
    !
    ! This could be done directly by using the lon, lat coordinates,
    ! but call cube_face_number_from_cart just so that there is one place
    ! to do the conversions and they are all consistant.
    !=================================================================
    !
    ! Passed Variables
    !-------------------
    type(spherical_polar_t),intent(in):: sphere
    integer                           :: face_no 
    !
    ! Local Values
    !-----------------
    type(cartesian3d_t):: cart   

    cart    = spherical_to_cart(sphere)
    face_no = cube_face_number_from_cart(cart) 

    ! End Function
    !-------------
    return
  end function cube_face_number_from_sphere
  !==================================================================


  !==================================================================
  subroutine cart2cubedspherexy(cart3d,face_no,cartxy)
    ! cart2cubedspherexy:
    !
    ! CE, need real (cartesian) xy coordinates on the cubed sphere
    !===============================================================
    !
    ! Passed Variables
    !----------------------
    type(cartesian3D_t),intent(in ):: cart3d
    integer            ,intent(in ):: face_no
    type(cartesian2d_t),intent(out):: cartxy   

    ! a (half length of a cube side) is supposed to be 1
    !-------------------------------------------------
    select case (face_no)
      case (1)
             cartxy%x =  cart3D%y/cart3D%x
             cartxy%y =  cart3D%z/cart3D%x
      case (2)
             cartxy%x = -cart3D%x/cart3D%y
             cartxy%y =  cart3D%z/cart3D%y
      case (3)
             cartxy%x =  cart3D%y/cart3D%x
             cartxy%y = -cart3D%z/cart3D%x
      case (4)
             cartxy%x = -cart3D%x/cart3D%y
             cartxy%y = -cart3D%z/cart3D%y
      case (5)       !bottom face
             cartxy%x = -cart3D%y/cart3D%z
             cartxy%y = -cart3D%x/cart3D%z
      case (6)        !top face
             cartxy%x =  cart3D%y/cart3D%z
             cartxy%y = -cart3D%x/cart3D%z
    end select

    ! End Routine
    !-------------
    return
  end subroutine cart2cubedspherexy
  !==================================================================


  !==================================================================
  subroutine sphere_tri_area( v1, v2, v3, area )
    ! sphere_tri_area:
    !
    !  input: v1(3),v2(3),v3(3)  cartesian coordinates of triangle
    !  output: area
    !  based on formulas in STRI_QUAD:
    !  http://people.sc.fsu.edu/~burkardt/f_src/stri_quad/stri_quad.html
    !===================================================================
    use SE_Constants, only : DD_PI
    !
    ! Passed Variables
    !-----------------
    type(cartesian3D_t) :: v1,v2,v3
    real(kind=real_kind):: area
    !
    ! Local Values
    !----------------
    real(kind=real_kind):: aa,bb,cc,al,bl,cl,sina,sinb,sinc,sins,a1,b1,c1
  
    ! compute great circle lengths
    !-------------------------------
    al = acos( v2%x * v3%x + v2%y * v3%y + v2%z * v3%z )
    bl = acos( v3%x * v1%x + v3%y * v1%y + v3%z * v1%z )
    cl = acos( v1%x * v2%x + v1%y * v2%y + v1%z * v2%z )

    ! compute angles
    !----------------
    sina = sin( (bl+cl-al)/2 )  ! sin(sl-al)
    sinb = sin( (al+cl-bl)/2 )  ! sin(sl-bl)
    sinc = sin( (al+bl-cl)/2 ) 
    sins = sin( (al+bl+cl)/2 )

#if 0
    ! apply Girard's theorem
    !-----------------------
    aa = 2*atan2(sqrt(sinb*sinc), sqrt(sins*sina)  )
    bb = 2*atan2(sqrt(sina*sinc), sqrt(sins*sinb) )
    cc = 2*atan2(sqrt(sina*sinb), sqrt(sins*sinc) )
    area = aa+bb+cc - DD_PI
#endif

    ! for small areas, formula above looses precision.  
    ! 2atan(x) + 2atan(1/x) = pi      
    ! 2atan(x) - pi = -2atan(1/x)
    !-----------------------------------------------
    aa = sqrt( (sinb*sinc) / (sins*sina) ) 
    bb = sqrt( (sina*sinc) / (sins*sinb) ) 
    cc = sqrt( (sina*sinb) / (sins*sinc) ) 

    a1 = 2*atan(aa)
    b1 = 2*atan(bb)
    c1 = 2*atan(cc)

    if((aa.gt.bb).and.(aa.gt.cc)) then
       a1 = -2*atan(1/aa)
    elseif(bb.gt.cc) then
       b1 = -2*atan(1/bb)
    else 
       c1 = -2*atan(1/cc)
    endif

    ! apply Girard's theorem
    !-----------------------
    area = a1+b1+c1  

    ! End Routine
    !-------------
    return
  end subroutine sphere_tri_area
  !==================================================================


  !==================================================================
  function surfareaxy(x1,x2,y1,y2) result(area)
    ! surfareaxy:
    !
    !CE, 5.May 2011
    !INPUT: Points in xy cubed sphere coordinates, counterclockwise
    !OUTPUT: corresponding area on the sphere
    !===============================================================
    !
    ! Passed Variables
    !---------------------
    real(kind=real_kind),intent(in):: x1, x2, y1, y2
    real(kind=real_kind)           :: area
    !
    ! Local Values
    !----------------
    real(kind=real_kind):: a1,a2,a3,a4

    ! cubed-sphere cell area, from Lauritzen & Nair MWR 2008
    ! central angles:
    ! cube face: -pi/4,-pi/4 -> pi/4,pi/4
    ! this formula gives 2   so normalize by 4pi/6 / 2 = pi/3
    ! use implementation where the nodes a counterclockwise (not as in the paper)
    !-----------------------------------------------------------------------------
    a1 = acos(-sin(atan(x1))*sin(atan(y1)))             
    a2 =-acos(-sin(atan(x2))*sin(atan(y1)))  
    a3 = acos(-sin(atan(x2))*sin(atan(y2)))              
    a4 =-acos(-sin(atan(x1))*sin(atan(y2)))         
    area = (a1+a2+a3+a4)

    ! End Function
    !-------------
    return
  end function surfareaxy
  !==================================================================

end module coordinate_systems_mod

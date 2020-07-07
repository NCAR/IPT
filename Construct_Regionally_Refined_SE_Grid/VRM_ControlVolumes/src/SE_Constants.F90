module SE_Constants
!=========================================================================
!
! Purpose: Implement constants needed for SE computations. Including
!          type definitions, physical constants, and parameters needed
!          to specify dimensions. 
!
! Revisions:
!   
!=========================================================================

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private

  ! Types
  !-------
  public:: int_kind
  public:: long_kind
  public:: log_kind
  public:: real_kind
  public:: longdouble_kind
  public:: MAX_STRING_LEN
  public:: MAX_FILE_LEN

  ! Physical Constants
  !--------------------
  public:: QQ_PI
  public:: DD_PI
  public:: rearth
  public:: g
  public:: omega
  public:: Rgas
  public:: Cp
  public:: p0
  public:: MWDAIR
  public:: Rwater_vapor
  public:: Cpwater_vapor
  public:: kappa
  public:: Rd_on_Rv
  public:: Cpd_on_Cpv
  public:: rrearth
  public:: Lc

  ! Spectral Element Dimensions
  !-----------------------------
  public:: np
  public:: nc
  public:: nip
  public:: nipm
  public:: nep
  public:: npdg
  public:: pcnst
  public:: qsize_d
  public:: ntrac_d
  public:: npsq
  public:: nlev
  public:: nlevp

  ! Set int/float kinds
  !---------------------
  integer(kind=4),parameter:: int_kind       = 4
  integer(kind=4),parameter:: long_kind      = 8
  integer(kind=4),parameter:: log_kind       = 4
  integer(kind=4),parameter:: real_kind      = 8
  integer(kind=4),parameter:: longdouble_kind= 8
  integer,parameter        :: MAX_STRING_LEN = 80
  integer,parameter        :: MAX_FILE_LEN   = 240

  ! Physical Constants
  !---------------------
  real(kind=longdouble_kind),parameter:: QQ_PI = 3.141592653589793238462643383279_longdouble_kind
  real(kind=real_kind)      ,parameter:: DD_PI = 3.141592653589793238462643383279_real_kind
  real(kind=real_kind)      ,parameter:: rearth       = 6.376D6    ! m
  real(kind=real_kind)      ,parameter:: g            = 9.80616D0  !  m s^-2
  real(kind=real_kind)      ,parameter:: omega        = 7.292D-5   !  radians/s
  real(kind=real_kind)      ,parameter:: Rgas         = 287.04D0
  real(kind=real_kind)      ,parameter:: Cp           = 1005.0D0
  real(kind=real_kind)      ,parameter:: p0           = 100000.0D0 !  surface pressure (mbar)
  real(kind=real_kind)      ,parameter:: MWDAIR       = 28.966D0
  real(kind=real_kind)      ,parameter:: Rwater_vapor = 461.50D0
  real(kind=real_kind)      ,parameter:: Cpwater_vapor= 1870.0D0   !ASC ANDII VERIFY PLS
  real(kind=real_kind)      ,parameter:: kappa        = Rgas/Cp
  real(kind=real_kind)      ,parameter:: Rd_on_Rv     = Rgas/Rwater_vapor  
  real(kind=real_kind)      ,parameter:: Cpd_on_Cpv   = Cp/Cpwater_vapor
  real(kind=real_kind)      ,parameter:: rrearth      = 1.0_real_kind/rearth 
  real(kind=real_kind)      ,parameter:: Lc           = 2.5D+6     ! multicloud J/Kg

  ! Dimensions 
  !------------------------------------
  integer,parameter:: np     = 4         ! NP
  integer,parameter:: nc     = 4         ! NC
  integer,parameter:: nip    = 3         ! number of interpolation values, works only for this
  integer,parameter:: nipm   = nip-1
  integer,parameter:: nep    = nipm*nc+1 ! number of points in an element  
  integer,parameter:: npdg   = 0         ! dg degree for hybrid cg/dg element  0=disabled 
  integer,parameter:: pcnst  = 25        ! number of advected constituents (including water vapor)
  integer,parameter:: qsize_d= pcnst
  integer,parameter:: ntrac_d= pcnst
  integer,parameter:: npsq   = np*np
  integer,parameter:: nlev   = 30
  integer,parameter:: nlevp  = nlev+1

end module SE_Constants

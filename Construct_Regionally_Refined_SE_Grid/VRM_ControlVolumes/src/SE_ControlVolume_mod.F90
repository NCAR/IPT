module SE_ControlVolume_mod
!===================================================================================
! 
! Purpose: Provide a datastructure (SEctrlvol_t) containing cubed-sphere Spectral 
!          Element (SE) control volumes data and a routines to initialize and save
!          its values for the given options.
!
! Author: Patrick Callaghan
!
! Description:
!     This module is an interface between the base grid information contained in 
!     the element_t datastructure from the 'cesm1_3_beta06' version of the CESM 
!     and the modules used to compute the corresponding control volumes. 
!     
!     The control volumes for the grid are computed using surfaces_mod.F90 
!     routines from HOMME.
!     The 'physics' control volume information is computed using 
!     fvm_control_volume_mod.F90 routines from CESM.
!     
!     Once computed, the control volume information can be written out to files
!     for later use in a variety of optional formats.
!
! VERSION: 1.0 
!
! Revisions: 
!   Jul 2018 - Original Version
!
! Usage: 
!    Calc_ControlVolumes(elem,SEopt,SEcv)
!            Given the initialized elem() datastructure and the options 
!            specifying resolution and related parameters (SEopt), return 
!            the iniitalized control volumes datastructure(SEcv). 
!
!    Create_SCRIP_file(SEcv,FileName)
!            Given the initialized control volume datastructure(SEcv), create a
!            SCRIP grid file with the given FileName.
!
!    Create_LATLON_file(SEcv,FileName)
!            Given the initialized control volume datastructure(SEcv), create a
!            LATLON grid file with the given FileName containing the connectivity
!            information for the unique gridpoints.
!
!    Create_PHYS_file(SEcv,FileName)
!            Given the initialized control volume datastructure(SEcv), create a
!            PHYS grid file with the given FileName containing control volumes 
!            for the physics grid.
!
!    Create_GRID_file(SEcv,FileName)
!            Given the initialized control volume datastructure(SEcv), create a
!            grid file with the given FileName containing all the SCRIP,LATLON, 
!            and PHYS grid information.
!
!===================================================================================
  ! Useful modules
  !----------------
  use netcdf
  use err_exit    ,only: iulog,endrun
  use element_mod ,only: element_t
  use SE_Options  ,only: SEoptions_t
  use SE_Constants,only: real_kind,log_kind,int_kind,np,nc,DD_PI,MAX_FILE_LEN,MAX_STRING_LEN

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private
  save

  public :: SEctrlvol_t

  public :: Calc_ControlVolumes
  public :: Create_SCRIP_file
  public :: Create_LATLON_file
  public :: Create_PHYS_file
  public :: Create_GRID_file

  private:: SEinitialized

  ! Type Definitions
  !-------------------
  type SEctrlvol_t
    character(len=MAX_STRING_LEN):: GridName
    integer                      :: global_ncol
    integer                      :: global_nc
    integer                      :: global_nsub
    integer                      :: maxvert
    integer(int_kind ),pointer   :: subelement_corners(:,:) => null()
    real   (real_kind),pointer   :: area(:)                 => null()
    real   (real_kind),pointer   :: lonp(:)                 => null()
    real   (real_kind),pointer   :: latp(:)                 => null()
    real   (real_kind),pointer   :: cv_lon(:,:)             => null()
    real   (real_kind),pointer   :: cv_lat(:,:)             => null()
!+++ARH
    real   (real_kind),pointer   :: rrfac(:)                => null()
!---ARH
  end type SEctrlvol_t

  ! Global Data
  !-------------
  logical:: SEinitialized = .false.

contains

  !=================================================================
  subroutine Calc_ControlVolumes(elem,I_SEopt,O_SEcv)
    !
    ! Calc_ControlVolumes: Interface for the calculation of control volume
    !                      information. 
    ! 
    !                      Store values into SEelem%cv datastructure for
    !                      future use.
    !==========================================================
    use element_mod ,only: element_t
    use SE_Options  ,only: SEoptions_t
    use surfaces_mod,only: cvlist,InitControlVolumesData,InitControlVolumes
    use dof_mod     ,only: CreateMetaData,UniqueCoords,UniquePoints
    integer,parameter:: NPTSMAX = 20
    ! 
    ! Passed Variables
    !---------------------
    type(element_t)  ,intent(inout):: elem(:)
    type(SEoptions_t),intent(in   ):: I_SEopt
    type(SEctrlvol_t),intent(out  ):: O_SEcv
    ! 
    ! Local Values
    !--------------
    character(len=MAX_STRING_LEN)   :: Name_Template
    integer,allocatable             :: subelement_corners(:,:)
    real(kind=real_kind),allocatable:: lonp(:)
    real(kind=real_kind),allocatable:: latp(:)
    real(kind=real_kind),allocatable::  lon_tmp(:,:,:)
    real(kind=real_kind),allocatable::  lat_tmp(:,:,:)
    real(kind=real_kind),allocatable:: area_tmp(:,:)
    real(kind=real_kind),allocatable:: lon3d (:,:)
    real(kind=real_kind),allocatable:: lat3d (:,:)
    real(kind=real_kind),allocatable:: area2d(:)
    integer                         :: global_nc,global_nsub,global_ncol
    integer                         :: ie,nxyp,st,en,maxvert,kk
!+++ARH
    integer                         :: ii,jj
    real(kind=real_kind),allocatable:: tmp_elemarea1D(:)
    real(kind=real_kind)            :: maxelemarea
!---ARH

    ! Sum up the # of unique gridpoints
    !-----------------------------------
    nxyp = 0
    do ie = 1,I_SEopt%nelem
      nxyp = nxyp + elem(ie)%idxP%NumUniquePts
    end do
    print *,' nxyp=',nxyp,' ncol=',I_SEopt%GlobalUniqueCols
    if(nxyp.ne.I_SEopt%GlobalUniqueCols) then
      call endrun('nxyp NE ncol')
    endif

    ! Set the number of physics points and subelements
    !---------------------------------------------------
    global_ncol = I_SEopt%GlobalUniqueCols
    global_nc   =  nc   * nc   *I_SEopt%nelem  ! total number of physics points
    global_nsub = (np-1)*(np-1)*I_SEopt%nelem  ! total number of subelements

!+++ARH
    ! Get element area (before boundary exchange)
    !-------------------------------------
    do ie=1,I_SEopt%nelem
      do ii = 1,np
        do jj = 1,np
          elem(ie)%tmp_elemarea(ii,jj) = sum(elem(ie)%spheremp(:,:))
        end do
      end do
    end do
!---ARH

    ! Compute Control Volumes
    !--------------------------
    if(.not. allocated(cvlist)) then
      print *,'Initializing for control volumes:'
      call InitControlVolumesData(I_SEopt)
!+++ARH, modified subroutine to compute boundary exchange for element area
      print *,'Computing control volumes:'
      call InitControlVolumes    (elem,1,I_SEopt%nelem)
!---ARH
    endif

!+++ARH
    ! Get refinement factor
    ! normalized by max element area ... should instead be average over coarse region
    !-------------------------------------
    maxelemarea = 0.d0
    do ie=1,I_SEopt%nelem
      maxelemarea = MAX(MAXVAL(elem(ie)%tmp_elemarea(:,:)),maxelemarea)
    end do

    do ie=1,I_SEopt%nelem
      elem(ie)%tmp_elemarea(:,:) = SQRT(maxelemarea/elem(ie)%tmp_elemarea(:,:))
    end do
!---ARH

    ! Allocate some workspace
    !-------------------------
    allocate(subelement_corners(global_nsub,NPTSMAX))
    allocate(lonp(nxyp))
    allocate(latp(nxyp))
    allocate( lon_tmp(np,np,NPTSMAX))
    allocate( lat_tmp(np,np,NPTSMAX))
    allocate(area_tmp(np,np        ))
    allocate( lon3d(nxyp,NPTSMAX))
    allocate( lat3d(nxyp,NPTSMAX))
    allocate(area2d(nxyp        ))
!+++ARH
    allocate(tmp_elemarea1D(nxyp))
!---ARH
    print *,'Space allocated'

    ! Calc subelement corners
    !-------------------------
    subelement_corners(:,:) = 0.d0
    call CreateMetaData(elem,subelement_corners(:,1:4))
    print *,'corners calculated'

    ! Get the lat/lon coordinates
    !-----------------------------
    st=1
    do ie=1,I_SEopt%nelem
      en = st + elem(ie)%idxP%NumUniquePts-1
      call UniqueCoords(elem(ie)%idxP,elem(ie)%spherep,latp(st:en),lonp(st:en))
      st=en+1
    end do
    print *,'latp/lonp done'

    ! Get the area values
    !---------------------
    st=1
    do ie=1,I_SEopt%nelem
      en = st + elem(ie)%idxP%NumUniquePts-1
      area_tmp(:,:) = 1.d0/elem(ie)%rspheremp(:,:)
      call UniquePoints(elem(ie)%idxP,area_tmp,area2d(st:en))
      st=en+1
    end do

    ! Get cv_lon/cv_lat values
    !--------------------------
    maxvert = 0
    do ie=1,I_SEopt%nelem
      maxvert = MAX(maxvert,MAXVAL(cvlist(ie)%nvert))
    end do
    if(NPTSMAX < maxvert) then
      call endrun('cv output requires NPTSMAX >= max number of vertex')
    endif
    print *,'maxvert=',maxvert,' NPTSMAX=',NPTSMAX

    st=1
    do ie=1,I_SEopt%nelem
      en = st + elem(ie)%idxP%NumUniquePts-1
      lon_tmp(:,:,:) = 0.d0
      lat_tmp(:,:,:) = 0.d0
      do kk=1,maxvert
        lon_tmp(:,:,kk) = cvlist(ie)%vert_latlon(kk,:,:)%lon
        lat_tmp(:,:,kk) = cvlist(ie)%vert_latlon(kk,:,:)%lat
      end do
      call UniquePoints(elem(ie)%idxP,NPTSMAX,lon_tmp,lon3d(st:en,:))
      call UniquePoints(elem(ie)%idxP,NPTSMAX,lat_tmp,lat3d(st:en,:))
      st=en+1
    end do
    print *,'cv_lat/cv_lon done'

!+++ARH
    ! Reshape rrfac (i.e., tmp_elemearea)
    !--------------------
    st=1
    do ie=1,I_SEopt%nelem
      en = st + elem(ie)%idxP%NumUniquePts-1
      call UniquePoints(elem(ie)%idxP,elem(ie)%tmp_elemarea,tmp_elemarea1D(st:en))
      st=en+1
    end do
!---ARH

    ! Get faceno values??????
    !--------------------


    ! Compute FVM Control Volumes
    !------------------------------

    ! Get phys_lat/phys_lon values
    !------------------------------

    ! Get phys_area values
    !----------------------

    ! Get phys_cv_lon/phys_cv_lat values
    !-------------------------------------

    ! Store values for later use
    !----------------------------
    O_SEcv%global_ncol = global_ncol
    O_SEcv%global_nc   = global_nc
    O_SEcv%global_nsub = global_nsub
    O_SEcv%maxvert     = maxvert
    allocate(O_SEcv%subelement_corners(global_nsub,maxvert))
    allocate(O_SEcv%area              (global_ncol)        )
    allocate(O_SEcv%lonp              (global_ncol)        )
    allocate(O_SEcv%latp              (global_ncol)        )
    allocate(O_SEcv%cv_lon            (global_ncol,maxvert))
    allocate(O_SEcv%cv_lat            (global_ncol,maxvert))

    O_SEcv%subelement_corners(:,1:maxvert) = subelement_corners(:,1:maxvert)
    O_SEcv%area              (:)           = area2d            (:)
    O_SEcv%lonp              (:)           = lonp              (:)
    O_SEcv%latp              (:)           = latp              (:)
    O_SEcv%cv_lon            (:,1:maxvert) = lon3d             (:,1:maxvert)
    O_SEcv%cv_lat            (:,1:maxvert) = lat3d             (:,1:maxvert)

!+++ARH
    allocate(O_SEcv%rrfac(global_ncol))
    O_SEcv%rrfac(:) = tmp_elemarea1D(:)
!---ARH

    if(I_SEopt%ne.eq.0) then
      O_SEcv%GridName = "Variable Resolution: "//trim(I_SEopt%MeshFile)
    elseif(I_SEopt%ne.lt.10) then
      Name_Template = "Uniform Resolution: neXnpY"
      write(Name_Template(23:23),'(i1)') I_SEopt%ne
      write(Name_Template(26:26),'(i1)') np
      O_SEcv%GridName = trim(Name_Template)
    elseif(I_SEopt%ne.lt.100) then
      Name_Template = "Uniform Resolution: neXXnpY"
      write(Name_Template(23:24),'(i2)') I_SEopt%ne
      write(Name_Template(27:27),'(i1)') np
      O_SEcv%GridName = trim(Name_Template)
    else
      Name_Template = "Uniform Resolution: neXXXnpY"
      write(Name_Template(23:25),'(i3)') I_SEopt%ne
      write(Name_Template(28:28),'(i1)') np
      O_SEcv%GridName = trim(Name_Template)
    endif

    SEinitialized = .true.

    ! My Momma always said... Clean up after yourself!
    !-------------------------------------------------
    deallocate(subelement_corners)
    deallocate(lonp)
    deallocate(latp)
    deallocate( lon_tmp)
    deallocate( lat_tmp)
    deallocate(area_tmp)
    deallocate( lon3d)
    deallocate( lat3d)
    deallocate(area2d)
!+++ARH
    deallocate(tmp_elemarea1D)
!---ARH

    ! End Routine
    !--------------
    return
  end subroutine Calc_ControlVolumes
  !==================================================================


  !=================================================================
  subroutine Create_SCRIP_file(I_SEcv,I_FileName)
    !
    ! Create_SCRIP_file: 
    !
    !==========================================================
    ! 
    ! Passed Variables
    !---------------------
    type(SEctrlvol_t)          ,intent(in):: I_SEcv
    character(len=MAX_FILE_LEN),intent(in):: I_FileName
    ! 
    ! Local Values
    !--------------
    integer             ,allocatable:: grid_dims      (:)
    integer             ,allocatable:: grid_imask     (:)
    real(kind=real_kind),allocatable:: grid_area      (:)
    real(kind=real_kind),allocatable:: grid_center_lon(:)
    real(kind=real_kind),allocatable:: grid_center_lat(:)
    real(kind=real_kind),allocatable:: grid_corner_lon(:,:)
    real(kind=real_kind),allocatable:: grid_corner_lat(:,:)
!+++ARH
    real(kind=real_kind),allocatable:: rrfac(:)
!---ARH

    real(kind=real_kind):: RTD,EPS
    integer             :: grid_size, grid_corners, grid_rank
    integer             :: ncid,retval,nn,ii,num_lat,num_lon
    integer             :: I_grid_size ,I_grid_corners ,I_grid_rank
    integer             :: I_grid_dims ,I_grid_imask   ,I_grid_area 
    integer             :: I_grid_center_lat ,I_grid_center_lon 
    integer             :: I_grid_corner_lat ,I_grid_corner_lon
!+++ARH
    integer             :: I_rrfac
!---ARH

    RTD = 180.d0/DD_PI

    ! Verify the control volume datastructure has been initialized
    !-------------------------------------------------------------
    if(.not.SEinitialized)then
      call endrun('Create_SCRIP_file() ERROR: SEcv is not initailized')
    endif

    ! set dimensions for SCRIP file
    !------------------------------
    grid_size    = I_SEcv%global_ncol
    grid_corners = I_SEcv%maxvert
    grid_rank    = 1

    ! Allocate space for the output values
    !---------------------------------------
    allocate(grid_dims      (grid_rank)             )
    allocate(grid_imask     (grid_size)             )
    allocate(grid_area      (grid_size)             )
    allocate(grid_center_lon(grid_size)             )
    allocate(grid_center_lat(grid_size)             )
    allocate(grid_corner_lon(grid_corners,grid_size))
    allocate(grid_corner_lat(grid_corners,grid_size))
!+++ARH
    allocate(rrfac(grid_size)                       )
!+++ARH
   
    ! load SCRIP output values
    !--------------------------
    grid_dims (:) = 1
    grid_imask(:) = 1

    grid_area      (:) = I_SEcv%area(:)
    grid_center_lon(:) = I_SEcv%lonp(:)*RTD
    grid_center_lat(:) = I_SEcv%latp(:)*RTD
!+++ARH
    rrfac          (:) = I_SEcv%rrfac(:)
!---ARH

    do nn=1,grid_corners
      grid_corner_lon(nn,:) = I_SEcv%cv_lon(:,nn)*RTD
      grid_corner_lat(nn,:) = I_SEcv%cv_lat(:,nn)*RTD
    end do

    !============================================
    ! Apply the 'fix' from convert_scrip_999.ncl
    !============================================
    ! Loop over SE gridpoints
    !--------------------------
    EPS     = 900.d0
    num_lat = 0
    num_lon = 0
    do ii=1,grid_size
      if(grid_corner_lat(grid_corners,ii).gt.EPS) then
        ! if the latitude of the last vertex contains a fill value, then
        ! loop thru and set all latitudes of fill-value verticies equal 
        ! to the latidude of the last valid vertex
        !----------------------------------------------------------------
        do nn=2,grid_corners
          if(grid_corner_lat(nn,ii).gt.EPS) then
            grid_corner_lat(nn,ii) = grid_corner_lat(nn-1,ii)
            num_lat = num_lat + 1
          endif
        end do
      endif
      if(grid_corner_lon(grid_corners,ii).gt.EPS) then
        ! if the longitude of the last vertex contains a fill value, then
        ! loop thru and set all longitudes of fill-value verticies equal 
        ! to the longitude of the last valid vertex
        !----------------------------------------------------------------
        do nn=2,grid_corners
          if(grid_corner_lon(nn,ii).gt.EPS) then
            grid_corner_lon(nn,ii) = grid_corner_lon(nn-1,ii)
            num_lon = num_lon + 1
          endif
        end do
      endif
    end do
    print *,'999-FIX: Number of verticies fixed: num_lat=',num_lat,' num_lon=',num_lon

    ! Try to create the file
    !-----------------------
    print *,' '
    print *,' Creating SCRIP File=',trim(I_FileName)
    print *,' '
    retval = nf90_create(trim(I_FileName),ior(NF90_CLOBBER,NF90_NETCDF4),ncid)
    if(retval.ne.0) then
      print *,' retval=',retval
      call endrun('ERROR trying to create SCRIP file')
    endif

    ! Define the dimensions
    !-------------------------
    retval = nf90_def_dim(ncid,'grid_size'   ,grid_size   ,I_grid_size   )
    retval = nf90_def_dim(ncid,'grid_corners',grid_corners,I_grid_corners)
    retval = nf90_def_dim(ncid,'grid_rank'   ,grid_rank   ,I_grid_rank   )

    ! Define the Variables
    !-----------------------
    retval = nf90_def_var(ncid,'grid_dims'      ,NF90_INT   ,(/I_grid_rank/),I_grid_dims      )
    retval = nf90_def_var(ncid,'grid_imask'     ,NF90_INT   ,(/I_grid_size/),I_grid_imask     )
    retval = nf90_def_var(ncid,'grid_area'      ,NF90_DOUBLE,(/I_grid_size/),I_grid_area      )
    retval = nf90_def_var(ncid,'grid_center_lat',NF90_DOUBLE,(/I_grid_size/),I_grid_center_lat)
    retval = nf90_def_var(ncid,'grid_center_lon',NF90_DOUBLE,(/I_grid_size/),I_grid_center_lon)
    retval = nf90_def_var(ncid,'grid_corner_lat',NF90_DOUBLE,(/I_grid_corners,I_grid_size/),I_grid_corner_lat)
    retval = nf90_def_var(ncid,'grid_corner_lon',NF90_DOUBLE,(/I_grid_corners,I_grid_size/),I_grid_corner_lon)
!+++ARH
    retval = nf90_def_var(ncid,'rrfac',NF90_DOUBLE,(/I_grid_size/),I_rrfac)
!---ARH

    retval = nf90_put_att(ncid,I_grid_area      ,"long_name","area weights")
    retval = nf90_put_att(ncid,I_grid_area      ,"units"    ,"radians^2")
    retval = nf90_put_att(ncid,I_grid_center_lat,"units"    ,"degrees")
    retval = nf90_put_att(ncid,I_grid_center_lon,"units"    ,"degrees")
    retval = nf90_put_att(ncid,I_grid_corner_lat,"units"    ,"degrees")
    retval = nf90_put_att(ncid,I_grid_corner_lon,"units"    ,"degrees")
!+++ARH
    retval = nf90_put_att(ncid,I_rrfac          ,"units"    ,"neXX/ne30")
!---ARH

    ! Global Attributes
    !--------------------
    retval = nf90_put_att(ncid,NF90_GLOBAL,"Grid",trim(I_SEcv%GridName))
    retval = nf90_put_att(ncid,NF90_GLOBAL,"Created by",'Gen_ControlVolumes.exe')
    retval = nf90_enddef(ncid)

    ! Write out values
    !------------------
    retval = nf90_put_var(ncid,I_grid_dims      ,grid_dims      )
    retval = nf90_put_var(ncid,I_grid_imask     ,grid_imask     )
    retval = nf90_put_var(ncid,I_grid_area      ,grid_area      )
    retval = nf90_put_var(ncid,I_grid_center_lat,grid_center_lat)
    retval = nf90_put_var(ncid,I_grid_center_lon,grid_center_lon)
    retval = nf90_put_var(ncid,I_grid_corner_lat,grid_corner_lat)
    retval = nf90_put_var(ncid,I_grid_corner_lon,grid_corner_lon)
!+++ARH
    retval = nf90_put_var(ncid,I_rrfac          ,rrfac          )
!---ARH
    retval = nf90_close(ncid)

    ! My Momma always said... Clean up after yourself!
    !-------------------------------------------------
    deallocate(grid_dims      )
    deallocate(grid_imask     )
    deallocate(grid_area      )
    deallocate(grid_center_lon)
    deallocate(grid_center_lat)
    deallocate(grid_corner_lon)
    deallocate(grid_corner_lat)
!+++ARH
    deallocate(rrfac)
!---ARH

    ! End Routine
    !--------------
    return
  end subroutine Create_SCRIP_file
  !==================================================================


  !=================================================================
  subroutine Create_LATLON_file(I_SEcv,I_FileName)
    !
    !  Create_LATLON_file: 
    !
    !==========================================================
    ! 
    ! Passed Variables
    !---------------------
    type(SEctrlvol_t)          ,intent(in):: I_SEcv
    character(len=MAX_FILE_LEN),intent(in):: I_FileName
    ! 
    ! Local Values
    !--------------
    real(kind=real_kind),allocatable:: area           (:)
    real(kind=real_kind),allocatable:: lat            (:)
    real(kind=real_kind),allocatable:: lon            (:)
    integer             ,allocatable:: element_corners(:,:)

    real(kind=real_kind):: RTD
    integer             :: ncol,ncorners,ncenters
    integer             :: ncid,retval,nn
    integer             :: I_ncol, I_ncorners,I_ncenters
    integer             :: I_lat , I_lon     ,I_area    ,I_element_corners

    RTD = 180.d0/DD_PI

    ! Verify the control volume datastructure has been initialized
    !-------------------------------------------------------------
    if(.not.SEinitialized)then
      call endrun('Create_LATLON_file() ERROR: SEcv is not initailized')
    endif

    ! set dimensions for the LATLON file
    !------------------------------------
    ncol     = I_SEcv%global_ncol
    ncorners = 4
    ncenters = I_SEcv%global_nsub

    ! Allocate space for the output values
    !---------------------------------------
    allocate(area           (ncol)             )
    allocate(lat            (ncol)             )
    allocate(lon            (ncol)             )
    allocate(element_corners(ncenters,ncorners))

    ! load SCRIP output values
    !--------------------------
    area(:) = I_SEcv%area(:)
    lon (:) = I_SEcv%lonp(:)*RTD
    lat (:) = I_SEcv%latp(:)*RTD
    do nn=1,ncenters
      element_corners(nn,1:4) = I_SEcv%subelement_corners(nn,1:4)
    end do

    ! Try to create the file
    !-----------------------
    print *,' '
    print *,' Creating LATLON File=',trim(I_FileName)
    print *,' '
    retval = nf90_create(trim(I_FileName),ior(NF90_CLOBBER,NF90_NETCDF4),ncid)
    if(retval.ne.0) then
      print *,' retval=',retval
      call endrun('ERROR trying to create LATLON file')
    endif

    ! Define the dimensions
    !-------------------------
    retval = nf90_def_dim(ncid,'ncol'    ,ncol    ,I_ncol    )
    retval = nf90_def_dim(ncid,'ncorners',ncorners,I_ncorners)
    retval = nf90_def_dim(ncid,'ncenters',ncenters,I_ncenters)

    ! Define the variables
    !----------------------
    retval = nf90_def_var(ncid,'lat'            ,NF90_DOUBLE,(/I_ncol/),I_lat     )
    retval = nf90_def_var(ncid,'lon'            ,NF90_DOUBLE,(/I_ncol/),I_lon     )
    retval = nf90_def_var(ncid,'area'           ,NF90_DOUBLE,(/I_ncol/),I_area    )
    retval = nf90_def_var(ncid,'element_corners',NF90_INT   ,(/I_ncenters,I_ncorners/),I_element_corners)

    retval = nf90_put_att(ncid,I_area,"long_name","area weights"    )
    retval = nf90_put_att(ncid,I_area,"units"    ,"radians^2"       )
    retval = nf90_put_att(ncid,I_lat ,"long_name","column latitude" )
    retval = nf90_put_att(ncid,I_lat ,"units"    ,"degrees_north"   )
    retval = nf90_put_att(ncid,I_lon ,"long_name","column longitude")
    retval = nf90_put_att(ncid,I_lon ,"units"    ,"degrees_east"    )

    ! Global Attributes
    !--------------------
    retval = nf90_put_att(ncid,NF90_GLOBAL,"Grid",trim(I_SEcv%GridName))
    retval = nf90_put_att(ncid,NF90_GLOBAL,"Created by",'Gen_ControlVolumes.exe')
    retval = nf90_enddef(ncid)

    ! Write out values
    !------------------
    retval = nf90_put_var(ncid,I_lat            ,lat            )
    retval = nf90_put_var(ncid,I_lon            ,lon            )
    retval = nf90_put_var(ncid,I_area           ,area           )
    retval = nf90_put_var(ncid,I_element_corners,element_corners)
    retval = nf90_close(ncid)

    ! My Momma always said... Clean up after yourself!
    !-------------------------------------------------
    deallocate(area           )
    deallocate(lat            )
    deallocate(lon            )
    deallocate(element_corners)

    ! End Routine
    !--------------
    return
  end subroutine Create_LATLON_file
  !==================================================================


  !=================================================================
  subroutine Create_PHYS_file(I_SEcv,I_FileName)
    !
    ! Create_PHYS_file:
    !
    !==========================================================
    ! 
    ! Passed Variables
    !---------------------
    type(SEctrlvol_t)          ,intent(in):: I_SEcv
    character(len=MAX_FILE_LEN),intent(in):: I_FileName
    ! 
    ! Local Values
    !--------------

    ! Verify the control volume datastructure has been initialized
    !-------------------------------------------------------------
    if(.not.SEinitialized)then
      call endrun('Create_PHYS_file() ERROR: SEcv is not initailized')
    endif

    print *,' '
    print *,' Create_PHYS_file() not implemented yet' 
    print *,' '

    ! End Routine
    !--------------
    return
  end subroutine Create_PHYS_file
  !==================================================================


  !=================================================================
  subroutine Create_GRID_file(I_SEcv,I_FileName)
    !
    ! Create_GRID_file: There is no name conflict, so why not
    !                   combine all of the LATLON, SCRIP, and PHYS
    !                   grid info into one file???
    !==========================================================
    ! 
    ! Passed Variables
    !---------------------
    type(SEctrlvol_t)          ,intent(in):: I_SEcv
    character(len=MAX_FILE_LEN),intent(in):: I_FileName
    ! 
    ! Local Values
    !--------------

    ! Verify the control volume datastructure has been initialized
    !-------------------------------------------------------------
    if(.not.SEinitialized)then
      call endrun('Create_GRID_file() ERROR: SEcv is not initailized')
    endif

    print *,' '
    print *,' Create_GRID_file() not implemented yet' 
    print *,' '

    ! End Routine
    !--------------
    return
  end subroutine Create_GRID_file
  !==================================================================

end module SE_ControlVolume_mod

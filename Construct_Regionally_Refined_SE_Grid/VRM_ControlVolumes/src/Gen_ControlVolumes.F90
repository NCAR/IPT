  program Gen_ControlVolumes
    ! Useful Modules
    !------------------
    use SE_Constants        ,only: MAX_FILE_LEN
    use SE_Options          ,only: SEoptions_t, Init_SEoptions
    use SE_Element_mod      ,only: SEelem_t,Init_SEelem
    use SE_ControlVolume_mod,only: SEctrlvol_t,Create_SCRIP_file,       &
                                   Create_LATLON_file,Create_PHYS_file, &
                                   Create_GRID_file
    implicit none

    ! Variables
    !-----------
    type(SEoptions_t)                   :: SEopt
    type(SEelem_t   ),target,allocatable:: SEelem(:)
    type(SEctrlvol_t)                   :: SEcv
    character(len=MAX_FILE_LEN)         :: SCRIP_file
    character(len=MAX_FILE_LEN)         :: LATLON_file
    character(len=MAX_FILE_LEN)         :: PHYS_file
    character(len=MAX_FILE_LEN)         :: GRID_file
    character(len=MAX_FILE_LEN)         :: NL_file
    integer          :: Narg
    character(len=32):: Arg

    ! Get namelist file from command line arguments
    !-----------------------------------------------
    Narg = command_argument_count()
    if(Narg.eq.1) then
      call get_command_argument(1,Arg)
      NL_file = trim(Arg)
    else
      print *,' '
      print *,'USAGE: Gen_ControlVolumes.exe   NameListFile'
      print *,' Please provide a namelist file'
      print *,' '
      stop
    endif
    print *,' '
    print *,'Gen_ControlVolumes.exe Reading namlist values from: ',NL_file
    print *,' '

    ! Init SE modules
    !----------------------
    call Init_SEoptions(SEopt,NL_file)
    call Init_SEelem   (SEopt,SEelem,SEcv)

    ! Write desired grid files with control volume values
    !------------------------------------------------------
    if(SEopt%create_SCRIP_file) then
      SCRIP_file = trim(SEopt%OutputPath)//trim(SEopt%SCRIP_filename)
      call Create_SCRIP_file(SEcv,SCRIP_file)
    endif

    if(SEopt%create_LATLON_file) then
      LATLON_file = trim(SEopt%OutputPath)//trim(SEopt%LATLON_filename)
      call Create_LATLON_file(SEcv,LATLON_file)
    endif

    if(SEopt%create_PHYS_file) then
      PHYS_file = trim(SEopt%OutputPath)//trim(SEopt%PHYS_filename)
      call Create_PHYS_file(SEcv,PHYS_file)
    endif

    if(SEopt%create_GRID_file) then
      GRID_file = trim(SEopt%OutputPath)//trim(SEopt%GRID_filename)
      call Create_GRID_file(SEcv,GRID_file)
    endif

  end program

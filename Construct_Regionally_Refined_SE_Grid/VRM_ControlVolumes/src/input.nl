!
!  set input_nl values:
!  ----------------------
!     GridPath = the path containing input EXODUS file and where the 
!                created SCRIP and LATLON files will be created.
!     GridName = The root name for the grid. Note that this is just 
!                the root filename, the program will append the rest 
!                of the name for the actual files.
!
!                The program will open the EXODUS file: 
!
!                          $(GRIDPATH)/$(GRIDNAME)_EXODUS.nc
!
!                 and create the files:
!
!                          $(GRIDPATH)/$(GRIDNAME)_np4_SCRIP.nc
!                          $(GRIDPATH)/$(GRIDNAME)_np4_LATLON.nc
!================================================================================

&input_nl
  GridPath = '/path/to/my/repo/grids/'
  GridName = 'GRIDNAME_neZZxRR'
/

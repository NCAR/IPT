# Makefile to build interpic on various platforms
# Note: If netcdf library is not built in the standard location, you must set the
# variables INC_NETCDF and LIB_NETCDF

RM = rm

.SUFFIXES:
.SUFFIXES: .o .f90 .F90  .f

# Check whether NetCDF library and include directories specified in environment
# or on make commandline.
ifeq ($(strip $(LIB_NETCDF)),)
  LIB_NETCDF := /usr/local/lib
endif
ifeq ($(strip $(INC_NETCDF)),)
  INC_NETCDF := /usr/local/include
endif

### INC_NETCDF := /usr/local/netcdf-pgi/include
### LIB_NETCDF := /usr/local/netcdf-pgi/lib


# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------

ifeq ($(UNAMEM),CRAY)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC = f90
FFLAGS = -64 -c -trapuv -I$(INC_NETCDF) -g -C -DEBUG:trap_uninitialized=ON
LDFLAGS = -64 -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC = f90
FFLAGS = -c -stackvar -f -I$(INC_NETCDF) -g
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC = xlf90
FFLAGS = -c -I$(INC_NETCDF) -qsuffix=cpp=F90 -WF,-DAIX -g -qfullpath
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)

# g95
#FC = g95 
#FFLAGS =  -c -I$(INC_NETCDF) -g -ftrace=full
#FC = f95 
#FFLAGS =  -c -I$(INC_NETCDF) 

# pgf90
### FC = pgf90

### FFLAGS =  -c -I$(INC_NETCDF) -g -Ktrap=fp -Mrecursive -Mbounds

# gfortran
##cub INC_NETCDF := /usr/lib64/gfortran/modules
##INC_NETCDF := /usr/lib64/gfortran/modules

FC = gfortran
#FFLAGS =  -c -I$(INC_NETCDF) -O3
FFLAGS =  -c -I$(INC_NETCDF) -g 

# lf95
#FC = lf95
#FFLAGS =  -c -I$(INC_NETCDF) -g --chk a,e,s,u --pca --trace --trap

# ifort
##FC = ifort
##FFLAGS =  -c -I$(INC_NETCDF) -O
##FFLAGS =  -c -I$(INC_NETCDF) -g -check all -fpe0 -traceback -ftz -convert big_endian -fp-model precise

###LDFLAGS = -L$(LIB_NETCDF) -lnetcdff -lnetcdf
#wisp LDFLAGS = -L$(LIB_NETCDF) -lnetcdff -lnetcdf  /lib64/libgfortran.so.3
##LDFLAGS = -L$(LIB_NETCDF) -lnetcdff -lnetcdf  /lib64/libgfortran.so.4
LDFLAGS = -L$(LIB_NETCDF) -lnetcdff -lnetcdf 
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------
EXEDIR = .

EXENAM01 = Gen_ControlVolumes.exe
OBJS01  := Gen_ControlVolumes.o \
	SE_ControlVolume_mod.o surfaces_mod.o SE_Element_mod.o mesh_mod.o \
	dof_mod.o derivative_mod.o mass_matrix_mod.o cube_mod.o element_mod.o edge_mod.o \
	spacecurve_mod.o metagraph_mod.o gridgraph_mod.o coordinate_systems_mod.o \
	quadrature_mod.o SE_Options.o SE_Constants.o err_exit.o thread_mod.o

.F90.o:
	$(FC) $(FFLAGS) $<

$(EXEDIR)/$(EXENAM01): $(OBJS01)
	$(FC) -O -o $@ $(OBJS01) $(LDFLAGS)

clean:
	$(RM) -f *.o *.mod *.exe


Gen_ControlVolumes.o: Gen_ControlVolumes.F90 \
	SE_ControlVolume_mod.o surfaces_mod.o SE_Element_mod.o mesh_mod.o \
	dof_mod.o derivative_mod.o mass_matrix_mod.o cube_mod.o element_mod.o edge_mod.o \
	spacecurve_mod.o metagraph_mod.o gridgraph_mod.o coordinate_systems_mod.o \
	quadrature_mod.o SE_Options.o SE_Constants.o err_exit.o thread_mod.o


SE_ControlVolume_mod.o: SE_ControlVolume_mod.F90 \
	surfaces_mod.o SE_Constants.o SE_Options.o element_mod.o dof_mod.o  err_exit.o
surfaces_mod.o: surfaces_mod.F90 \
	 SE_Constants.o SE_Options.o edge_mod.o element_mod.o coordinate_systems_mod.o \
	cube_mod.o quadrature_mod.o err_exit.o
SE_Element_mod.o: SE_Element_mod.F90 \
	SE_ControlVolume_mod.o surfaces_mod.o element_mod.o coordinate_systems_mod.o \
	gridgraph_mod.o derivative_mod.o \
	quadrature_mod.o SE_Constants.o SE_Options.o \
	thread_mod.o cube_mod.o spacecurve_mod.o mesh_mod.o metagraph_mod.o \
	mass_matrix_mod.o dof_mod.o
mesh_mod.o: mesh_mod.F90 \
	SE_Constants.o SE_Options.o err_exit.o coordinate_systems_mod.o \
	gridgraph_mod.o spacecurve_mod.o cube_mod.o element_mod.o 
dof_mod.o: dof_mod.F90 \
	SE_Constants.o SE_Options.o element_mod.o quadrature_mod.o edge_mod.o \
	coordinate_systems_mod.o
derivative_mod.o: derivative_mod.F90 \
	element_mod.o quadrature_mod.o SE_Constants.o SE_Options.o
mass_matrix_mod.o: mass_matrix_mod.F90 \
	SE_Constants.o SE_Options.o element_mod.o quadrature_mod.o edge_mod.o
cube_mod.o: cube_mod.F90 \
	SE_Constants.o SE_Options.o coordinate_systems_mod.o err_exit.o element_mod.o \
	gridgraph_mod.o spacecurve_mod.o
element_mod.o: element_mod.F90 \
	edge_mod.o coordinate_systems_mod.o gridgraph_mod.o SE_Constants.o SE_Options.o
edge_mod.o: edge_mod.F90 \
	SE_Constants.o SE_Options.o thread_mod.o err_exit.o coordinate_systems_mod.o
spacecurve_mod.o: spacecurve_mod.F90 \
	SE_Options.o gridgraph_mod.o err_exit.o 
metagraph_mod.o: metagraph_mod.F90 \
	SE_Constants.o gridgraph_mod.o
gridgraph_mod.o: gridgraph_mod.F90 \
	SE_Constants.o SE_Options.o err_exit.o
coordinate_systems_mod.o: coordinate_systems_mod.F90 \
	SE_Constants.o err_exit.o
quadrature_mod.o: quadrature_mod.F90 \
	SE_Constants.o
SE_Options.o: SE_Options.F90 \
	SE_Constants.o err_exit.o
SE_Constants.o: SE_Constants.F90
err_exit.o: err_exit.F90
thread_mod.o: thread_mod.F90

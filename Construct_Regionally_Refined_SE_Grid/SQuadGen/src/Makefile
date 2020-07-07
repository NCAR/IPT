# Compiler and flags
#CC= icpc
CC= g++

# NETCDF library directories
NETCDF_INCLUDEDIR=./include
NETCDF_LIBDIR=./lib
#NETCDF_INCLUDEDIR=/opt/local/include
#NETCDF_LIBDIR=/opt/local/lib
#NETCDF_INCLUDEDIR= /usr/local/netcdf-c-4.6.1-f-4.4.4-gcc-g++-gfortran-8.1.0/include
#NETCDF_LIBDIR= /usr/local/netcdf-c-4.6.1-f-4.4.4-gcc-g++-gfortran-8.1.0/lib
#NETCDF_INCLUDEDIR=/glade/u/apps/ch/opt/netcdf/4.4.1.1/intel/17.0.1/include
#NETCDF_LIBDIR=/glade/u/apps/ch/opt/netcdf/4.4.1.1/intel/17.0.1/lib

# Library files to include
LDFILES= -lnetcdf -lnetcdf_c++

# Local files
FILES= SQuadGen.cpp \
	   GridElements.cpp \
	   CSRefinementMap.cpp \
       CubedSphereGrid.cpp \
	   IcosahedralFlagGrid.cpp \
	   RefineGrid.cpp \
       RefinementTemplateCUBIT.cpp \
	   RefinementTemplateLOWCONN.cpp \
	   RefinementTemplateLOWCONNOLD.cpp \
       SpringDynamics.cpp \
       Tessellate.cpp\
       lodepng.cpp

# Load system-specific defaults
CFLAGS+= -I$(NETCDF_INCLUDEDIR)
LDFLAGS+= -L$(NETCDF_LIBDIR)

include Make.defs

##
## Build instructions
##
SQuadGen: $(FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)
	mv SQuadGen ..

##
## Clean
##
clean:
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

SQuadGen.o: /usr/include/png.h /usr/include/pnglibconf.h
SQuadGen.o: /usr/include/pngconf.h /usr/include/limits.h
SQuadGen.o: /usr/include/bits/libc-header-start.h /usr/include/features.h
SQuadGen.o: /usr/include/stdc-predef.h /usr/include/sys/cdefs.h
SQuadGen.o: /usr/include/bits/wordsize.h /usr/include/bits/long-double.h
SQuadGen.o: /usr/include/gnu/stubs.h /usr/include/bits/posix1_lim.h
SQuadGen.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
SQuadGen.o: /usr/include/bits/posix2_lim.h /usr/include/stdio.h
SQuadGen.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
SQuadGen.o: /usr/include/bits/types/__FILE.h /usr/include/bits/types/FILE.h
SQuadGen.o: /usr/include/libio.h /usr/include/_G_config.h
SQuadGen.o: /usr/include/bits/types/__mbstate_t.h
SQuadGen.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
SQuadGen.o: /usr/include/setjmp.h /usr/include/bits/setjmp.h
SQuadGen.o: /usr/include/bits/types/__sigset_t.h /usr/include/time.h
SQuadGen.o: /usr/include/bits/time.h /usr/include/bits/types/clock_t.h
SQuadGen.o: /usr/include/bits/types/time_t.h
SQuadGen.o: /usr/include/bits/types/struct_tm.h
SQuadGen.o: /usr/include/bits/types/struct_timespec.h
SQuadGen.o: /usr/include/bits/types/clockid_t.h
SQuadGen.o: /usr/include/bits/types/timer_t.h
SQuadGen.o: /usr/include/bits/types/struct_itimerspec.h
SQuadGen.o: /usr/include/bits/types/locale_t.h
SQuadGen.o: /usr/include/bits/types/__locale_t.h /usr/include/netcdfcpp.h
SQuadGen.o: /usr/include/ncvalues.h /usr/include/netcdf.h
SQuadGen.o: /usr/include/errno.h /usr/include/bits/errno.h
SQuadGen.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
SQuadGen.o: /usr/include/asm-generic/errno.h
SQuadGen.o: /usr/include/asm-generic/errno-base.h GridElements.h Exception.h
SQuadGen.o: CubedSphereGrid.h IcosahedralFlagGrid.h CSRefinementMap.h
SQuadGen.o: DataMatrix3D.h DataVector.h RefineGrid.h
SQuadGen.o: RefinementTemplateCUBIT.h RefinementTemplateLOWCONN.h
SQuadGen.o: RefinementTemplateLOWCONNOLD.h SpringDynamics.h Tessellate.h
SQuadGen.o: CommandLine.h

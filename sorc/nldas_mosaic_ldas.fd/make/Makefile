# Set up special characters


null  :=
space := $(null) $(null)
doctool :=../../utils/docsgen.sh

# Check for directory in which to put executable
ifeq ($(MODEL_EXEDIR),$(null))
MODEL_EXEDIR := .
endif

# Check for name of executable
ifeq ($(EXENAME),$(null))
EXENAME := nldas_mosaic_ldas
endif

# Check if SPMD is defined in "misc.h"
# Ensure that it is defined and not just "undef SPMD" set in file
ifeq ($(SPMD),$(null))
  SPMDSET := $(shell /bin/grep SPMD misc.h)
  ifneq (,$(findstring define,$(SPMDSET)))
    SPMD := TRUE
  else
    SPMD := FALSE
  endif
endif

# Determine platform
UNAMES := $(shell uname -s)
UMACHINE := $(shell uname -m)

ifeq ($(UNAMES),Linux)
INC_NETCDF :=$(NETCDF_FFLAGS)
LIB_NETCDF :=$(LD_LIBRARY_PATH)
# we cannot find MPI lib and include on Tide 
# (Checking with Geroge V. on 8 March 2013)
#LIB_MPI := /usr/lib
#INC_MPI := /usr/include
LIB_MPI :=
INC_MPI := 
ESMF_DIR   := /usrx/local/esmf-3.1.0rp5
LIB_ESMF   := $(ESMF_DIR)/lib/libO/Linux.intel.64.intelmpi.default
MOD_ESMF   := $(ESMF_DIR)/mod/modO/Linux.intel.64.intelmpi.default

endif

# Load dependency search path.
dirs := . $(shell cat Filepath)
# Set cpp search path, include netcdf
cpp_dirs := $(dirs) $(INC_NETCDF) $(INC_MPI)
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
VPATH    := $(subst $(space),:,$(VPATH))

#------------------------------------------------------------------------
# Primary target: build the model
#------------------------------------------------------------------------
all: $(MODEL_EXEDIR)/$(EXENAME)

# Get list of files and determine objects and dependency files
FIND_FILES = $(wildcard $(dir)/*.F $(dir)/*.f $(dir)/*.F90 $(dir)/*.c)
FILES      = $(foreach dir, $(dirs),$(FIND_FILES))
SOURCES   := $(sort $(notdir $(FILES)))
DEPS      := $(addsuffix .d, $(basename $(SOURCES)))
OBJS      := $(addsuffix .o, $(basename $(SOURCES)))
DOCS      := $(addsuffix .tex, $(basename $(SOURCES)))

$(MODEL_EXEDIR)/$(EXENAME): $(OBJS)
	 $(FC) -o $@ $(OBJS) $(FOPTS) $(LDFLAGS)
debug: $(OBJS)
        echo "FFLAGS: $(FFLAGS)"
        echo "LDFLAGS: $(LDFLAGS)"
        echo "OBJS: $(OBJS)"

#***********************************************************************
#********** Architecture-specific flags and rules***********************
#***********************************************************************
#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UMACHINE),i686)
ESMF_ARCH = linux_absoft

ifeq ($(ESMF_ARCH),linux_pgi)

FC          := $(MPI_PREFIX)/bin/mpif90
CC          := $(MPI_PREFIX)/bin/mpicc
#FC          := pgf90
CPP         := /lib/cpp

CFLAGS      := $(cpp_path) -c -O2 -DPGI 
FFLAGS      := $(cpp_path) -c -O2 -DHIDE_SHR_MSG -DNO_SHR_VMATH -DPGI -Mbyteswapio -r4 -i4 -Mpreprocess 

LDFLAGS     := -lmpich -lm -Mlfs

endif
ifeq ($(ESMF_ARCH),linux_ifc)

#FC          := $(MPI_PREFIX)/bin/mpif90
FC          := /opt/intel_fc_80/bin/ifort
CPP         := /lib/cpp

CFLAGS      := $(cpp_path) -c -O2 -DIFC

#FFLAGS      = $(cpp_path) -c -I$(MOD_ESMF)/$(ESMF_ARCH) -DHIDE_SHR_MSG -DNO_SHR_VMATH -O

FFLAGS      = $(cpp_path) -c -O -nomixed_str_len_arg -names lowercase -convert big_endian -assume byterecl -DHIDE_SHR_MSG -DNO_SHR_VMATH -DIFC

#LDFLAGS     = -L$(LIB_ESMF)/$(ESMF_ARCH) -lesmf -lmpich
LDFLAGS     = 

endif
ifeq ($(ESMF_ARCH),linux_absoft)

#FC          := $(MPI_PREFIX)/bin/mpif90
#CC          := $(MPI_PREFIX)/bin/mpicc
FC          := f90
CC          := gcc
CPP         := /lib/cpp

# non-opendap, no netcdf, no esmf
CFLAGS      := $(cpp_path) -c -O2 -Wall -DABSOFT -DLITTLE_ENDIAN
FFLAGS      := $(cpp_path) -c -O2 -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -YDEALLOC=ALL -DHIDE_SHR_MSG -DNO_SHR_VMATH -DABSOFT -DLITTLE_ENDIAN

# non-opendap, no esmf
#CFLAGS      := $(cpp_path) -c -O2 -Wall -DABSOFT -DLITTLE_ENDIAN
#FFLAGS      := $(cpp_path) -c -O2 -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -YDEALLOC=ALL -p$(INC_NETCDF) -DHIDE_SHR_MSG -DNO_SHR_VMATH -DABSOFT -DLITTLE_ENDIAN

# non-opendap
#CFLAGS      := $(cpp_path) -c -O2 -Wall -DABSOFT -DLITTLE_ENDIAN
#FFLAGS      := $(cpp_path) -c -O2 -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -YDEALLOC=ALL -p$(MOD_ESMF)/$(ESMF_ARCH) -p$(INC_NETCDF) -DHIDE_SHR_MSG -DNO_SHR_VMATH -DABSOFT -DLITTLE_ENDIAN

# opendap, no netcdf, no esmf
#CFLAGS      := $(cpp_path) -c -O2 -Wall -DABSOFT -DLITTLE_ENDIAN -DOPENDAP
#FFLAGS      := $(cpp_path) -c -O2 -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -YDEALLOC=ALL -DHIDE_SHR_MSG -DNO_SHR_VMATH -DABSOFT -DLITTLE_ENDIAN -DOPENDAP

# opendap, farmer-dog-bones, no netcdf, no esmf
#CFLAGS      := $(cpp_path) -c -O2 -Wall -DABSOFT -DLITTLE_ENDIAN -DOPENDAP -DFARMER_DOG_BONES
#FFLAGS      := $(cpp_path) -c -O2 -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -YDEALLOC=ALL -DHIDE_SHR_MSG -DNO_SHR_VMATH -DABSOFT -DLITTLE_ENDIAN -DOPENDAP -DFARMER_DOG_BONES

# opendap
#CFLAGS      := $(cpp_path) -c -O2 -Wall -DABSOFT -DLITTLE_ENDIAN -DOPENDAP
#FFLAGS      := $(cpp_path) -c -O2 -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -YDEALLOC=ALL -p$(MOD_ESMF)/$(ESMF_ARCH) -p$(INC_NETCDF) -DHIDE_SHR_MSG -DNO_SHR_VMATH -DABSOFT -DLITTLE_ENDIAN -DOPENDAP

# debugging with opendap
#CFLAGS      := $(cpp_path) -c -g -DABSOFT -DLITTLE_ENDIAN -DOPENDAP
#FFLAGS      := $(cpp_path) -c -O1 -g -Rb -Rc -Rs -Rp -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -YDEALLOC=ALL -p$(MOD_ESMF)/$(ESMF_ARCH) -p$(INC_NETCDF) -DHIDE_SHR_MSG -DNO_SHR_VMATH -DABSOFT -DLITTLE_ENDIAN -DOPENDAP

# debugging without opendap
#CFLAGS      := $(cpp_path) -c -g -DABSOFT -DLITTLE_ENDIAN
#FFLAGS      := $(cpp_path) -c -O1 -g -Rb -Rc -Rs -Rp -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -YDEALLOC=ALL -p$(MOD_ESMF)/$(ESMF_ARCH) -p$(INC_NETCDF) -DHIDE_SHR_MSG -DNO_SHR_VMATH -DABSOFT -DLITTLE_ENDIAN

# debugging without opendap, esmf, netcdf
##CFLAGS      := $(cpp_path) -c -gdwarf -DABSOFT -DLITTLE_ENDIAN
#CFLAGS      := $(cpp_path) -c -g -DABSOFT -DLITTLE_ENDIAN
#FFLAGS      := $(cpp_path) -c -O1 -g -m0 -Rb -Rc -Rs -Rp -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -YDEALLOC=ALL -DHIDE_SHR_MSG -DNO_SHR_VMATH -DABSOFT -DLITTLE_ENDIAN

# profiling -- don't forget the -P in LDFLAGS
#CFLAGS      := $(cpp_path) -c -O2 -DABSOFT -DLITTLE_ENDIAN -DOPENDAP -pg
#FFLAGS      := $(cpp_path) -c -O2 -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -YDEALLOC=ALL -p$(MOD_ESMF)/$(ESMF_ARCH) -p$(INC_NETCDF) -DHIDE_SHR_MSG -DNO_SHR_VMATH -DABSOFT -DLITTLE_ENDIAN -DOPENDAP -P

#LDFLAGS     := -L$(LIB_ESMF)/$(ESMF_ARCH) -L$(LIB_NETCDF) -lnetcdf -lesmf -lmpich -lU77 -lm
#LDFLAGS     := -L$(LIB_NETCDF) -lnetcdf -lU77 -lm
#LDFLAGS     := -lmpich -lU77 -lm
LDFLAGS     := -lU77 -lm
endif

ifeq ($(ESMF_ARCH),linux_lf95)

#FC          := $(MPI_PREFIX)/bin/mpif90
FC          := lf95
CPP         := /lib/cpp
CFLAGS      := $(cpp_path) -c -O -DUSE_GCC -DLAHEY -DLITTLE_ENDIAN
FFLAGS      := $(cpp_path) -c -O -DHIDE_SHR_MSG -DLINUX -DNO_SHR_VMATH -I$(MOD_ESMF)/$(ESMF_ARCH) -DLAHEY -DLITTLE_ENDIAN
#LDFLAGS     := -L$(LIB_ESMF)/$(ESMF_ARCH) -lesmf -L$(LIB_MPI) -lmpich -s --staticlink
#LDFLAGS     := -s --staticlink
LDFLAGS     := -s 
 
endif

# Library directories
#LIB_DIR  = ../lib/pc-32/$(ESMF_ARCH)/
LIB_DIR  = ../w3lib/
HDFLIBDIR  = $(LIB_DIR)
GFIOLIBDIR = $(LIB_DIR)
CPPFLAGS    := -P
PSASINC     := 
#FOPTS = $(LIB_DIR)bacio_32_pclinux $(LIB_DIR)w3lib_32_pclinux 
FOPTS = $(LIB_DIR)libw3.a 
# WARNING: -mp and -g together cause wrong answers

# WARNING: - Don't run hybrid on SGI (that's what the -= -mp is all about)

ifeq ($(SPMD),TRUE)
#  FFLAGS  -= -mp
#  FFLAGS  += -macro_expand
#  FFLAGS  += -I$(INC_MPI) -macro_expand

#  LDFLAGS += -L$(LIB_MPI) -lmpi
else
  FFLAGS  += -DHIDE_MPI
endif

.SUFFIXES:
.SUFFIXES: .F90 .c .o

.F90.o:
	$(FC) $(FFLAGS) $<
.c.o:
	$(CC) $(cpp_path) $(CFLAGS) $<

endif

ifeq ($(UMACHINE),ia64)
FC          := ifort
#FC          := ifort
#FC          := $(MPI_PREFIX)/bin/mpif90
CC          := /usr/bin/icc
CPP         := /lib/cpp

CFLAGS      := $(cpp_path) -c -I $(INC_MPI) -DIFC
#FFLAGS      := $(cpp_path) -c -O -DHIDE_SHR_MSG -DNO_SHR_VMATH -DIFC

#CFLAGS      := $(cpp_path) -c -DIFC -I $(INC_MPI)
FFLAGS      := $(cpp_path) -c -traceback -fpe0 -nomixed_str_len_arg -names lowercase -convert big_endian -assume byterecl -DHIDE_SHR_MSG -DNO_SHR_VMATH -DIFC

#LDFLAGS     := -Vaxlib -lmpi
LDFLAGS     := -Vaxlib

# Library directories
#LIB_DIR  = ../lib/pc-32/$(ESMF_ARCH)/
LIB_DIR  = /ldas7/wd20yx/Grib/
LIB_DIR1 = /usr/lib/
HDFLIBDIR  = $(LIB_DIR)
GFIOLIBDIR = $(LIB_DIR)
CPPFLAGS    := -P
PSASINC     := 
#FOPTS = $(LIB_DIR)bacio_32_pclinux $(LIB_DIR)w3lib_32_pclinux 
FOPTS = $(LIB_DIR1)libip_d.a $(LIB_DIR1)libbacio.a $(LIB_DIR)libw3.a
# WARNING: -mp and -g together cause wrong answers

# WARNING: - Don't run hybrid on SGI (that's what the -= -mp is all about)

ifeq ($(SPMD),TRUE)
#  FFLAGS  -= -mp
#  FFLAGS  += -macro_expand
#  FFLAGS  += -I$(INC_MPI) -macro_expand

#  LDFLAGS += -L$(LIB_MPI) -lmpi
else
  FFLAGS  += -DHIDE_MPI
endif

.SUFFIXES:
.SUFFIXES: .F90 .c .o

.F90.o:
	$(FC) $(FFLAGS) $<
.c.o:
	$(CC) $(cpp_path) $(CFLAGS) $<

endif

ifeq ($(UMACHINE),x86_64)
FC          := ifort
CC          := /usr/bin/gcc
CPP         := /lib/cpp

CFLAGS      := $(cpp_path) -c -O2 -DIFC
FFLAGS      := $(cpp_path) -c -O1 -fpe0 -nomixed_str_len_arg -names lowercase -132 -convert big_endian -assume byterecl -DHIDE_SHR_MSG -DNO_SHR_VMATH -DIFC

LDFLAGS     := 

# Library directories
LIB_DIR  = /nwprod/lib/
LIB_DIR1 = /usr/lib/
HDFLIBDIR  = $(LIB_DIR)
GFIOLIBDIR = $(LIB_DIR)
CPPFLAGS    := -P
PSASINC     := 

 FOPTS = ${G2_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4} ${IP_LIB4} ${SP_LIB4}  -openmp ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}

#FOPTS = $(LIB_DIR)libip_d.a $(LIB_DIR)libbacio_4.a $(LIB_DIR)libw3nco_4.a

# WARNING: -mp and -g together cause wrong answers

# WARNING: - Don't run hybrid on SGI (that's what the -= -mp is all about)

ifeq ($(SPMD),TRUE)
#  FFLAGS  -= -mp
#  FFLAGS  += -macro_expand
#  FFLAGS  += -I$(INC_MPI) -macro_expand

#  LDFLAGS += -L$(LIB_MPI) -lmpi
else
  FFLAGS  += -DHIDE_MPI
endif

.SUFFIXES:
.SUFFIXES: .F90 .c .f .o

.F90.o:
	$(FC) $(FFLAGS) $<
.f.o:
	$(FC) $(FFLAGS) $<
.c.o:
	$(CC) $(cpp_path) $(CFLAGS) $<

endif


RM := rm
# Add user defined compiler flags if set, and replace FC if USER option set.
FFLAGS  += $(USER_FFLAGS)
ifneq ($(USER_FC),$(null))
FC := $(USER_FC)
endif

clean:
	$(RM) -f *.o *.d *.mod *.stb  $(MODEL_EXEDIR)/$(EXENAME)

realclean:
	$(RM) -f *.o *.d *.mod *.stb  $(MODEL_EXEDIR)/$(EXENAME)
doc:
	$(doctool) 
#------------------------------------------------------------------------
#!!!!!!!!!!!!!!!!DO NOT EDIT BELOW THIS LINE.!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#------------------------------------------------------------------------
# These rules cause a dependency file to be generated for each source
# file.  It is assumed that the tool "makdep" (provided with this
# distribution in clm2/tools/makdep) has been built and is available in
# the user's $PATH.  Files contained in the clm2 distribution are the
# only files which are considered in generating each dependency.  The
# following filters are applied to exclude any files which are not in
# the distribution (e.g. system header files like stdio.h).
#
#  1) Remove full paths from dependencies. This means gnumake will not break
#     if new versions of files are created in the directory hierarchy
#     specified by VPATH.
#
#  2) Because of 1) above, remove any file dependencies for files not in the
#     clm2 source distribution.
#
# Finally, add the dependency file as a target of the dependency rules.  This
# is done so that the dependency file will automatically be regenerated
# when necessary.
#
#     i.e. change rule
#       make.o : make.c make.h
#       to:
#       make.o make.d : make.c make.h
#------------------------------------------------------------------------
DEPGEN := ./MAKDEP/makdep -s F
%.d : %.c
	@echo "Building dependency file $@"
	@$(DEPGEN) -f $(cpp_path) $< > $@
%.d : %.f
	@echo "Building dependency file $@"
	@$(DEPGEN) -f $(cpp_path) $< > $@
%.d : %.F90
	@echo "Building dependency file $@"
	@$(DEPGEN) -f $(cpp_path) $<  > $@
%.d : %.F
	@echo "Building dependency file $@"
	@$(DEPGEN) -f $(cpp_path) $< > $@
#
# if goal is clean or realclean then don't include .d files
# without this is a hack, missing dependency files will be created
# and then deleted as part of the cleaning process
#
INCLUDE_DEPS=TRUE
ifeq ($(MAKECMDGOALS), realclean)
 INCLUDE_DEPS=FALSE
endif
ifeq ($(MAKECMDGOALS), clean)
 INCLUDE_DEPS=FALSE
endif

ifeq ($(INCLUDE_DEPS), TRUE)
-include $(DEPS)
endif

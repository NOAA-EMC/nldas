# Filename=  $GRIB_ENV/config.os	12/20/97 Nakajima, SAIC MRY
# Include  this at beginning of all GRIB library related Makefiles.
###################################################################
SHELL		= /bin/sh

#.......................................................
# General Variables
#       LIBRARY 	name of the GRIB library (leave as is)
#       INCLUDES	location for the include files (leave as is)
#       BIN	        location for the executable files (leave as is)
#       RM              command to delete files
#	PRINT		command to send files to printer
#.......................................................
LIBRARY		= $(GRIB_ENV)/lib/libgrib.a
#LIBRARY		= $(GRIB_ENV)/lib/grib.a
INCLUDES 	= $(GRIB_ENV)/include
BIN             = $(GRIB_ENV)/bin
RM 		= /bin/rm -f
PRINT		= /bin/lp

#.......................................................
# The DEBUG flag can be left blank, set to the standard -g, or set to
# -DVERBOSE for lots of additional diagnostic printing from the library.
#.......................................................
DEBUG		=
#DEBUG		= -g -DVERBOSE 

#.......................................................
# Change the two path variables as required to match the vis5d 
# installation on your system.
# Inclusion of VIS5D in gribsimp is controlled at the gribsimp Makefile
# level.  The GRIB library can be compiled independently.
# Leave blank if not available.
#       V5DSRCPATH	location of the Vis5D source files
#       V5DINCPATH	location of the Vis5D include files
#.......................................................
#V5DSRCPATH	= 
#V5DINCPATH	=
V5DSRCPATH      = /usr/local/vis5d-4.3/src
V5DINCPATH      = /usr/local/vis5d-4.3/src

#....................................................
# Prompt the user to make a correct selection
#
default: 
	@echo ""
	@echo "No arguments given at 'make' or 'Install'."
	@echo "Valid Arguments:"
	@echo "    sun_gcc         SunOS (sun4c) version 4.1.3 (v2.6)"
	@echo "    sun_acc         SunOS (sun4c) v4.1.3  (acc v.SC3.0.1)"
	@echo "    sgi_cc          SGI (IP28) IRIX64 v6.5"
	@echo "    v5d_sgi_cc      SGI (IP28) IRIX64 v6.5, with vis5d"

#############################################################################
# PLATFORM & OPERATING SYSTEM SPECIFIC VARIABLES
# Compiler related variables:
#	CC		C-compiler command, used in Makefile
#	CFLAG		options for C-compiler 
#	LDFLAG		options for linker
#	RANLIB		used for building Library on SunOS
#       V5DON           Inclusion of VIS5D in gribsimp is controlled at the 
#			gribsimp Makefile level;  leave as is
#############################################################################
#
# SGI (IP28) IRIX64 version 6.5 , using 'cc' (version 7.2.1)
#
sgi_cc:
	${MAKE} target1  \
	"CC		= cc"  \
	"CFLAG  	= -n32 -xansi $(DEBUG) -I$(INCLUDES)" \
	"LDFLAG 	= -lm" \
	"RANLIB		=" \
	"V5DON  	=" \
	"TGT		= $@"

#
# SGI (IP28) IRIX64 version 6.5 , using 'cc' (version 7.2.1)
# 				  with vis5d extensions included
#
v5d_sgi_cc:
	${MAKE} target2 \
	"CC     	= cc"  \
	"CFLAG  	= -n32 -xansi $(DEBUG) -I$(INCLUDES)" \
	"LDFLAG 	= -lm" \
	"RANLIB		=" \
	"V5DON		= -DV5D_ON -I$(V5DINCPATH)" \
	"TGT		= $@"
#
# -----------------------------------------------------------
# SunOS (sun4c) version 4.1.3  using 'acc' (version SC3.0.1) 
# Must explicitly invoke 'ranlib' under SunOS
#
sun_acc:
	${MAKE} target1  \
	"CC		= acc" \
	"CFLAG		= -Xa $(DEBUG) -I$(INCLUDES)" \
	"LDFLAG		=" \
	"RANLIB		= ranlib $(LIBRARY)" \
	"V5DON  	=" \
	"TGT		= $@"

# -----------------------------------------------------------
# SunOS (sun4c) version 4.1.3   using 'gcc' (version 2.6)
# Must explicitly invoke 'ranlib' under SunOS
#
sun_gcc:
	${MAKE} target1  \
	"CC		= gcc" \
	"CFLAG         	= -ansi $(DEBUG) -I$(INCLUDES)" \
	"LDFLAG		= /usr/lib/libm.a" \
	"RANLIB		= ranlib $(LIBRARY)" \
	"V5DON  	=" \
	"TGT		= $@"

# -----------------------------------------------------------
# Linux version 2.2.14-6.0   using 'gcc' (version egcs-2.91.66)
# 
#
linux_gcc:
	${MAKE} target1  \
	"CC     = gcc"  \
	"CFLAG  = -ansi -Wall -Wmissing-prototypes $(DEBUG) -I$(INCLUDES)" \
	"LDFLAG = -lm" \
	"RANLIB =" \
	"V5DON  =" \
	"TGT	= $@"

# -----------------------------------------------------------
# ***Creating your Own Operating System Configuration***
# Choose the appropriate template depending on whether or not you
# wish to compile with vis5d extensions
# Use Template #1 for the standard compilation for your system
# Template #2 is optional and is used to compile gribsimp with vis5d
#
# Define the macros for the Makefiles:
#  TargetName 
#	is the name you will use to compile the code later
#	ie (Install TargetName or make TargetName)
#
#  ${MAKE} target1  or ${MAKE} target2  
#	Leave this line AS IS
#
#  CC 	is the name of the C-compiler 
#
#  CFLAG is the compiler options;  
#	It is recommended that you use the ANSI-C option on your compiler
#
#  LDFLAG  is the loader options
#
#  RANLIB generates an index to the contents of an archive;  
#	Define it only if your archive/library maintainer (ar) 
#	doesn't call ranlib automatically.  Else leave blank;
#
#  V5DON tells Makefile in src/gribsimp whether to compile with
# 	or without the vis5d extensions;  
#	Does not affect the Library directly;
#	Set to "-DV5D_ON -I$(V5DINCPATH)" to include the Vis5d extensions
#	(V5DINCPATH & V5DSRCPATH must already be defined, and the
#       Vis5d source and include files must exist in these paths);
#       else leave blank;
#
#  TGT is name of the current target.  Leave line AS IS.
# --------------------------------------------------------------
#
#*** TEMPLATE #1:  Compile without 'vis5d' extensions:
#TargetName:
#	${MAKE} target1  \
#	"CC     = cc"  \
#	"CFLAG  = $(DEBUG) -I$(INCLUDES)" \
#	"LDFLAG =" \
#	"RANLIB =" \
#	"V5DON  =" \
#	"TGT	= $@"
 
#*** TEMPLATE #2:  Compile with 'vis5d' extensions:
#Vis5dTargetName:
#	${MAKE} target2 \
#	"CC     = cc"  \
#	"CFLAG  = $(DEBUG) -I$(INCLUDES)" \
#	"LDFLAG =" \
#	"RANLIB =" \
#	"V5DON  = -DV5D_ON -I$(V5DINCPATH)" \
#	"TGT	= $@"


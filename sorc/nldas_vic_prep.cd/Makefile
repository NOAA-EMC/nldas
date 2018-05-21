# Makefile for ncep_forcing_pp

# Greg O'Donnell    May, 2000
# Modified: Justin Sheffield (justin@princeton.edu)
#
# Requires the operating system (OS) as a command line argument
# NOTE: this works for Gnu Make and may not work for other versions
#
# Command line arguments: you can switch certain parts of the code on/off by
# using the LOCATION and FORCING arguments to the make command.
# e.g. make LOCATION=NCEP FORCING=REALTIME
# LOCATION refers to NCEP or elsewhere, such as PU (Princeton University).
# This is required because the pathnames of the forcing files are different at NCEP
# FORCING refers to the forcing dataset (REALTIME or RETRO). This is used
# to switch on/off code that applies to only to a certain dataset (e.g. using
# EDAS solar radiation only for July 1999 to Feb 2000 for the realtime data).
#

CC			= icc
CFLAGS_BASE		= -D${LOCATION} -D${FORCING} -I. -I./Libraries/grib_v3.0.0/include -I/usr/include
CFLAGS_irix  		= $(CFLAGS_BASE) -O3 -n32 -fullwarn
CFLAGS_linux  		= $(CFLAGS_BASE) -O3 -Wall -Wno-unused -Wuninitialized
CFLAGS_linux_ia64 	= $(CFLAGS_BASE) -g -Wall -Wno-unused 
CFLAGS_linux_alpha      = $(CFLAGS_BASE) -O3 -Wall -Wno-unused -mieee
CFLAGS_sunos      	= $(CFLAGS_BASE) -O
# this line may cause problems for old versions of make
CFLAGS 			= $(CFLAGS_$(OS))
#CFLAGS 			= $(CFLAGS_linux)
LIBDIR 			= ./Libraries/grib_v3.0.0/lib
LIBS			= $(LIBDIR)/libgrib.a -lm

HDRS =  preproc.h

OBJS = 	preproc.o	\
	initialize_global.c	\
	get_global_param.c	\
	get_ldas_global_param.c	\
	get_force_type.c	\
	make_dmy_ldas.o		\
	make_fname.o		\
	open_file.o		\
	process_masks.o		\
	write_grib.o		\
	read_grib.o		\
	nrerror.o		\
	qc_data.o		\
	process_precip.o	\
	process_dswrf.o	\
	process_inst2ave.o

SRCS = $(OBJS:%.o=%.c) 

#EXE = VIC-LDAS-preproc-$(FORCING).$(OS)
EXE = nldas_vic_prep

all:
	@if test "$(OS)" = ''; then \
	  make -s os_error; \
	else \
	  make model; \
	fi

default:
	make model

clean::
	/bin/rm -f *.o core $(EXE) log

model: $(OBJS)
	@if test "$(OS)" = ''; then \
	  make -s os_error; \
	else \
	  $(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LIBS); \
	fi

test:
	@if test "$(OS)" = ''; then \
	  make -s os_error; \
        fi
	@echo OS = $(OS)
	@echo CFLAGS = $(CFLAGS)
	@echo LIBS = $(LIBS)	
os_error:
	@echo "Error: you must specify the operating system as an argument"
	@echo "e.g. make OS=irix|linux|linux_alpha|sunos LOCATION=NCEP|PU FORCING=RETRO|REALTIME|NARR_RETRO"


# -------------------------------------------------------------
# depend
# -------------------------------------------------------------
depend: .depend
.depend:	$(SRCS) $(HDRS)
	$(CC) $(CFLAGS) -MM $(SRCS) > $@

clean::
	\rm -f .depend	     

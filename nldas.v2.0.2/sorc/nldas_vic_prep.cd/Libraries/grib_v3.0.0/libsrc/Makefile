###############################################################
# Makefile for building GRIB Decoder & Encoder Library
# Targets (as defined in config.os):
#   sgi_cc	: builds Library using SGI config
#   v5d_sgi_cc  : builds Lib using SGI config w/vis5d extensions
#   sun_acc	: builds Library using SUN config
#   v5d_sun_acc : builds Lib using SUN config w/vis5d extensions
#   sun_gcc	: builds Library using SUN config
#   v5d_sun_gcc : builds Lib using SUN config w/vis5d extensions
#   clean	: cleans objects and library
#   clobber	: cleans corefile, objects and library
#   * DONOT do 'make target1' or 'make target2' 
###############################################################
include $(GRIB_ENV)/config.os
 
# Library Functions
CFILES  =  	FTP_getfile.c \
		apply_bitmap.c \
		display_gribhdr.c \
		gbyte.c \
		grib_dec.c \
		grib_enc.c \
		grib_seek.c \
		gribgetbds.c \
		gribgetbms.c \
		gribgetgds.c \
		gribgetpds.c \
		gribhdr2file.c \
		gribputbds.c \
		gribputgds.c \
		gribputpds.c \
		hdr_print.c \
		init_dec_struct.c \
		init_enc_struct.c \
		init_gribhdr.c \
		init_struct.c \
		ld_dec_lookup.c \
		ld_enc_input.c \
		ld_enc_lookup.c \
		ld_grib_origctrs.c \
		make_default_grbfn.c \
		make_grib_log.c \
		map_lvl.c \
		map_parm.c \
		pack_spatial.c \
		prt_inp_struct.c \
		upd_child_errmsg.c \
	        prt_badmsg.c

OFILES  =  	FTP_getfile.o \
		apply_bitmap.o \
		display_gribhdr.o \
		gbyte.o \
		grib_dec.o \
		grib_enc.o \
		grib_seek.o \
		gribgetbds.o \
		gribgetbms.o \
		gribgetgds.o \
		gribgetpds.o \
		gribhdr2file.o \
		gribputbds.o \
		gribputgds.o \
		gribputpds.o \
		hdr_print.o \
		init_dec_struct.o \
		init_enc_struct.o \
		init_gribhdr.o \
		init_struct.o \
		ld_dec_lookup.o \
		ld_enc_input.o \
		ld_enc_lookup.o \
		ld_grib_origctrs.o \
		make_default_grbfn.o \
		make_grib_log.o \
		map_lvl.o \
		map_parm.o \
		pack_spatial.o \
		prt_inp_struct.o \
		upd_child_errmsg.o \
	        prt_badmsg.o

## --------------------------------------------------------------------------
## Include file dependencies
## --------------------------------------------------------------------------
DPRINTS	= $(OFILES)
GRIB	= $(OFILES)
GRIBFUNCS= $(OFILES)
GRIB_LOOKUP= grib_enc.o ld_dec_lookup.o ld_enc_lookup.o ld_grib_origctrs.o \
	make_default_grbfn.o make_grib_log.o map_lvl.o map_parm.o \
	prt_inp_struct.o
INPUT=  grib_enc.c gribputbds.c gribputgds.c gribputpds.c init_enc_struct.c \
	init_enc_struct.c init_struct.c ld_enc_input.c make_default_grbfn.c \
	map_lvl.c map_parm.c
ISDB= 	gbyte.o pack_spatial.o

$(DPRINTS):	$(INCLUDES)/dprints.h
$(GRIB):	$(INCLUDES)/grib.h
$(GRIBFUNCS):	$(INCLUDES)/gribfuncs.h
$(GRIB_LOOKUP):	$(INCLUDES)/grib_lookup.h
$(INPUT):	$(INCLUDES)/input.h
$(ISDB):	$(INCLUDES)/isdb.h

##################################################################
target1 target2:		$(LIBRARY)
###################################################################
$(LIBRARY):	$(OFILES)
	ar r $@ $(OFILES)
	$(RANLIB)

$(OFILES):	$(GRIB_ENV)/config.os

.c.o:
	$(CC) -c $(CFLAG) $*.c

print_all:
	$(PRINT) Makefile $(LIBINCLUDES) $(CFILES)
clean:
	rm -f *.o
clobber:
	rm -f *.o core $(LIBRARY)

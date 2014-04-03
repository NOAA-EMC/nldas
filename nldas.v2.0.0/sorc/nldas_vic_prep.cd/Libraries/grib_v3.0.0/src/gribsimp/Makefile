###################################################################
#  Makefile for program 'gribsimp'
#  Invoke with either:
#     'make sgi_cc', or 'make v5d_sgi_cc', or 'make sun_acc', or 
#     'make sun_gcc'.  See $GRIB_ENV/config.os for avail choices;
#  NOTE:  'make target1' & 'make target2' will NOT work if since
#     it won't load the macros in config.os under the specific
#     Operating System configuration.  
#
#  If the library has been built w/o vis5d extensions, and you
#  want to include vis5d extensions then all you need to rebuild
#  is gribsimp.  libsrc/*.c do not use V5D so don't need to touch it.
#  Make sure 'config.os' contains all of the V5D environment required, 
#  then from dir $GRIB_ENV/src/gribsimp/, do:  
#	make clean		<< don't forget this step !
#	make v5d_xxx_xx
#
#  Note that the Library will automatically be recompiled if file 
#  'config.os' has been edited since the last Install.
###################################################################

include $(GRIB_ENV)/config.os


PROG=   $(BIN)/gribsimp

SRC= 	gribsimp.c h_e_time.c reg_trfm.c \
	ld_grad_msg.c grad_boundary_box.c make_grad_files.c

V5DSRC= e_h_time.c ld_v5d_msg.c make_v5d_file.c

OBJ= 	gribsimp.o h_e_time.o reg_trfm.o \
	ld_grad_msg.o grad_boundary_box.o make_grad_files.o

#
# SPeCial object files required for including Vis5d extension
#
SPCOBJ= ld_v5d_msg.o make_v5d_file.o e_h_time.o
V5DOBJ= $(V5DSRCPATH)/binio.o $(V5DSRCPATH)/v5d.o


HDRS= 	$(INCLUDES)/isdb.h $(INCLUDES)/grib.h \
 	$(INCLUDES)/grib_lookup.h $(INCLUDES)/dprints.h \
	$(INCLUDES)/gribfuncs.h $(INCLUDES)/grads.h

V5DHDRS=  $(INCLUDES)/gsv5d.h  \
	$(V5DINCPATH)/binio.h $(V5DINCPATH)/v5d.h

target1:
	@echo "**********   Target 1 of gribsimp   **********"
	@make CKLIB
	@make PLAINVERS
	@echo "done"

target2:
	@echo "**********   Target 2 of gribsimp   **********"
	@make CKLIB
	@make V5DVERS
	@echo "done"


CKLIB:
	@echo "*********     Checking Library        *********"
	@(cd $(GRIB_ENV)/libsrc;make $(TGT))


PLAINVERS: 	$(OBJ) $(LIBRARY)
	@echo "*********  Linking gribsimp (no V5d extension)  *********"
	$(CC) $(CFLAG) $(OBJ) $(LIBRARY) -o $(PROG) $(LDFLAG)


V5DVERS:	$(OBJ) $(SPCOBJ) $(LIBRARY)
	@echo "*********  Linking gribsimp WITH V5d extension  *********"
	$(CC) $(CFLAG) $(V5DON) $(OBJ) $(SPCOBJ) $(V5DOBJ) $(LIBRARY) $(LDFLAG) -o $(PROG)

$(OBJ):		$(HDRS) $(LIBRARY)

$(SPCOBJ):	$(HDRS) $(V5DHDRS) $(LIBRARY)

#
# Include the V5DON macro from 'config.os' explicitly 
#
.c.o:		$(HDRS)
	$(CC) -c $(V5DON) $(CFLAG) $*.c
print_all:
	$(PRINT)  Makefile $(HDRS) $(V5DHDRS) $(SRC) $(V5DSRC)
clean:
	$(RM) *.o $(PROG)
clobber:
	$(RM) core *.o $(PROG)

###########################################################
#  Makefile for program 'decoder_ex'
#  call 'make' with 1 of the targets defined in 'config.os'
#  such as sgi_cc, v5d_sgi_cc, sun_acc, sun_gcc ...
###########################################################
include $(GRIB_ENV)/config.os

PROG    =  $(BIN)/decoder_ex
SRC	=  decoder_ex.c
OBJS	=  decoder_ex.o
HDRS	=  $(INCLUDES)/gribfuncs.h

target1 target2:
	@echo "*********     Checking Library        *********"
	(cd $(GRIB_ENV)/libsrc;make $(TGT))
	(cd $(GRIB_ENV)/src/decoder_ex; make $(PROG))

$(PROG): 	$(OBJS)
	@echo ""
	@echo "*********     Linking Program         *********"
	cd $(GRIB_ENV)/src/decoder_ex
	$(CC) $(CFLAG) $(OBJS) $(LIBRARY) -o $(PROG) $(LDFLAG)
	@echo "done"


$(OBJS):  $(HDRS) $(LIBRARY)

.c.o:
	$(CC) -c $(CFLAG) $*.c

print_all:
	$(PRINT) Makefile $(HDRS) $(SRC) 

clean:
	$(RM) *.o $(PROGS)

clobber:
	$(RM) core *.o $(PROGS)

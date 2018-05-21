###########################################################
#  Makefile for program 'getgribieee'
#  call 'make' with 1 of the targets defined in 'config.os'
#  such as sgi_cc, v5d_sgi_cc, sun_acc, sun_gcc ...
###########################################################
include $(GRIB_ENV)/config.os

PROG    =  $(BIN)/getgribieee
SRC	=  getgribieee.c
OBJS	=  getgribieee.o
HDRS	=  $(INCLUDES)/gribfuncs.h

target1 target2:
	@echo "*********     Checking Library        *********"
	(cd $(GRIB_ENV)/libsrc;make $(TGT))
	(cd $(GRIB_ENV)/src/getgribieee; make $(PROG))

$(PROG):    $(OBJS)
	@echo "*********     Linking Program         *********"
	$(CC) $(CFLAG) $(OBJS) $(LIBRARY) -o $(PROG) $(LDFLAG)
	@echo "done"


$(OBJS):    $(HDRS) $(LIBRARY)

.c.o:
	@echo "Compiling $@"
	$(CC) -c $(CFLAG) $*.c

print_all:
	$(PRINT) Makefile $(HDRS) $(SRC) 

clean:
	@$(RM) *.o $(PROGS)

clobber:
	@$(RM) core *.o $(PROGS)

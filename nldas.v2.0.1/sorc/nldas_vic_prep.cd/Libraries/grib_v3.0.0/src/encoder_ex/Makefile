###################################################################
#  Makefile for program 'encoder_ex'
#  Invoke with either:
#     'make sgi_cc', or 'make v5d_sgi_cc', or 'make sun_acc', or 
#     'make sun_gcc'.  See $GRIB_ENV/config.os for avail choices;
#  NOTE:  'make target1' & 'make target2' will NOT work if since
#     it won't load the macros in config.os under the specific
#     Operating System configuration.  
###################################################################

include $(GRIB_ENV)/config.os

PROGS= 	$(BIN)/encoder_ex1 $(BIN)/encoder_ex2 $(BIN)/encoder_ex3
OBJS= 	$(BIN)/encoder_ex1.o $(BIN)/encoder_ex2.o $(BIN)/encoder_ex3.o
SRC=   	encoder_ex1.c encoder_ex2.c encoder_ex3.c
HDRS= 	$(INCLUDES)/grib_lookup.h $(INCLUDES)/dprints.h $(INCLUDES)/gribfuncs.h


target1 target2:  
	(make CKLIB)
	(make $(BIN)/encoder_ex1)
	(make $(BIN)/encoder_ex2)
	(make $(BIN)/encoder_ex3)

CKLIB:
	@echo "*********     Checking Library        *********"
	(cd $(GRIB_ENV)/libsrc;make $(TGT))

$(BIN)/encoder_ex1:	encoder_ex1.o $(LIBRARY)
	@echo "*********	Linking encoder_1	*********"
	$(CC) $(CFLAG) encoder_ex1.o $(LIBRARY) -o $(BIN)/encoder_ex1 $(LDFLAG)
	@echo "done"

$(BIN)/encoder_ex2:	encoder_ex2.o $(LIBRARY)
	@echo "*********	Linking encoder_2	*********"
	$(CC) $(CFLAG) encoder_ex2.o $(LIBRARY) -o $(BIN)/encoder_ex2 $(LDFLAG)
	@echo "done"

$(BIN)/encoder_ex3:	encoder_ex3.o $(LIBRARY)
	@echo "*********	Linking encoder_ex3	*********"
	$(CC) $(CFLAG) encoder_ex3.o $(LIBRARY) -o $(BIN)/encoder_ex3 $(LDFLAG)
	@echo "done"

$(OBJS):	$(HDRS) $(LIBRARY)

.c.o:
	$(CC) -c $(CFLAG) $*.c

print_all:
	$(PRINT)  Makefile $(HDRS) $(SRC)
clean:
	$(RM) *.o $(PROGS)
clobber:
	$(RM) core *.o $(PROGS)

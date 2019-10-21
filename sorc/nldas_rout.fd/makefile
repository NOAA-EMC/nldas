OBJS = 	rout.o  \
	rout_it.o \
	readcntl.o \
	READ_REAL_SPATIAL.o \
	READ_INT_SPATIAL.o \
	CALC_LENGTH.o\
	READ_RUNOFF.o\
	READ_GRIB.o\
        grib2_wrt_g2func.o\
	calc_area.o\
        MAKEPDS.o\
        GRIBOUT.o
# Tunable parameters 
#
# COMPILER      Compiler name
# FFLAGS	Flags to the fortran compiler
# LIBS		A list of libraries to use 
# PROGRAM	Name of the executable

COMPILER=ifort
INC = $(G2_INC4)
LIBS=   ${G2_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4} ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}

PROGRAM=nldas_rout
FFLAGS= -O3 -extend-source -convert big_endian -I$(INC) -c

$(PROGRAM): $(OBJS)
	    $(COMPILER) $(OBJS) $(LIBS) -o $(PROGRAM)

clean:
	/bin/rm -rf $(OBJS) $(EXE)

rout.o: rout.f
	$(COMPILER) $(FFLAGS) rout.f -o rout.o

rout_it.o: rout_it.f
	$(COMPILER) $(FFLAGS) rout_it.f -o rout_it.o

readcntl.o: readcntl.f
	$(COMPILER) $(FFLAGS) readcntl.f -o readcntl.o

READ_REAL_SPATIAL.o: READ_REAL_SPATIAL.f
	$(COMPILER) $(FFLAGS) READ_REAL_SPATIAL.f -o READ_REAL_SPATIAL.o

READ_INT_SPATIAL.o: READ_INT_SPATIAL.f
	$(COMPILER) $(FFLAGS) READ_INT_SPATIAL.f -o READ_INT_SPATIAL.o

CALC_LENGTH.o: CALC_LENGTH.f
	$(COMPILER) $(FFLAGS) CALC_LENGTH.f -o CALC_LENGTH.o

READ_RUNOFF.o: READ_RUNOFF.f
	$(COMPILER) $(FFLAGS) READ_RUNOFF.f -o READ_RUNOFF.o

READ_GRIB.o: READ_GRIB.f
	$(COMPILER) $(FFLAGS) READ_GRIB.f -o READ_GRIB.o

calc_area.o: calc_area.f
	$(COMPILER) $(FFLAGS) calc_area.f -o calc_area.o

MAKEPDS.o: MAKEPDS.f
	$(COMPILER) $(FFLAGS) MAKEPDS.f -o MAKEPDS.o

GRIBOUT.o: GRIBOUT.f
	$(COMPILER) $(FFLAGS) GRIBOUT.f -o GRIBOUT.o

grib2_wrt_g2func.o: grib2_wrt_g2func.f
	$(COMPILER) $(FFLAGS) grib2_wrt_g2func.f -o grib2_wrt_g2func.o




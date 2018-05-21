OBJS = 	precforce.o  \
	update.o \
        gridcpc.o 
# Tunable parameters 
#
# COMPILER      Compiler name
# FFLAGS	Flags to the fortran compiler
# LIBS		A list of libraries to use 
# PROGRAM	Name of the executable

COMPILER=ifort

INC =$(G2_INC4)

LIBS=${G2_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4} ${IP_LIB4} ${SP_LIB4} -openmp ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}

PROGRAM=precforce
FFLAGS= -O3 -extend-source -convert big_endian -I$(INC) -c

$(PROGRAM): $(OBJS)
	    $(COMPILER) $(OBJS) $(LIBS) -o $(PROGRAM)

clean:
	/bin/rm -rf $(OBJS) $(EXE)

precforce.o: precforce.f
	$(COMPILER) $(FFLAGS) precforce.f -o precforce.o

update.o: update.f
	$(COMPILER) $(FFLAGS) update.f -o update.o

gridcpc.o: gridcpc.f
	$(COMPILER) $(FFLAGS) gridcpc.f -o gridcpc.o



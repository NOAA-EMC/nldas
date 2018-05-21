OBJ =   obj/main_sac.o  	\
	obj/READCNTL.o		\
	obj/READ_PAR.o		\
	obj/READ_INITIAL.o	\
	obj/CHECK_INITIAL.o	\
	obj/READ_FORCING.o	\
	obj/READ_NOAH_GRIB.o	\
        obj/ZENANGL.o           \
	obj/CHECK_FORCING.o	\
	obj/QDATAP.o		\
	obj/GRIBOUT.o		\
        obj/grib2_wrt_g2func.o	\
	obj/CALC_WEIGTHS.o	\
        obj/SKIP_COMMENTS.o     \
        obj/E.o 	    	\
        obj/DQSDT.o     	\
        obj/DQS.o	     	\
	obj/ex_sac1.o		\
	obj/exsnow19.o		\
	obj/sac1.o		\
	obj/pack19.o		\
	obj/updt19.o		\
	obj/adjc19.o		\
	obj/aeco19.o		\
	obj/aesc19.o		\
	obj/melt19.o		\
	obj/rout19.o		\
	obj/snowpack.o		\
	obj/sndepth.o		\
	obj/snew.o		\
	obj/snowt.o		\
	obj/penman1.o		\
	obj/zero19.o

# Tunable parameters 
#
# COMPILER      Compiler name
# FFLAGS	Flags to the fortran compiler
# LIBS		A list of libraries to use 
# PROGRAM	Name of the executable

COMPILER = ifort -convert big_endian
INC =           /nwprod/lib/incmod/g2_4
#LIBS=   /nwprod/lib/libip_d.a \
#        /nwprod/lib/libbacio_4.a \
#        /nwprod/lib/libw3nco_4.a
LIBS=   -L /nwprod/lib -lg2_4 -lw3nco_4 -lw3emc_4 -lbacio_4 -ljasper -lpng -lz
                            
PROGRAM= nldas_sac_ldas
FFLAGS= -O3 -extend-source -I$(INC) -c 

# These are the
# actual rules which make uses to build corr.
# 
# The first line says that the program depends on the object
# files.  If the program is out of date or does not
# exist then the objects have to be linked together using
# the compiler (linker)
#
# The subsequent lines describe the dependencies for each
# object file.  Each object file depends on its source
# file and 0 or more include files.  If any include file
# or the source file have been modified, the source
# file is recompiled 

$(PROGRAM): $(OBJ)
	    $(COMPILER) $(OBJ) $(LIBS) -o $(PROGRAM)

clean:
	/bin/rm -f obj/*.o $(PROGRAM)
obj/GRIBOUT.o: DRIVER/GRIBOUT.f
	$(COMPILER) $(FFLAGS) DRIVER/GRIBOUT.f -o obj/GRIBOUT.o

obj/grib2_wrt_g2func.o:  DRIVER/grib2_wrt_g2func.f
	$(COMPILER) $(FFLAGS) DRIVER/grib2_wrt_g2func.f -o obj/grib2_wrt_g2func.o

obj/E.o: DRIVER/E.f
	$(COMPILER) $(FFLAGS) DRIVER/E.f -o obj/E.o

obj/DQSDT.o: DRIVER/DQSDT.f
	$(COMPILER) $(FFLAGS) DRIVER/DQSDT.f -o obj/DQSDT.o

obj/DQS.o: DRIVER/DQS.f
	$(COMPILER) $(FFLAGS) DRIVER/DQS.f -o obj/DQS.o

obj/CALC_WEIGTHS.o: DRIVER/CALC_WEIGTHS.f
	$(COMPILER) $(FFLAGS) DRIVER/CALC_WEIGTHS.f -o obj/CALC_WEIGTHS.o

obj/SKIP_COMMENTS.o: DRIVER/SKIP_COMMENTS.f
	$(COMPILER) $(FFLAGS) DRIVER/SKIP_COMMENTS.f -o obj/SKIP_COMMENTS.o

obj/main_sac.o: DRIVER/main_sac.f DRIVER/main_sac.h
	$(COMPILER) $(FFLAGS) DRIVER/main_sac.f -o obj/main_sac.o

obj/QDATAP.o: DRIVER/QDATAP.f
	$(COMPILER) $(FFLAGS) DRIVER/QDATAP.f -o obj/QDATAP.o

obj/READCNTL.o: DRIVER/READCNTL.f
	$(COMPILER) $(FFLAGS) DRIVER/READCNTL.f -o obj/READCNTL.o

obj/READ_PAR.o: DRIVER/READ_PAR.f
	$(COMPILER) $(FFLAGS) DRIVER/READ_PAR.f -o obj/READ_PAR.o

obj/READ_FORCING.o: DRIVER/READ_FORCING.f
	$(COMPILER) $(FFLAGS) DRIVER/READ_FORCING.f -o obj/READ_FORCING.o

obj/READ_NOAH_GRIB.o: DRIVER/READ_NOAH_GRIB.f
	$(COMPILER) $(FFLAGS) DRIVER/READ_NOAH_GRIB.f -o obj/READ_NOAH_GRIB.o

obj/ZENANGL.o: DRIVER/ZENANGL.f
	$(COMPILER) $(FFLAGS) DRIVER/ZENANGL.f -o obj/ZENANGL.o

obj/READ_INITIAL.o: DRIVER/READ_INITIAL.f
	$(COMPILER) $(FFLAGS) DRIVER/READ_INITIAL.f -o obj/READ_INITIAL.o

obj/CHECK_INITIAL.o: DRIVER/CHECK_INITIAL.f
	$(COMPILER) $(FFLAGS) DRIVER/CHECK_INITIAL.f -o obj/CHECK_INITIAL.o

obj/CHECK_FORCING.o: DRIVER/CHECK_FORCING.f
	$(COMPILER) $(FFLAGS) DRIVER/CHECK_FORCING.f -o obj/CHECK_FORCING.o

obj/zero19.o: SOURCE/zero19.f
	$(COMPILER) $(FFLAGS) SOURCE/zero19.f -o obj/zero19.o

obj/snowpack.o: SOURCE/SNOWPACK.f
	$(COMPILER) $(FFLAGS) SOURCE/SNOWPACK.f -o obj/snowpack.o

obj/sndepth.o: SOURCE/SNDEPTH.f
	$(COMPILER) $(FFLAGS) SOURCE/SNDEPTH.f -o obj/sndepth.o

obj/snew.o: SOURCE/SNEW.f
	$(COMPILER) $(FFLAGS) SOURCE/SNEW.f -o obj/snew.o

obj/snowt.o: SOURCE/SNOWT.f
	$(COMPILER) $(FFLAGS) SOURCE/SNOWT.f -o obj/snowt.o

obj/rout19.o: SOURCE/rout19.f
	$(COMPILER) $(FFLAGS) SOURCE/rout19.f -o obj/rout19.o

obj/melt19.o: SOURCE/melt19.f
	$(COMPILER) $(FFLAGS) SOURCE/melt19.f -o obj/melt19.o

obj/adjc19.o: SOURCE/adjc19.f
	$(COMPILER) $(FFLAGS) SOURCE/adjc19.f -o obj/adjc19.o

obj/aeco19.o: SOURCE/aeco19.f
	$(COMPILER) $(FFLAGS) SOURCE/aeco19.f -o obj/aeco19.o

obj/aesc19.o: SOURCE/aesc19.f
	$(COMPILER) $(FFLAGS) SOURCE/aesc19.f -o obj/aesc19.o

obj/updt19.o: SOURCE/updt19.f
	$(COMPILER) $(FFLAGS) SOURCE/updt19.f -o obj/updt19.o

obj/pack19.o: SOURCE/PACK19.f
	$(COMPILER) $(FFLAGS) SOURCE/PACK19.f -o obj/pack19.o

obj/sac1.o: SOURCE/sac1.f
	$(COMPILER) $(FFLAGS) SOURCE/sac1.f -o obj/sac1.o

obj/exsnow19.o: SOURCE/exsnow19.f
	$(COMPILER) $(FFLAGS) SOURCE/exsnow19.f -o obj/exsnow19.o

obj/penman1.o: SOURCE/penman1.f
	$(COMPILER) $(FFLAGS) SOURCE/penman1.f -o obj/penman1.o

obj/ex_sac1.o: SOURCE/ex_sac1.f
	$(COMPILER) $(FFLAGS) SOURCE/ex_sac1.f -o obj/ex_sac1.o

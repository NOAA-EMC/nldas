# VIC postprocessor Makefile
# Last Changed: 13 March 2013 by Youlong Xia <Youlong.Xia@noaa.gov>
#   make OS=linux_ia64 to compile post-process code
FC = ifort -132
FFLAGS_BASE		= -I. 
## for linux x86_64 machine
FFLAGS_linux_ia64 	= $(FFLAGS_BASE) -O3
LIBRARY = ${G2_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4} ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}

HDRS = 	postproc.h           

OBJS =  main.o		\
        usage.o       	\
        get_options.o 	\
        get_mask.o    	\
        get_vars.o    	\
        read_bin.o	\
        make_dates.o  	\
	gribout.o	\
	grib2_wrt_g2func_vic.o	\
	daybefore.o	

SRCS = $(OBJS:%.o=%.f) 

EXE = nldas_vic_post

all:
	@if test "$(OS)" = ''; then \
	  make os_error; \
	else \
	  make model; \
	fi
default:
	@make model

clean::
	@/bin/rm -f *.o core log *~	
model: $(OBJS)
	@if test "$(OS)" = ''; then \
	  make os_error; \
	else \
	  $(FC) -o $(EXE) $(OBJS) $(FFLAGS) $(LIBRARY); \
	fi	

# -------------------------------------------------------------
# tags
# so we can find our way around
# -------------------------------------------------------------
tags:	TAGS
TAGS:	$(SRCS) $(HDRS)
	etags $(SRCS) $(HDRS)
clean::
	\rm -f TAGS	       


# -------------------------------------------------------------
# depend
# -------------------------------------------------------------
depend: .depend
.depend:	$(SRCS) $(HDRS)
	$(FC) $(FFLAGS) -MM $(SRCS) > $@

clean::
	\rm -f .depend	     

$(OBJS): $(HDRS)

LIB	= libw3.a

# Make sure one of the following options appears in your CFLAGS 
# variable to indicate which system you are on.
# -DLINUX, -DSGI, -DHP, -DCRAY90, -DIBM4, -DIBM8, -DLINUXF90 

# OPTIONS FOR IBM-SP 
#F77	= xlf
#FFLAGS	= -O3 -q64 -qnosave -qarch=auto
#CFLAGS	= -DIBM4 -q64
#ARFLAGS	= -X64 -r
#CC	= cc

# OPTIONS FOR LINUX LAHEY
#F77	= lf95
#FFLAGS	= 
#CFLAGS	= -O -DLINUX
#CC	= cc

# OPTIONS FOR LINUX PGI
#F77	= pgf90
#FFLAGS	= -Mbyteswapio -r4 -i4 -Mpreprocess 
#CFLAGS	= -O -DLINUX
#CC	= cc

# OPTIONS FOR LINUX ABSOFT
#F77	= f90
#FFLAGS	= -YEXT_NAMES=LCS -s -B108 -YCFRL=1
#CFLAGS	= -O -DLINUX
#CC	= gcc
#ARFLAGS	= r

# OPTIONS FOR SGI
#F77     = f90 -64
#CFLAGS  = -O -64 -DSGI
#CC      = cc
#ARFLAGS = r

# OPTIONS FOR INTEL FORTRAN COMPILER IA64
#F77     = efc
#CFLAGS      := -DLINUX
#FFLAGS      := -c -arch pn4 -tune pn4 -nomixed_str_len_arg -names lowercase -convert big_endian -assume byterecl -DHIDE_SHR_MSG -DNO_SHR_VMATH -DIFC
#CC      = gcc
#ARFLAGS = r


# OPTIONS FOR INTEL FORTRAN COMPILER
#F77     = /opt/intel_fc_80/bin/ifort
F77     = ifort
FFLAGS  = -nomixed_str_len_arg -names lowercase
CFLAGS  = -O -DLINUX
CC      = gcc
ARFLAGS = r

# OPTIONS FOR INTEL FORTRAN COMPILER IA64
#F77     = efc
#FFLAGS  = 
#CFLAGS  = -O -DLINUX
#CC      = gcc
#ARFLAGS = r

# OPTIONS FOR VPP5000
#F77     = frt
#FFLAGS  = -X9 -Ad -Sw -Am -Of -Wv,-Of
#CFLAGS  = -O -K4 -DVPP5000
#CC      = vcc

.SUFFIXES: .o .f .F .c
#
# *****************************************************************
#
OBJ_MOD =  bacio_module.o

OBJS	= \
getgb.o     getgbmh.o   putgbex.o   w3fi73.o     baciof.o \
getgb1r.o   getgbmp.o   putgbn.o    w3fi74.o \
getgb1re.o  getgbp.o    r63w72.o    w3fi75.o \
getgb1s.o   getgi.o     sbyte.o     w3fi76.o \
getgbe.o    getgir.o    sbytes.o    w3fi82.o \
getgbeh.o   idsdef.o    skgb.o      w3fi83.o \
getgbem.o   iw3jdn.o    w3difdat.o  w3fs26.o \
getgbemh.o  ixgb.o      w3doxdat.o  w3log.o \
getgbemp.o  lengds.o    w3fi01.o    w3movdat.o \
getgbens.o  pdsens.o    w3fi58.o    w3reddat.o \
getgbep.o   pdseup.o    w3fi59.o    w3tagb.o \
errmsg.o    getgbex.o   putgb.o     w3fi63.o    \
gbytes.o                getgbexm.o  putgbe.o    w3fi68.o    \
gbytes_char.o           getgbh.o    putgben.o   w3fi71.o \
getbit.o                getgbm.o    putgbens.o  w3fi72.o \
errexit.o fparsei.o fparser.o instrument.o start.o summary.o w3utcdat.o \
w3fs21.o w3locdat.o

OBJS_CC= bacio.v1.3.o mova2i.o

SRC	= $(OBJ_MOD:.o=.f) $(OBJS:.o=.f) $(OBJS_CC:.o=.c)
#
# *****************************************************************
#

#$(LIB): $(OBJ_MOD) $(OBJS) $(OBJS_CC) 
#	ar $(ARFLAGS) -ruv $(LIB) $(OBJ_MOD) $(OBJS) $(OBJS_CC)

$(LIB): $(OBJ_MOD) $(OBJS) $(OBJS_CC) 
	ar $(ARFLAGS) $(LIB) $(OBJ_MOD) $(OBJS) $(OBJS_CC)

clean:
	rm -f $(OBJ_MOD) $(OBJS) $(OBJS_CC) *.mod *.a

.F.o:
	$(F77) $(FFLAGS) -c -d $< 
.f.o:
	$(F77) $(FFLAGS) -c $< 
.c.o:
	$(CC) $(CFLAGS) $(DEFS) -c $< 

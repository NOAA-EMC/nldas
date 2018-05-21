sorc=convert_cpcops_cpcus_cpcglobe.f
pgrm=cpc_precip_convert

lib=${IP_LIB4} ${SP_LIB4} ${BACIO_LIB4} ${W3NCO_LIB4} 
 
FC=ifort
FFLAGS=-openmp -extend_source -convert big_endian

${pgrm}:	${sorc}
	${FC} ${FFLAGS} -o ${pgrm} ${sorc} ${lib}

clean:
	rm ${pgrm}

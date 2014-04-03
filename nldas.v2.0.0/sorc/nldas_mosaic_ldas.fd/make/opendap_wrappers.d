opendap_wrappers.o opendap_wrappers.d : opendap_wrappers.F90
opendap_wrappers.o : spmdMod.o
opendap_wrappers.o : opendap_module.o

wrap_mpi.o wrap_mpi.d : wrap_mpi.F90
wrap_mpi.o : misc.h
wrap_mpi.o : precision.o
wrap_mpi.o : mpishorthand.o

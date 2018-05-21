      
      logical VERBOSE, DEBUG
      parameter (VERBOSE=.TRUE., DEBUG=.FALSE.)

      integer EXIT_FAILURE
      integer MISSING
      integer MAXRECS, MAXVARS
      real TOL, UNDEF
      parameter (EXIT_FAILURE=1)
      parameter (MISSING=-1)
      parameter (MAXRECS=26292, MAXVARS=48)
      parameter (TOL=1E-10, UNDEF=1e+20)

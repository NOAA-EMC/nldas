C-----------------------------------------------------------------------
C get_vars: 	reads in the list of variables from the variables 
C		database file
C-----------------------------------------------------------------------

      subroutine get_vars(vars_file, var_name, var_mult,
     &                    var_offset, var_rate, nvars)

      implicit none

      include "postproc.h"

C inputs variables
      character*256 vars_file

C input/output variables
      integer nvars
      character*40 var_name(MAXVARS)
      real var_mult(MAXVARS), var_offset(MAXVARS)
      logical*1 var_rate(MAXVARS)

C local variables
      integer ivar, status
      integer idum1, idum2
      character*40 str, cdum1, cdum2, cdum3, cdum4, 
     &                  cdum5, cdum6, cdum7

C off we go
      if(DEBUG) write(*,'(a)') 'get_vars'

C open the variables database file
      open(unit=10,file=vars_file, status='old', iostat=status)
      if( status .gt. 0 ) then
        write(*,*) 'ERROR: Can''t open vars file: ',trim(vars_file)
        stop
      endif
      if(VERBOSE) then
        write(*,'(a,a)') 'Opened vars file: ',trim(vars_file)
      endif     

C skip header
      read(10,*)
      read(10,*)
     
C read in the variable attributes and calculate the total number of vars
      ivar=1
   10 read(10,*,end=999) var_name(ivar), 
     & cdum1, idum1, idum2, cdum2, cdum3, cdum4,
     & var_mult(ivar), var_offset(ivar), cdum7

C translate the rate into T or F
        if(cdum7 .eq. 'NONE' .or. cdum7 .eq. 'F' 
     &                       .or. cdum7 .eq. 'FALSE') then
          var_rate(ivar) = .FALSE.
        else if(cdum7 .eq. 'T' 
     &                       .or. cdum7 .eq. 'TRUE') then
          var_rate(ivar) = .TRUE.
        else
          write(*,*) 'ERROR: Unknown rate type ',cdum7,
     &               ' in vars file ',vars_file
        endif

C write the variable attributes
        if(DEBUG) then
          write(*,*) 'var # ',ivar
          write(*,*) 'name = ',var_name(ivar)
          write(*,*) 'grib name = ',cdum1
          write(*,*) 'grib id = ',idum1
          write(*,*) 'grib level = ',idum2
          write(*,*) 'desc = ',cdum2
          write(*,*) 'grib units = ',cdum3
          write(*,*) 'units = ',cdum4
          write(*,*) 'mult = ',var_mult(ivar)
          write(*,*) 'offset = ',var_offset(ivar) 
          write(*,*) 'rate = ',var_rate(ivar)
        endif

C       increase the total number of vars
        ivar=ivar+1

      goto 10

  999 continue
      close (10)
      nvars = ivar - 1

C check that our array sizes are ok
      if(nvars .gt. MAXVARS) then
        write(*,*) 'ERROR: nvars (',nvars,') > MAXVARS (',MAXVARS,')'
        stop
      endif

C write out the variable names
      if(VERBOSE) write(*,'(a,i)') 'nvars = ',nvars
      if(VERBOSE) then
        do ivar=1,nvars
          write(*,'(i2,x,a20,f12.3,f12.3,l4)') ivar, var_name(ivar), 
     &           var_mult(ivar), var_offset(ivar), var_rate(ivar)
        enddo
      endif

      end

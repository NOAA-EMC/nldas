C-----------------------------------------------------------------------
C get_options: 	reads in the options from the global control file 
C-----------------------------------------------------------------------

      subroutine get_options(ctrl_file, nrecs, dt, vic_output,
     &           startyear, startmonth, startday, starthour,
     &           endyear, endmonth, endday, endhour, 
     &           statehour, lats_table, vars_table,
     &           mask_file, pmask_file, pmask_value, pp_output)
      
      implicit none

      include 'postproc.h'

c input vars
      character*256 ctrl_file
      
c input/output vars
      character*256 vic_output, mask_file, pmask_file, 
     &              pp_output, lats_table, vars_table
      integer pmask_value, statehour
      integer startyear, startmonth, startday, starthour
      integer endyear, endmonth, endday, endhour
      integer nrecs, dt

C local vars
      character*256 cmdstr,optstr,optvalue,dummy,dummy2,getstr
      character*1 char
      integer status, strlen

C off we go
      if(DEBUG) write(*,'(a)') 'get_options'
      
C initialise some variables
      nrecs=MISSING
      endyear=MISSING
      endmonth=MISSING
      endday=MISSING
      mask_file = ''
      pmask_file = ''
      lats_table = ''
      pp_output = ''
      vic_output = ''
      vars_table = ''

C open the model control file
      open(10, file=ctrl_file, status='old', iostat=status)
      if( status .gt. 0 ) then
        write(*,*) 'ERROR: Can''t open control file: ', trim(ctrl_file)
        stop
      endif
      if(VERBOSE) then
        write(*,'(a,a)') 'Opened control file: ',trim(ctrl_file)
      endif
      
C read a line
   10 read(10, '(a)', end=999) cmdstr     
   
C check for comment lines
      if(cmdstr(1:1).ne.'#' .and. cmdstr(1:1).ne.'') then

C extract the option name from the line
        optstr = getstr(cmdstr, 1)
        optvalue = getstr(cmdstr, 2)

        if(DEBUG) then
          write(*,*) 'cmdstr: ',trim(cmdstr)
          write(*,*) 'optstr: ',trim(optstr)
          write(*,*) 'optvalue: ', trim(optvalue)
        endif
        
C extract the option value from the line
        if(optstr .eq. 'TIME_STEP') then
           read(optvalue, *) dt
        else if(optstr .eq. 'STARTYEAR') then
           read(optvalue, *) startyear       
        else if(optstr .eq. 'STARTMONTH') then
           read(optvalue, *) startmonth       
        else if(optstr .eq. 'STARTDAY') then
           read(optvalue, *) startday       
        else if(optstr .eq. 'STARTHOUR') then
           read(optvalue, *) starthour       
        else if(optstr .eq. 'ENDYEAR') then
           read(optvalue, *) endyear       
        else if(optstr .eq. 'ENDMONTH') then
           read(optvalue, *) endmonth       
        else if(optstr .eq. 'ENDDAY') then
           read(optvalue, *) endday       
        else if(optstr .eq. 'ENDHOUR') then
           read(optvalue, *) endhour       
        else if(optstr .eq. 'RESULT_DIR') then
           vic_output = optvalue      
        else if(optstr .eq. 'NRECS') then
           read(optvalue, *) nrecs       
        else if(optstr .eq. 'MASK') then
           mask_file = optvalue
        else if(optstr .eq. 'PROCESSOR_MASK') then
           pmask_file = optvalue   
        else if(optstr .eq. 'PROCESSOR') then
           read(optvalue, *) pmask_value        
        else if(optstr .eq. 'PP_OUTPUT') then
           pp_output = optvalue     
        else if(optstr .eq. 'LATS_TABLE') then
           lats_table  = optvalue
        else if(optstr .eq. 'VARS_TABLE') then
           vars_table = optvalue        
        else if(optstr .eq. 'ENDHOUR') then
           read(optvalue, *) endhour        
        else if(optstr .eq. 'STATEHOUR') then
           read(optvalue, *) statehour       
        else if(optstr .eq. 'TIME_STEP') then
           read(optvalue, *) dt       
        else
          if(DEBUG) write(*,'(a,a)') 'Unknown option: ', trim(optstr)
        endif    
        
      endif
        
      goto 10  
      
  999 continue
      close (10)

C check all the options are present
      if(pmask_value .eq. 0) then
         write(*,*) 'ERROR: Must specify PROCESSOR in global file'
         stop 1
      else if(mask_file .eq. '') then
         write(*,*) 'ERROR: Must specify MASK in global file'
         stop 1
      else if(pmask_file .eq. '') then
         write(*,*) 'ERROR: Must specify PROCESSOR_MASK in global file'
         stop 1
      else if(lats_table .eq. '') then
         write(*,*) 'ERROR: Must specify LATS_TABLE in global file'
         stop 1
      else if(vars_table .eq. '') then
         write(*,*) 'ERROR: Must specify VARS_TABLE in global file'
         stop 1
      else if(pp_output .eq. '') then
         write(*,*) 'ERROR: Must specify PP_OUTPUT in global file'
         stop 1
      else if(vic_output .eq. '') then
         write(*,*) 'ERROR: Must specify RESULT_DIR in global file'
         stop 1
      else if(endhour .lt. 0) then
         write(*,*) 'ERROR: Must specify ENDHOUR in global file'
         stop 1
      endif

      end

C-----------------------------------------------------------------------
C getstr: 	gets the next string from within a larger string 
C		NOTE: doesn't recognise tabs so replace tabs with
C		spaces using sed
C-----------------------------------------------------------------------
      function getstr(str, which)

      implicit none

      character*256 str, getstr
      integer which, i, j, word

C check the input args
      if(which .le. 0) then
        write(*,*) 'Error: String number in function getstr ',
     &             'must be > 0'
        stop
      endif

C initialise
      i = 1
      j = 1

C loop until we get to our word
      do word=1,which

C strip leading blanks
        do
          if(str(i:i) .ne. ' ') exit
          i = i + 1
        enddo

C get non-blanks
        do
          if(str(i:i) .eq. ' ') exit

C if it's our word then make a copy
          if(word .eq. which) then
            getstr(j:j) = str(i:i)
            j = j + 1
          endif

          i = i + 1
        enddo

      enddo

C fill rest of string with blanks
      do i=j,256
        getstr(i:i) = ' '
      enddo

      end

C-----------------------------------------------------------------------
C strlen: 	calculates the length of a string 
C-----------------------------------------------------------------------
      function strlen(str)

      implicit none

      character*256 str
      integer       strlen

      strlen = 1
 100  continue
      if (str(strlen:strlen) .ne. ' ') then
         strlen = strlen + 1
      else
         goto 200
      end if
      goto 100
 200  continue

      strlen = strlen - 1

      end

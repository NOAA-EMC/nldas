C-----------------------------------------------------------------------
C make_dates: 	generates the date structure (y,m,d,h) for each timestep
C 		in the simulation and stores them in arrays
C-----------------------------------------------------------------------

      subroutine make_dates(nrecs, 
     &   global_startyear, global_startmonth, global_startday,
     &   global_starthour,
     &   global_endyear, global_endmonth, global_endday,
     &   global_endhour, dt,
     &   dmy_hour, dmy_day, dmy_month, dmy_year, dmy_day_in_year)

      implicit none
      
      include 'postproc.h'

C input variables
      integer nrecs
      integer global_startyear, global_startmonth
      integer global_startday, global_starthour
      integer global_endyear, global_endmonth
      integer global_endday, global_endhour
      integer dt

C input/output variables
      integer dmy_hour(MAXRECS), dmy_day(MAXRECS), dmy_month(MAXRECS)
      integer dmy_year(MAXRECS), dmy_day_in_year(MAXRECS)

C local variables
      integer startyear, startmonth, startday, starthour
      integer endyear, endmonth, endday, endhour
 
      integer hour, day, month, year
      integer tmphour, tmpday, tmpmonth, tmpyear, tmpjday
      logical leapyr, DONE
      integer ii, offset, jday

      integer days(12)
      data days /31,28,31,30,31,30,31,31,30,31,30,31/

C let's go
      if(DEBUG) write(*,'(a)') 'make_dates'

C initialise
      hour = global_starthour
      year = global_startyear
      day = global_startday
      month = global_startmonth

C check if user defined end date or number of records
      if(nrecs .lt. 0) then

        if(global_endyear.lt.0 .or. global_endmonth.lt.0 .or.
     &     global_endday.lt.0 .or. global_endhour.lt.0) then
          write(*,*) 'The model file MUST define EITHER the number'
          write(*,*) 'of records to simulate (NRECS), or the year '
          write(*,*) '(ENDYEAR), month (ENDMONTH), day (ENDDAY) and '
          write(*,*) 'hour (ENDHOUR) of the last full simulation day'
          stop
        endif

        endday = global_endday
        endmonth = global_endmonth
        endyear = global_endyear
        endhour = global_endhour

        if(leapyr(endyear)) then
          days(2)=29
        else 
          days(2)=28
        endif

        DONE = .FALSE.
        ii   = 0

        tmpyear  = year
        tmpmonth = month
        tmpday   = day
        tmphour  = hour

C calculate the number of recs
        do 
          if(DONE) exit

          call get_next_time_step(tmpyear,tmpmonth,tmpday,tmphour,
     &                       tmpjday, dt)
          ii = ii + 1
          if(tmpyear .eq. endyear) then
            if(tmpmonth .eq. endmonth) then
              if(tmpday .eq. endday) then
                if(tmphour .eq. endhour) then
                  ii = ii + 1
                  DONE = .TRUE.
                endif
              endif
            endif
          endif

        enddo
        nrecs = ii;
      
      else
        offset = 0
        tmphour = hour
        do
          if(tmphour .eq. 0) exit
          tmphour = tmphour + dt
          offset = offset + 1
          if(tmphour .ge. 24) then
            tmphour = 0
          endif
        enddo

      endif

      if(nrecs .gt. MAXRECS) then
        write(*,*) 'ERROR: nrecs (',nrecs,') > MAXRECS (',MAXRECS,')'
        stop
      endif

c create the date structure for each timestep
      jday = day
      if(leapyr(year)) then
        days(2) = 29
      else 
        days(2) = 28
      endif

      do ii=1,month
        jday = jday + days(ii)
      enddo
  
      DONE = .FALSE.
      ii = 1

      do
        if(DONE) exit 

        dmy_hour(ii) = hour
        dmy_day(ii) = day
        dmy_month(ii) = month
        dmy_year(ii) = year
        dmy_day_in_year(ii) = jday

        call get_next_time_step(year,month,day,hour,jday,dt)

        ii = ii + 1
        if(ii .gt. nrecs) DONE=.TRUE.

      enddo

C write out the data
      if(DEBUG) then
        write(*,'(a)') 'Date structure:'
        do ii=1,nrecs
          write(*,'(i4,x,i2,x,i2,x,i2,x,i3)') 
     &          dmy_year(ii),dmy_month(ii),dmy_day(ii),
     &          dmy_hour(ii),dmy_day_in_year(ii)
        enddo
      endif
  
      end
      
C-----------------------------------------------------------------------
C leapyr: 	returns TRUE if the year is a leap year and FALSE otherwise
C-----------------------------------------------------------------------

      function leapyr (year)

      implicit none

C input variables
      integer year

C function return variable
      logical leapyr

      if(mod(year,400).eq.0) then
        leapyr = .TRUE.
      elseif((mod(year,4).eq.0).and.(mod(year,100).ne.0)) then
        leapyr = .TRUE.
      else 
        leapyr = .FALSE.
      endif

      end

C-----------------------------------------------------------------------
C get_next_time_step: 	returns TRUE if the year is a leap year and FALSE otherwise
C-----------------------------------------------------------------------

      subroutine get_next_time_step(year, month, day, hour,
     &                            jday, dt)

      implicit none

C input variables
      integer dt

C input/output variables
      integer hour, day, month, year, jday

C local variables
      logical leapyr
      integer days(12)
      data days /31,28,31,30,31,30,31,31,30,31,30,31/

C go to the next timestep  
      hour = hour + dt

C check for next day
      if(hour .ge. 24) then
        hour=0
        day = day + 1
        jday = jday + 1
    
        if(leapyr(year)) then
          days(2) = 29
        else 
          days(2) = 28
        endif    

C check for next month
        if(day .gt. days(month)) then
          day = 1
          month = month + 1

C check for next year
          if(month .eq. 13) then
	    month = 1
	    jday  = 1
	    year = year + 1
          endif
        endif
      endif

      end



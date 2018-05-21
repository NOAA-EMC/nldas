      function daybefore (today)

        implicit none
        character*8 today, daybefore
        integer year, month, day
        logical leapyr

        read(today, '(I4,I2,I2)') year, month, day
        day = day - 1
        if(day .eq. 0) then
          month = month - 1
          if(month .eq. 0) then
            year = year - 1
            month = 12
          endif
          if(month .eq. 1) day = 31
          if(month .eq. 2) then
            if(leapyr(year)) then
              day = 29
            else 
              day = 28
            endif
          endif
          if(month .eq. 3) day = 31
          if(month .eq. 4) day = 30
          if(month .eq. 5) day = 31
          if(month .eq. 6) day = 30
          if(month .eq. 7) day = 31
          if(month .eq. 8) day = 31
          if(month .eq. 9) day = 30
          if(month .eq. 10) day = 31
          if(month .eq. 11) day = 30
          if(month .eq. 12) day = 31
        endif
        write(daybefore, '(I4.4,I2.2,I2.2)') year, month, day

      end function daybefore



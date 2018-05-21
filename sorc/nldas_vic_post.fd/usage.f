C-----------------------------------------------------------------------
C usage: 	shows a program usage message and exitss
C-----------------------------------------------------------------------

      subroutine usage()
      
      character*40 str

      call getarg(0, str)
      write(*,*) 'Usage: ', str, ' <model_control_file>'
      stop 
      
      end



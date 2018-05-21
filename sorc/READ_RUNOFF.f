      SUBROUTINE READ_RUNOFF(nldas,VAR1,filname,lugb)

      implicit none

      integer nldas, lugb
      real VAR1(nldas)

      character filname*100
      
      OPEN(lugb, FILE=filname, FORM="UNFORMATTED")
      READ(lugb) VAR1

      END    

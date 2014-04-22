      SUBROUTINE READ_INT_SPATIAL(fn,nx,ny,nz,temp)

      IMPLICIT NONE

      INTEGER nx, ny, nz, i, j, k
      INTEGER temp(nx,ny,nz)
      CHARACTER*100 fn
      INTEGER ios
      
      OPEN(80,FILE = TRIM(fn),IOSTAT=ios,
     &     STATUS = 'OLD', FORM='unformatted')
      IF (ios .NE. 0) THEN
         WRITE(*,*) 'File '//TRIM(fn)//
     &        ' does not exist, Subroutine READ_INT_SPATIAL'
CJZ  Added Error Checking
         CALL ERREXIT(ios)
         STOP
      END IF
      READ(80) (((temp(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      CLOSE(80)

C      WRITE(*,*) (((temp(i,j,k),i=1,4),j=1,10),k=1,nz)

      END

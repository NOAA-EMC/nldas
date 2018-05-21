C-----------------------------------------------------------------------
C read_bin: 	reads data from VIC binary output files and converts
C		units as necessary
C-----------------------------------------------------------------------

      subroutine read_bin(funit,nx,ny,mask,pmask,pmask_value,
     &                    rec,data,mult,offset,rate,dt,nrecs)
      
      implicit none

      include "postproc.h"

C input args
      integer funit, nx, ny, rec, dt, nrecs
      integer pmask(nx,ny), pmask_value
      logical*1 mask(nx,ny)
      real    mult, offset
      logical*1 rate

C input/output args
      real    data(nx,ny)

C local args
      integer i, j, record, cell 

      if(DEBUG) write(*,'(a)') 'read_bin'

C read in the data (for mask points only)
      record=rec
      cell = 0
      do j=1,ny
        do i=1,nx
          cell = cell + 1
          if(mask(i,j) .and. pmask(i,j).eq.pmask_value) then
            read(funit, rec=record, err=999) data(i,j)
            if(DEBUG) 
     &        write(*,*) cell, rec, record, i, j, data(i,j)

C convert units as required
            if(mult .ne. -999) data(i,j) = data(i,j)*mult
            if(offset .ne. -999) data(i,j) = data(i,j)+offset
            if(rate) data(i,j) = data(i,j)/(dt*3600.0)
c            if(DEBUG)
c     &        write(*,*) data(i,j)

C skip to the next cell
            record=record+nrecs

          else
            data(i,j) = UNDEF
          endif
        enddo
      enddo

      return

  999 write(*,*) 'error in read_bin: ',record,rec
      stop
      
      end
      

C-----------------------------------------------------------------------
C get_mask: 	reads in the spatial mask from the arc ascii grid file
C		and calculates the number of active cells
C-----------------------------------------------------------------------

      subroutine get_mask(nx, ny, mask_file, pmask_file,
     &                            mask, pmask, pmask_value)

      implicit none

      include 'postproc.h'

C inputs variables
      integer nx, ny, pmask_value
      character*256 mask_file, pmask_file

C input/output variables
      integer pmask(nx,ny)
      logical*1 mask(nx,ny)

C local variables
      integer ncols, nrows, ncells, nodata
      real xllcorner, yllcorner, cellsize
      character*40 str
      integer i,j, status
      integer data(nx,ny)

C the start
      if(DEBUG) write(*,'(a)') 'get_mask'

C open the mask file
cccc      open(unit=10,file=mask_file, status='old', iostat=status)
cccc      if( status .gt. 0 ) then
cccc        write(*,*) 'ERROR: Can''t open mask file: ',trim(mask_file)
cccc        stop
cccc      endif
cccc      if(VERBOSE) then
cccc        write(*,'(a,a)') 'Opened mask file: ',trim(mask_file)
cccc      endif

C read the header data
cccc      if(DEBUG) write(*,'(a)') 'Mask header:'
cccc      read(10,*) str, ncols 
cccc      if(DEBUG) write(*,10) str, ncols
cccc      read(10,*) str, nrows 
cccc      if(DEBUG) write(*,10) str, nrows
cccc      read(10,*) str, xllcorner 
cccc      if(DEBUG) write(*,11) str, xllcorner
cccc      read(10,*) str, yllcorner 
cccc      if(DEBUG) write(*,11) str, yllcorner
cccc      read(10,*) str, cellsize 
cccc      if(DEBUG) write(*,11) str, cellsize
cccc      read(10,*) str, nodata 
cccc      if(DEBUG) write(*,10) str, nodata
   10 format(a20,i20)
   11 format(a20,f20.5)

C check that our array sizes are ok
cccc      if(nx .ne. ncols .or. ny .ne. nrows) then
cccc        write(*,*) 'ERROR: mask nx (',nx,') != ncols (',ncols,')'
cccc        write(*,*) 'ERROR: mask ny (',ny,') != nrows (',nrows,')'
cccc       stop
cccc      endif

C read the mask data
cccc      ncells = 0
cccc      do j=ny,1,-1
cccc        read(10,*) (data(i,j), i=1,nx)
cccc      enddo
cccc      do j=ny,1,-1
cccc        do i=1,nx
cccc          if(data(i,j) .gt. 0) then
cccc            mask(i,j)=.TRUE.
cccc           ncells = ncells + 1
cccc          else
cccc            mask(i,j)=.FALSE.
cccc          endif
cccc        enddo
cccc      enddo
      
C close the mask file
cccc      close (10)

C write the mask data
cccc      if(DEBUG) then
cccc        write(*,'(a)') 'Mask array:'
cccc        do j=ny,1,-1
cccc          write(*,*) (mask(i,j), i=1,nx)
cccc        enddo
cccc      endif
cccc      if(VERBOSE) write(*,'(a,i)') 'mask ncells = ',ncells

C open the processor mask file
      open(unit=10,file=pmask_file, status='old', iostat=status)
      if( status .gt. 0 ) then
        write(*,*) 
     &  'ERROR: Can''t open processor mask file: ',trim(pmask_file)
        stop
      endif
      if(VERBOSE) then
        write(*,'(a,a)') 'Opened processor mask file: ',trim(pmask_file)
      endif

C read the header data
      if(DEBUG) write(*,'(a)') 'Processor mask header:'
      read(10,*) str, ncols 
      if(DEBUG) write(*,10) str, ncols
      read(10,*) str, nrows 
      if(DEBUG) write(*,10) str, nrows
      read(10,*) str, xllcorner 
      if(DEBUG) write(*,11) str, xllcorner
      read(10,*) str, yllcorner 
      if(DEBUG) write(*,11) str, yllcorner
      read(10,*) str, cellsize 
      if(DEBUG) write(*,11) str, cellsize
      read(10,*) str, nodata 
      if(DEBUG) write(*,10) str, nodata

C check that our array sizes are ok
      if(nx .ne. ncols .or. ny .ne. nrows) then
        write(*,*) 'ERROR: pmask nx (',nx,') != ncols (',ncols,')'
        write(*,*) 'ERROR: pmask ny (',ny,') != nrows (',nrows,')'
        stop
      endif

C read the processor mask data
      do j=ny,1,-1
        read(10,*) (pmask(i,j), i=1,nx) 
C also store as logical mask
        do i=1,nx
          if(pmask(i,j).eq.pmask_value) then
            mask(i,j) = .TRUE.
          else
            mask(i,j) = .FALSE.
          endif
        enddo
      enddo

C close the processor mask file
      close (10)

C write the processor mask data
      if(DEBUG) then
        write(*,'(a)') 'Processor mask array:'
        do j=ny,1,-1
          write(*,*) (pmask(i,j), i=1,nx)
        enddo
      endif

C get total number of active cells
      ncells=0
      do j=1,ny
        do i=1,nx
C note where the mask and pmask don't match
          if(mask(i,j) .and. (pmask(i,j).eq.0)) then
            write(*,*) 'NO mask match ',i,j,mask(i,j),pmask(i,j)
          endif
          if(pmask(i,j).eq.pmask_value) then
            ncells=ncells+1
          endif
        enddo
      enddo
      if(VERBOSE) write(*,'(a,i)') 'processor mask ncells = ',ncells

      end

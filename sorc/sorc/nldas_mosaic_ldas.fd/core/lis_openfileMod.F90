!-------------------------------------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
! Released October 2005
!
! See SOFTWARE DISTRIBUTION POLICY for software distribution policies
!
! The LIS source code and documentation are in the public domain,
! available without fee for educational, research, non-commercial and
! commercial purposes.  Users may distribute the binary or source
! code to third parties provided this statement appears on all copies and
! that no charge is made for such copies.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See COPYRIGHT.TXT for copyright details.
!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: lis_openfileMod.F90
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines for opening data files.
!   
! !REVISION HISTORY: 
!  08Apr04    James Geiger Initial Specification
! 
!EOP
module lis_openfileMod
  use lisdrv_module, only : lis
  use lis_indices_module 

#if ( defined OPENDAP )
   logical, parameter :: use_opendap_server = .true.
#else
   logical, parameter :: use_opendap_server = .false.
#endif

contains
!BOP
! !ROUTINE: lis_set_filename 
! This routine overwrites the path for a GDS run
! 
! !INTERFACE:
subroutine lis_set_filename(file,time_offset,prefix_only)
#if ( defined OPENDAP )
   use opendap_module, only : opendap_data_prefix, ciam 
#endif
  character(len=*), intent(inout) :: file
  character(len=*), optional      :: time_offset
  logical, optional               :: prefix_only
#if ( defined OPENDAP )
  if ( use_opendap_server ) then
     if ( present(prefix_only) ) then
       if ( prefix_only ) then
        file = trim(opendap_data_prefix)//'/'//trim(adjustl(ciam))//'/'// &
               trim(adjustl(file))
        endif
     elseif ( PRESENT(time_offset) ) then
        file = trim(opendap_data_prefix)//'/'// &
             trim(adjustl(ciam))//'/'//"var_"//time_offset//".bin"
     else
        file = trim(opendap_data_prefix)//'/'// &
             trim(adjustl(ciam))//'/'//"var.bin"
     endif
  endif
#endif
end subroutine lis_set_filename
!BOP
! !ROUTINE: lis_open_file 
!
! !DESCRIPTION: 
! This routine is a generic open routine.  It parses its optional input
! arguments and builds an approriate open call.  It also determines
! whether or not data must be retrieve via a GraDS-DODS data server (GDS).
! If so, it calls the specified GDS script.
!
! !INTERFACE: 
subroutine lis_open_file(unit, file, form, status, access, recl, script, time_offset)

   implicit none

!INPUT PARAMETERS:
   integer,          intent(in) :: unit
   character(len=*), intent(in) :: file
   character(len=*), optional   :: form
   character(len=*), optional   :: status
   character(len=*), optional   :: access
   integer,          optional   :: recl
   character(len=*), optional   :: script
   character(len=*), optional   :: time_offset

!LOCAL VARIABLES:
   integer                      :: ios
   character(len=11)            :: form_use
   character(len=7)             :: status_use
   character(len=10)            :: access_use
   character(len=15)            :: script_use
   character(len=4)             :: cunit
   character(len=80)            :: file_tmp
!EOP

   ! Do not over-write the name of the file to open.  This original name
   ! may be used elsewhere in the code; e.g. lis%p%elevfile.
   ! Do all file name process on a copy.
   file_tmp = file

   ! If optional values are not assigned by caller, then set default values.
   if ( .not. PRESENT(form) ) then
      form_use ='unformatted'
   elseif ( trim(adjustl(form)) == 'unformatted' .or. &
            trim(adjustl(form)) == 'formatted'          ) then
      form_use = trim(adjustl(form))
   endif

   if ( .not. PRESENT(status) ) then
      status_use = 'old'
   elseif ( trim(adjustl(status)) == 'old'     .or. &
            trim(adjustl(status)) == 'new'     .or. &
            trim(adjustl(status)) == 'replace' .or. &
            trim(adjustl(status)) == 'unknown'        ) then
      status_use = trim(adjustl(status))
   endif
   if ( .not. PRESENT(access) ) then
      if(lis%d%gridDesc(9) .eq. 0.01) then 
         access_use = 'direct'
      else
         access_use = 'sequential'
      endif
   elseif ( trim(adjustl(access)) == 'sequential' .or. &
            trim(adjustl(access)) == 'direct'            ) then
      access_use = trim(adjustl(access))
   endif
!   if ( .not. PRESENT(recl) ) then
!      recl = 4
!   endif
   if ( .not. PRESENT(script) ) then
      script_use = 'none'
   else
      script_use = trim(adjustl(script))
   endif

   ! If script exists, retrieve data through GrADS-DODS server
   ! (if necessary)
   if ( use_opendap_server ) then
      if ( script_use /= 'none' ) then
         if(.not.PRESENT(time_offset)) then 
            call retrieve_data(file_tmp, script_use)
         else
            call retrieve_data(file_tmp, script_use, time_offset)
         endif
      endif
   endif
   ! Open the file
   call lis_log_msg('MSG: lis_open_file -- Opening '//trim(file_tmp))
   if ( access_use == 'sequential' ) then
      open(unit=unit, file=file_tmp, form=form_use, status=status_use, &
           access=access_use, IOSTAT=ios)
   else
      open(unit=unit, file=file_tmp, form=form_use, status=status_use, &
           access=access_use, recl=recl, IOSTAT=ios)
   endif

   ! Check the status of the open call
   write(cunit,'(i4)') unit
   if ( ios /= 0 ) then
      call lis_log_msg('ERR: lis_open_file -- Cannot open file '&
                       //trim(file_tmp)//' on unit '//adjustl(cunit))
      call endrun
   else
      call lis_log_msg('MSG: lis_open_file -- Successfully opened '&
                       //trim(file_tmp)//' on unit '//adjustl(cunit))
   endif
 
   return

end subroutine lis_open_file
!BOP
! !ROUTINE: lis_read_file 
!
! !DESCRIPTION: 
! This routine is a generic read routine.  It parses its optional input
! arguments and builds an approriate read call. 
!
! !INTERFACE: 
subroutine lis_read_file(unit, array)
  implicit none
  integer,          intent(in) :: unit
  real,          intent(inout) :: array(lis_nc_data, lis_nr_data)

  integer :: line1, line2, line
  integer :: c,r, glnc, glnr

#if ( defined OPENDAP )
  line = 1
  do r=1,lis_nr_data
     do c=1,lis_nc_data
        read(unit,rec=line) array(c,r)
        line = line + 1
     enddo
  enddo
#else
     line1 = nint((lis%d%gridDesc(4)-lis%d%gridDesc(44))/lis%d%gridDesc(9))+1
     line2 = nint((lis%d%gridDesc(5)-lis%d%gridDesc(45))/lis%d%gridDesc(10))+1
     do r=1,lis%d%lnr
        do c=1,lis%d%lnc
           glnc = line2+c-1
           glnr = line1+r-1
           line = (glnr-1)*nint(lis%d%gridDesc(42))+glnc
           read(unit,rec=line) array(c,r)
        enddo
     enddo
#endif
end subroutine lis_read_file

!BOP
! !ROUTINE: retrieve_data 
!
! !DESCRIPTION: 
! This routine retrieves data from a GDS.  It will make 3 attempts to
! retrieve data from the server.  If the data cannot be retrieved, this
! routine aborts by calling endrun.
!
! !INTERFACE: 
subroutine retrieve_data(file, script, time_offset)
!EOP

   character(len=*), intent(inout) :: file
   character(len=*), intent(in)    :: script
   character(len=*), optional      :: time_offset

#if ( defined OPENDAP )

   logical :: exists
   integer :: try
   character(len=80)               :: file_tmp

   exists = .false.
   try = 1

   do
      if ( .not. exists .and. try < 4 ) then ! keep trying to retrieve file

         file_tmp = trim(file)
         if(.not.PRESENT(time_offset)) then 
            call retrieve_script(file_tmp,script)
         else
            call retrieve_script(file_tmp, script, time_offset)
         endif

         inquire(FILE=file_tmp, EXIST=exists)
         try = try + 1
      else
         if ( .not. exists ) then ! error, could not retrieve the file
            call lis_log_msg('ERR: lis_open_file -- '// &
                             'Could not retrieve data file '//trim(file_tmp))
            call endrun
         else ! got it, break the do loop
            exit
         endif
      endif
   enddo

   file = trim(file_tmp)
#endif

   return

end subroutine retrieve_data

!BOP
! !ROUTINE: retrieve_script 
!
! !DESCRIPTION: 
! This routine makes the system call that executes the
! GrADS script that retrieves data from a GDS.
!
! !INTERFACE: 
subroutine retrieve_script(file, script, time_offset)
!EOP
#if ( defined OPENDAP )
   use lisdrv_module,        only : lis
   use geosdomain_module,    only : geosdrv
   use opendap_module,       only : ciam, cdom,             &
                                    cparm_slat, cparm_nlat, &
                                    cparm_wlon, cparm_elon
   use agrmetopendap_module, only : agrmet_time_index
   use geosopendap_module,   only : cgeos_slat, cgeos_nlat, get_geos_index
   use gdasopendap_module,   only : cgdas_slat, cgdas_nlat, &
                                    cgdas_wlon, cgdas_elon

   character(len=*), intent(inout) :: file
   character(len=*), intent(in)    :: script
   character(len=*), optional      :: time_offset
   character(len=4)                :: ctime_index
   character(len=16)               :: tmp_script

   call lis_log_msg('MSG: lis_open_file -- Retrieving data via '// &
                    trim(script)//' script')

   if ( trim(adjustl(script)) == "getgeos.pl" ) then

      if ( time_offset == '1' ) then
         write(ctime_index, '(i4)') get_geos_index(0)
      else ! order == 2
         write(ctime_index, '(i4)') get_geos_index(3)
      endif

      call lis_set_filename(file,prefix_only=.true.)

      call system("opendap_scripts/"//trim(adjustl(script))//" "// &
                   ciam//" "//                                     & 
                   trim(file)//" "//                               &
                   ctime_index//" "//cgeos_slat//" "//cgeos_nlat)

   elseif ( trim(adjustl(script)) == "getagrmet_lw.pl" .or. &
            trim(adjustl(script)) == "getagrmet_sw.pl" ) then

      write(ctime_index, '(i4)') agrmet_time_index

      call lis_set_filename(file,prefix_only=.true.)

      call system("opendap_scripts/"//trim(adjustl(script))//" "// &
                  ciam//" "//                                      &
                  trim(file)//" "//                                &
                  ctime_index)

   elseif ( trim(adjustl(script)) == "getgdas.pl" ) then


      call lis_set_filename(file,prefix_only=.true.)

      call system("opendap_scripts/"//trim(adjustl(script))//" "// &
                  ciam//" "//                                      & 
                  trim(file)//" "//                                & 
                  cgdas_slat//" "//cgdas_nlat//" "//               &
                  cgdas_wlon//" "//cgdas_elon)


   elseif ( trim(adjustl(script)) == "getelev.pl" ) then

      select case ( lis%f%force )
         case ( 1 )
            tmp_script = 'getgdas_diff.pl'
         case ( 2 )
            ! When geosdrv%gridchange == .true. it means that LIS must make
            ! the GEOS3 to GEOS4 grid change when the time is right.
            ! When geosdrv%gridchange == .false. it means that LIS has made
            ! the GEOS3 to GEOS4 grid change.
            ! Thus when geosdrv%gridchange == .true., LIS is still using the
            ! GEOS3 elevation difference data; when it is .false., LIS is
            ! switching to the GEOS4 elevation difference data.
            if ( geosdrv%gridchange ) then 
               tmp_script = 'getgeos3_diff.pl'
            else
               tmp_script = 'getgeos4_diff.pl'
            endif
         case default
            call lis_log_msg('ERR: retrieve_script -- Cannot determine '// &
                             'necessary elevation difference file')
            call endrun
      endselect

      call lis_set_filename(file)
      call system("opendap_scripts/"//trim(adjustl(tmp_script))//" "// &
                  ciam//" "//                                          &
                  trim(file)//" "//cdom//" "//                         &
                  cparm_slat//" "//cparm_nlat//" "//                   &
                  cparm_wlon//" "//cparm_elon)

   elseif(PRESENT(time_offset)) then 

      call lis_set_filename(file,time_offset=time_offset)
      call system("opendap_scripts/"//trim(adjustl(script))//" "// &
                  ciam//" "//                                      &
                  trim(file)//" "//cdom//" "//time_offset//" "//   & 
                  cparm_slat//" "//cparm_nlat//" "//               &
                  cparm_wlon//" "//cparm_elon)
   else

      call lis_set_filename(file)
      call system("opendap_scripts/"//trim(adjustl(script))//" "// &
                  ciam//" "//                                      &
                  trim(file)//" "//cdom//" "//                     &
                  cparm_slat//" "//cparm_nlat//" "//               &
                  cparm_wlon//" "//cparm_elon)
   endif
#endif
   return
end subroutine retrieve_script

!BOP
!
! !ROUTINE: create_output_directory
!
! !DESCRIPTION:  
!  Create the output directory for the output data files.
!
! !REVISION HISTORY:
! 06 Jun 2005: James Geiger; Initial Specification
! 
! !INTERFACE:
subroutine create_output_directory(dir_name)
! !USES:
   use lisdrv_module, only : lis
  
   implicit none 
! !ARGUMENTS:
   character(len=40), optional :: dir_name
!EOP
   character(len=10) :: cdate
   character(len=4)  :: mname
   character(len=40) :: out_dname

   out_dname = trim(lis%o%odir)//'/EXP'//trim(adjustl(lis%o%expcode))//'/'

   select case ( lis%d%lsm )
      case ( 1 )
         mname = 'NOAH'
      case ( 2 )
         mname = 'CLM2'
      case ( 3 )
         mname = 'VIC'
      case default
         call lis_log_msg('ERR: create_output_directory -- '// &
                          'Unrecognized lis%d%lsm value')
         call endrun 
   endselect
   out_dname = trim(out_dname)//mname//'/'
   
   write(unit=cdate, fmt='(i4.4)') lis%t%yr
   out_dname = trim(out_dname)//trim(cdate)//'/'

   write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') lis%t%yr, lis%t%mo, lis%t%da
   out_dname = trim(out_dname)//trim(cdate)

   if ( present(dir_name) ) then
      dir_name = trim(out_dname)
   else
      call lis_log_msg('MSG: create_output_directory -- mkdir -p '// &
                       trim(out_dname))
      call system('mkdir -p '//trim(out_dname))
   endif

end subroutine create_output_directory
  

!BOP
!
! !ROUTINE: create_output_filename
!
! !DESCRIPTION:  
!  Create the file name for the output data files.
!
! !REVISION HISTORY:
! 06 Jun 2005: James Geiger; Initial Specification
! 
! !INTERFACE:
subroutine create_output_filename(fname, model_name, writeint)
! !USES:
   use lisdrv_module, only : lis
  
   implicit none 
! !ARGUMENTS:
   character(len=*), intent(out)          :: fname
   character(len=*), intent(in), optional :: model_name ! needed for gswp run
   real, intent(in), optional             :: writeint ! output writing interval
                                                      ! e.g., noahdrv%writeintn
                                                      ! this value is needed if
                                                      ! you are doing a gswp run
!EOP
   character(len=8)        :: date
   character(len=10)       :: time
   character(len=5)        :: zone
   integer, dimension(8)   :: values
 
   character(len=10)       :: cdate

   integer                 :: curr_mo = 0
   character(len=40)       :: dname
   character(len=80), save :: out_fname

   if ( lis%o%wout == 4 ) then  ! GSWP-2 style output
      if ( .not. present(model_name) ) then
         call lis_log_msg('ERR: create_output_filename -- You must supply '// &
                          'the model name.')
         call endrun
      endif
      if ( .not. present(writeint) ) then
         call lis_log_msg('ERR: create_output_filename -- You must supply '// &
                          'the write interval.')
         call endrun
      endif

      if ( curr_mo /= lis%t%mo ) then 
         curr_mo = lis%t%mo
   
         call date_and_time(date, time, zone, values)
   
         out_fname = trim(model_name)//'LIS_exp'//  &
                     trim(adjustl(lis%o%expcode))// &
                     '_'//date
   
         if ( writeint == 24 ) then 
            out_fname = trim(out_fname)//'_d'
         else
            out_fname = trim(out_fname)//'_h'
         endif
   
         write(unit=cdate, fmt='(i4.4, i2.2)') lis%t%yr, curr_mo
         out_fname = trim(out_fname)//trim(cdate)
         out_fname = trim(out_fname)//'.nc'
      endif
   else
      call create_output_directory(dir_name=dname)
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2, i2.2)') lis%t%yr, lis%t%mo, &
                                                        lis%t%da, lis%t%hr
      out_fname = trim(dname)//'/'//cdate

      select case ( lis%o%wout )
         case ( 1 )
            out_fname = trim(out_fname)//'.gs4r'
         case ( 2 )
            out_fname = trim(out_fname)//'.grb'
         case ( 3 )
            out_fname = trim(out_fname)//'.nc'
         case default
            call lis_log_msg('ERR: create_output_filename -- '// &
                             'Unrecognized lis%o%wout value')
            call endrun 
      endselect

   endif

   fname = trim(out_fname)
        
end subroutine create_output_filename

!BOP
!
! !ROUTINE: create_restart_filename
!
! !DESCRIPTION:  
!  Create the file name for the restart data files.
!
! !REVISION HISTORY:
! 29 Jun 2005: James Geiger; Initial Specification
! 
! !INTERFACE:
subroutine create_restart_filename(fname,extn)
! !USES:
   use lisdrv_module, only : lis
  
   implicit none 
! !ARGUMENTS:
   character(len=*), intent(out) :: fname
   character(len=*), intent(in)  :: extn

!EOP
   character(len=10)       :: cdate

   character(len=40)       :: dname
   character(len=80)       :: out_fname

   call create_output_directory(dir_name=dname)
   write(unit=cdate, fmt='(i4.4, i2.2, i2.2, i2.2)') lis%t%yr, lis%t%mo, &
                                                     lis%t%da, lis%t%hr
   out_fname = trim(dname)//'/LIS.E'//trim(adjustl(lis%o%expcode))// &
               '.'//cdate//trim(extn)

   fname = trim(out_fname)
        
end subroutine create_restart_filename

!BOP
!
! !ROUTINE: create_stats_filename
!
! !DESCRIPTION:  
!  Create the file name for the stats files.
!
! !REVISION HISTORY:
! 29 Jun 2005: James Geiger; Initial Specification
! 
! !INTERFACE:
subroutine create_stats_filename(fname, name)
! !USES:
   use lisdrv_module, only : lis
  
   implicit none 
! !ARGUMENTS:
   character(len=*), intent(out) :: fname
   character(len=*), intent(in)  :: name

!EOP
   character(len=80) :: out_fname

   out_fname = trim(lis%o%odir)//'/EXP'//trim(lis%o%expcode)//'/'//trim(name)

   fname = trim(out_fname)
        
end subroutine create_stats_filename

end module lis_openfileMod

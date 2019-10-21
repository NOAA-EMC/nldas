      SUBROUTINE READ_FORCING(nldas,f,tair,spfh,psfc,uwind,vwind,
     &   lwrd,edasprec,cape,pevap,prec,swrd,lb,filnam,lugb)
C__________________________________________________________________
C  abstract:  Sample code for reading hourly LDAS forcing file in
C  GRIB format, unpacking, and storing parameters (13 forcing grids)
C  in separate binary arrays.  Intended for use in LDAS drivers.       
C
C  author:  Marshall
C  org:  NCEP/EMC/MMB
C  date: 12 March 1999
C
C  modified by D. Lohmann for 2-D NOAH input
C  modified by Y. Xia for new NARR NLDAS forcing dataset 

        implicit none

        integer i, n, j, lugi, jpds(25), jgds(22), kpds(25), kgds(22),
     +          iret, jret, kf, k, nldas, lugb, parmid(11) 

        real f(nldas), tair(nldas), spfh(nldas), psfc(nldas), 
     +       uwind(nldas), vwind(nldas),lwrd(nldas),pevap(nldas),
     +       cape(nldas),prec(nldas),swrd(nldas),edasprec(nldas)
        character filnam*150
        character clugb*2

        logical*1 lb(nldas)

C  CONSTANTS, DATA STATEMENTS, ETC...

        data parmid /11,51,1,33,34,205,153,157,228,61,204/

C  OPEN INPUT FILE AND SEARCH FOR GRIB MESSAGE SPECIFIED BY
C  jpds(5) AND jpds(2)
c  when read grib file not owned by you use baopenr
c       call baopen (lugb, filnam, jret)
        call baopenr (lugb, filnam, jret)
        do i = 1, 11
c  Dag Lohmann, April 2002
C  j has to be set to -i to jump -1-j fields in the input file
          j = -i
          lugi = 0 
          jpds = -1
          jpds(5) = parmid(i)

          call getgb(lugb, lugi, nldas, j, jpds, jgds, kf, k,
     +                kpds, kgds, lb, f, iret)
c         print*, 'i= ', i, ' iret= ', iret, ' kpds(5) = ', kpds(5),
c    +            ' jret = ', jret, ' lugb= ', lugb
C          WRITE(*,*) ' kpds = ', kpds, ' jpds = ', jpds

          IF (iret .NE. 0) THEN
             OPEN(80, file='stop.file')
             WRITE(80,*) 'SAC model stopped in Reading Forcing'
             WRITE(*,*) 'IRET=', iret
             WRITE(*,*) 'STOP SAC MODEL in Reading Forcing'
             STOP
          END IF

C  COPY RETURNED DATA INTO ITS RESPECTIVE BINARY ARRAY.  NOTE
C  DATA ARE ASSIGNED TO PROPER ARRAY BASED ON LOOP INDEX (i).  

          if (i.eq.1) then
             do n = 1, nldas
                tair(n) = f(n)
             end do
          else if (i.eq.2) then
             do n = 1, nldas
                spfh(n) = f(n)
             end do
          else if (i.eq.3) then
             do n = 1, nldas
                psfc(n) = f(n)
             end do
          else if (i.eq.4) then
             do n = 1, nldas
                uwind(n) = f(n)
             end do
          else if (i.eq.5) then
             do n = 1, nldas
                vwind(n) = f(n)
             end do
          else if (i.eq.6) then
             do n = 1, nldas
                lwrd(n) = f(n)
             end do
          else if (i.eq.7) then
             do n = 1, nldas
                edasprec(n) = f(n)
             end do
          else if (i.eq.8) then
             do n = 1, nldas
                cape(n) = f(n)
             end do
          else if (i.eq.9) then
             do n = 1, nldas
                pevap(n) = f(n)
             end do
          else if (i.eq.10) then
           do n = 1, nldas
                   prec(n) = f(n)
                end do
          else if (i.eq.11) then
             do n = 1, nldas
                swrd(n) = f(n)
             end do
          end if
        end do
        call baclose(lugb,jret)

        end     


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
      SUBROUTINE PARTITION ( & 
                           NCH,ITYP,RZEXC, & 
                           RZEQ,VGWMAX,CDCR1,CDCR2, & 
                           PSIS,BEE,poros,WPWET, & 
                           ars1,ars2,ars3,ara1,ara2,ara3,ara4, & 
                           arw1,arw2,arw3,arw4,BUG, & 
                           srfexc,catdef, & 
                           AR1, AR2, AR4, srfmx, srfmn, & 
                           SWSRF1,SWSRF2,SWSRF4, & 
                           cor & 
                          )

      IMPLICIT NONE
!**** CHIP HEADER FILE
!****
      INTEGER  NTYPS, FRSTCH, MEMFAC
      INTEGER  NLAY, SFCLY, ROOTLY, RECHLY

      REAL  ZERO, ONE, PIE
      REAL ALHE, ALHS, ALHM, TF, STEFAN, RGAS, SHW, SHI, RHOW, GRAV
      REAL EPSILON, NOSNOW, ASUM

      PARAMETER (NTYPS = 10, FRSTCH = 1, MEMFAC = 5)
      PARAMETER (NLAY = 3)
      PARAMETER (SFCLY = 1, ROOTLY = SFCLY + 1, RECHLY = ROOTLY + 1)

      PARAMETER (ZERO = 0., ONE = 1., PIE = 3.14159265)
      PARAMETER (ALHE = 2.4548E6, ALHS = 2.8368E6, ALHM = ALHS-ALHE)
      PARAMETER (TF = 273.16)
      PARAMETER (STEFAN = 5.669E-8)
      PARAMETER (RGAS = .286*1003.5)
      PARAMETER (SHW = 4200., SHI = 2060.)
      PARAMETER (RHOW = 1000.)
      PARAMETER (GRAV = 9.81)
      PARAMETER (EPSILON = 18.01/28.97)
      PARAMETER (NOSNOW = 0.)

      INTEGER NCH,N
      INTEGER ITYP(NCH)
      REAL RZEXC(NCH),SRFEXC(NCH),CATDEF(NCH),AR1(NCH),AR2(NCH), &
          AR3(NCH),AR4(NCH), &
!c,VGPH1X(NCH),VGPH2X(NCH),
          RZEQ(NCH),VGWMAX(NCH),CDCR1(NCH),CDCR2(NCH)
      REAL SWSRF1(NCH),SWSRF2(NCH),SWSRF3(NCH),SWSRF4(NCH),cor(nch)

      REAL ars1(NCH),ars2(NCH),ars3(NCH),ara1(NCH),ara2(NCH),ara3(NCH), &
          ara4(NCH),arw1(NCH),arw2(NCH),arw3(NCH),arw4(NCH)

      REAL BEE(NCH),PSIS(NCH),WPWET(NCH),srfmx(nch),srfmn(nch), &
         poros(nch)
      REAL A150,W150,WMIN,AX,WMNEW,WRZ,TERM1,TERM2,TERM3,AREA0,AREA1, &
         AREA2,AREA3,AREA4,ASCALE,WILT,D1,D2,CDI,DELTA1,DELTA2, &
         DELTA4,MULTAR,CATDEFX,RZEQX,RZEQW,FACTOR,X0,RZEQY,CATDEFW, &
         AR1W
      LOGICAL*1 LSTRESS, BUG

      DATA LSTRESS/.FALSE./

      DO N=1,NCH

        WILT=WPWET(N)
        WRZ=RZEXC(N)/VGWMAX(N)
        CATDEFX=AMIN1( CATDEF(N) , CDCR1(N) )

! CDI DEFINES IF THE SHAPE PARAMETER IS ADJUSTED IN ONE OR TWO SEGMENTS
        if (ara1(n) .ne. ara3(n)) then
            cdi=(ara4(n)-ara2(n))/(ara1(n)-ara3(n))
          else
            cdi=0.
          endif

        AR1(N)= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFX) &
                /(1.+ars2(n)*CATDEFX+ars3(n)*CATDEFX*CATDEFX))) 

        if (CATDEFX .ge. cdi) then
            ax=ara3(n)*CATDEFX+ara4(n)
          else
            ax=ara1(n)*CATDEFX+ara2(n)
          endif

        WMIN=AMIN1(1.,AMAX1(0.,arw4(n)+(1.-arw4(n))* &
                (1.+arw1(n)*CATDEFX) &
                /(1.+arw2(n)*CATDEFX+arw3(n)*CATDEFX*CATDEFX)))

!**** CRITICAL VALUE 1: AVERAGE MOISTURE IN ROOT ZONE AT WMIN
!**** ASSOCIATED WITH CATDEF.
        RZEQX=(WMIN-1.-(2./AX))*EXP(-AX*(1.-WMIN)) + WMIN + (2./AX)

!**** CRITICAL VALUE 2: AVERAGE MOISTURE IN ROOT ZONE WHEN WMIN
!**** IS EXACTLY AT WILTING POINT.
        RZEQW=(WILT-1.-(2./AX))*EXP(-AX*(1.-WILT)) + WILT + (2./AX)

!**** SITUATION 1: CATDEF LE CDCR1
        IF(CATDEF(N) .LE. CDCR1(N)) THEN
          RZEQY=RZEQX+WRZ
          WMNEW=WMIN+WRZ
          AREA1=(1.+AX-AX*WMIN)*EXP(-AX*(1.-WMIN))
          AREA2=(1.+AX-AX*WMNEW)*EXP(-AX*(1.-WMNEW))
          IF(WMNEW .GE. WILT) THEN
            AR1(N)=AR1(N)+AREA2-AREA1
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF
          IF(WMNEW .LT. WILT) THEN
            AREA3=(1.+AX-AX*WILT)*EXP(-AX*(1.-WILT))
            AR1(N)=AR1(N)+AREA3-AREA1
            AR2(N)=1.-AR1(N)
            FACTOR=(RZEQX+WRZ-WILT)/(RZEQW-WILT)
            AR1(N)=AR1(N)*FACTOR
            AR2(N)=AR2(N)*FACTOR
            AR4(N)=1.-FACTOR
            ENDIF
          ENDIF

!**** SITUATION 2: CATDEF GT CDCR1
!        IF(CATDEF(N) .GT. CDCR1(N)) THEN
!          AR2(N)=1.-AR1(N)
!          FACTOR=(CDCR2(N)-CATDEF(N)-RZEXC(N))/(CDCR2(N)-CDCR1(N))
!          FACTOR=AMIN1(1., AMAX1(0., FACTOR))
!          X0=(RZEQW-WILT)/(RZEQX-WILT)
!          IF(FACTOR.LT.X0) THEN
!            AR1(N)=AR1(N)*FACTOR
!            AR2(N)=AR2(N)*FACTOR
!            AR4(N)=1.-FACTOR
!            ENDIF
!          ENDIF

!**** SITUATION 2: CATDEF GT CDCR1
        IF(CATDEF(N) .GT. CDCR1(N)) THEN
          FACTOR=(CDCR2(N)-CATDEF(N))/(CDCR2(N)-CDCR1(N))
          RZEQY=WILT+(RZEQX-WILT)*FACTOR+WRZ

          IF(RZEQY .LT. WILT) THEN
!           IF(RZEQY .LT. WILT-.001) THEN
!             WRITE(*,*) 'RZEXC(N) WAY TOO LOW!: N=',N,'RZEXC=',RZEXC(N)
!           ELSE
!             WRITE(*,*) 'RZEXC(N) TOO LOW: N=',N
!           END IF
            RZEQY=WILT

            ENDIF

          IF(RZEQY .GE. RZEQX) THEN  ! RZEXC BRINGS MOISTURE ABOVE CDCR1 POINT
            WMNEW=WMIN+(RZEQY-RZEQX)
            AREA1=(1.+AX-AX*WMIN)*EXP(-AX*(1.-WMIN))
            AREA2=(1.+AX-AX*WMNEW)*EXP(-AX*(1.-WMNEW))
            AR1(N)=AR1(N)+(AREA2-AREA1)
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF

          IF(RZEQY .LT. RZEQX .AND. RZEQY .GE. RZEQW) THEN
            CATDEFW=CDCR2(N)+ &
                ((RZEQW-WILT)/(RZEQX-WILT))*(CDCR1(N)-CDCR2(N))
            AR1W= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFW) &
                /(1.+ars2(n)*CATDEFW+ars3(n)*CATDEFW*CATDEFW)))
            FACTOR=(RZEQY-RZEQW)/(RZEQX-RZEQW)
            AR1(N)=AR1W+FACTOR*(AR1(N)-AR1W)
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF

          IF(RZEQY .LT. RZEQW) THEN
            CATDEFW=CDCR2(N)+ &
                ((RZEQW-WILT)/(RZEQX-WILT))*(CDCR1(N)-CDCR2(N))
            AR1W= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFW) &
                /(1.+ars2(n)*CATDEFW+ars3(n)*CATDEFW*CATDEFW)))
            AR1(N)=AR1W
            AR2(N)=1.-AR1(N)
            FACTOR=(RZEQY-WILT)/(RZEQW-WILT)
            AR1(N)=AR1(N)*FACTOR
            AR2(N)=AR2(N)*FACTOR
            AR4(N)=1.-FACTOR
            ENDIF

          ENDIF


        SWSRF1(N)=1.
        SWSRF2(N)=AMIN1(1., AMAX1(0., RZEQY))
        SWSRF4(N)=AMIN1(1., AMAX1(0., WILT))

!**** EXTRAPOLATION OF THE SURFACE WETNESSES

! 1st step: surface wetness in the unstressed fraction without considering
!           the surface excess; we just assume an equilibrium profile from 
!           the middle of the root zone to the surface.

        IF (SWSRF2(N) .NE. 0.) THEN 
          SWSRF2(N)=((SWSRF2(N)**(-BEE(N))) - (.5/PSIS(N))) &
                **(-1./BEE(N))
          ENDIF
! if swsrf2 is 0. in the RZ, it is also 0. at the surface (and it corresponds 
! to ar2=0.)

        SWSRF4(N)=AMIN1(1.,AMAX1(0.,SWSRF4(N)))
        IF (SWSRF4(N) .NE. 0.) THEN
           SWSRF4(N)=((SWSRF4(N)**(-BEE(N))) - (.5/PSIS(N))) &
                **(-1./BEE(N))
          ENDIF
! if swsrf4 is 0. in the RZ, it is also 0. at the surface(and it corresponds 
! to ar4=0.)

! srfmx is the maximum amount of water that can be added to the surface layer
! The choice of defining SWSRF4 like SWSRF2 needs to be better examined.
        srfmx(n)=ar2(n)*(1.-swsrf2(n))*(20.*poros(n))
        srfmx(n)=srfmx(n)+ar4(n)*(1.-swsrf4(n))*(20.*poros(n))
!**** For calculation of srfmn, assume surface moisture associated with
!**** AR1 is constantly replenished by water table.
        srfmn(n)=-(ar2(n)*swsrf2(n)+ar4(n)*swsrf4(n))*(20.*poros(n))

        if(srfexc(n).gt.srfmx(n)) then
          cor(n)=srfmx(n)-srfexc(n)
          catdef(n)=catdef(n)-cor(n)
          srfexc(n)=srfmx(n)
        else if(srfexc(n).lt.srfmn(n)) then
          cor(n)=srfexc(n)-srfmn(n)
          catdef(n)=catdef(n)-cor(n)
          srfexc(n)=srfmn(n)
        else
          cor(n)=0.
        endif
          
        SWSRF2(N)=SWSRF2(N)+SRFEXC(N)/(20.*poros(n)*(1.-ar1(n))+1.e-20)
        SWSRF2(N)=AMIN1(1., AMAX1(1.E-10, SWSRF2(N)))
        swsrf4(n)=swsrf4(n)+srfexc(n)/(20.*poros(n)*(1.-ar1(n))+1.e-20)
        SWSRF4(N)=AMIN1(1., AMAX1(wilt, SWSRF4(N)))
        if(swsrf4(n).eq.wilt) swsrf4(n)=1.e-10
        if(swsrf4(n).gt.wilt .and. swsrf4(n).lt.wilt+.005) &
            swsrf4(n)=1.e-10+(wilt+.005)*(swsrf4(n)-wilt)/.005  ! ramping

        IF (AR1(N) .ge. 1.-1.E-5) then
          AR1(N)=1.
          AR2(N)=0.
          AR4(N)=0.
          SWSRF2(N)=1.
          SWSRF4(N)=wilt
        ELSE IF (AR1(N) .lt. 0.) then
          if(ar1(n) .lt. -1.e-2) write(*,*) 'ar1(n) too low: n=',n, &
                                           'ar1=',ar1(n)
          AR1(N)=0.
          ENDIF

        AR2(N)=AMAX1(0., AMIN1(1., AR2(N)))
        AR4(N)=AMAX1(0., AMIN1(1., AR4(N)))
        ASUM=AR1(N)+AR2(N)+AR4(N)
        IF (ASUM .LT. .999 .OR. ASUM .GT. 1.001) then
            WRITE(*,*) 'Areas do not add to 1, sum=',ASUM
  
            write(*,*) ' ityp=',ityp(n),' rzexc=',rzexc(n),' rzeq=', &
                rzeq(n),' vgwmax=',vgwmax(n),' cdcr1=',cdcr1(n), &
                ' cdcr2=',cdcr2(n),' psis=',psis(n),' bee=',bee(n),&
                ' poros=',poros(n),' wpwet=',wpwet(n),' ars1=',ars1(n),&
                ' ars2=',ars2(n),' ars3=',ars3(n),' ara1=',ara1(n),&
                ' ara2=',ara2(n),' ara3=',ara3(n),' arw1=',arw1(n),&
                ' arw2=',arw2(n),' arw3=',arw3(n),' arw4=',arw4(n),&
                ' srfexc=',srfexc(n),' catdef=',catdef(n)

          end if
        ENDDO
           
      RETURN
      END







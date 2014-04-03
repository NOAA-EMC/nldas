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
      SUBROUTINE RZEQUIL ( &
                         NCH,ITYP,CATDEF,VGWMAX,CDCR1,CDCR2,WPWET, &
                         ars1,ars2,ars3,ara1,ara2,ara3,ara4,&
                         arw1,arw2,arw3,arw4,&
                         RZEQ &
                        )

      IMPLICIT NONE
      INTEGER NCH,N
      INTEGER ITYP(NCH)
      REAL CATDEF(NCH),VGWMAX(NCH),RZEQ(NCH),CDCR1(NCH),CDCR2(NCH)
      real ars1(NCH),ars2(NCH),ars3(NCH),ara1(NCH),ara2(NCH),ara3(NCH), &
          ara4(NCH),arw1(NCH),arw2(NCH),arw3(NCH),arw4(NCH),WPWET(NCH)
      REAL AX,WMIN,ASCALE,TERM1,TERM2,ONEMW,WT,AREA1,AREA1X,cdi,wilt, &
         catdefx,wmnew,factor

! ----------------------------------------------------------------------

      DO N=1,NCH

        WILT=WPWET(N)
        CATDEFX=AMIN1( CATDEF(N) , CDCR1(N) )

! CDI DEFINES IF THE SHAPE PARAMETER IS ADJUSTED IN ONE OR TWO SEGMENTS
        if (ara1(n) .ne. ara3(n)) then
            cdi=(ara4(n)-ara2(n))/(ara1(n)-ara3(n))
          else
            cdi=0.
          endif


        AREA1X= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFX) &
                /(1.+ars2(n)*CATDEFX+ars3(n)*CATDEFX*CATDEFX))) 
            
        if (CATDEFX .ge. cdi) then
            ax=ara3(n)*CATDEFX+ara4(n)
          else
            ax=ara1(n)*CATDEFX+ara2(n)
          endif

        WMIN=AMIN1(1.,AMAX1(0.,arw4(n)+(1.-arw4(n))* &
                (1.+arw1(n)*CATDEFX) &
                /(1.+arw2(n)*CATDEFX+arw3(n)*CATDEFX*CATDEFX)))


        WMNEW=WMIN

!**** FOR POSITIVE EXCESS, ALLOW SATURATED AREA TO INCREASE
!**** (BY AREA2-AREA1).  FOR NEGATIVE EXCESS, ASSUME THAT SATURATED
!**** FRACTION IS MAINTAINED BY WATER TABLE DISTRIBUTION --- DON'T
!**** REDUCE SATURATED AREA.

        AREA1=(1.+AX-AX*WMNEW)*EXP(-AX*(1.-WMNEW))

!**** CRITICAL VALUE 1: AVERAGE MOISTURE IN ROOT ZONE AT WMIN
!**** ASSOCIATED WITH CATDEF.
        RZEQ(N)=(WMIN-1.-(2./AX))*EXP(-AX*(1.-WMIN)) + WMIN + (2./AX)
!        RZEQW=(WILT-1.-(2./AX))*EXP(-AX*(1.-WILT)) + WILT + (2./AX)

        IF(CATDEF(N) .GT. CDCR1(N)) THEN
          FACTOR=(CDCR2(N)-CATDEF(N))/(CDCR2(N)-CDCR1(N))
          RZEQ(N)=WILT+(RZEQ(N)-WILT)*FACTOR
          ENDIF

! scaling:    
        RZEQ(N)=AMIN1(1.,AMAX1(0.,RZEQ(N)))
        RZEQ(N)=RZEQ(N)*VGWMAX(N)

      ENDDO

      RETURN
      END

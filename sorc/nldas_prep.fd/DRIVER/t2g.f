!=========================================================================
!
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
!  D                                      L  Land Modeling and Data 
!  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
!  S                                      A  This is the GSFC-LDAS Code.
!  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
!
!   GSFC - NCEP - OH - Princeton - Washington - Rutgers
!
!=========================================================================
! t2g.f and g2t.f: 
!
! DESCRIPTION:
!  Transfer variables between grid and tile space.
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial Code
!=========================================================================


      SUBROUTINE T2GI(T,G,NC,NR,NCH,FGRD,COL,ROW)
      IMPLICIT NONE
      INTEGER NCH,NC,NR,C,R,I,COL(NCH),ROW(NCH)
      INTEGER T(NCH),G(NC,NR)
      REAL FGRD(NCH)

      DO C=1,NC
       DO R=1,NR
        G(C,R)=0
       ENDDO
      ENDDO

      DO I=1,NCH
       G(COL(I),ROW(I))=G(COL(I),ROW(I))+T(I)*FGRD(I)
      ENDDO
      RETURN
      END

      SUBROUTINE T2GR(T,G,NC,NR,NCH,FGRD,COL,ROW)
      IMPLICIT NONE
      INTEGER NCH,NC,NR,C,R,I,COL(NCH),ROW(NCH)
      REAL T(NCH),G(NC,NR)
      REAL FGRD(NCH)

      DO C=1,NC
       DO R=1,NR
        G(C,R)=0.0
       ENDDO
      ENDDO

      DO I=1,NCH
       G(COL(I),ROW(I))=G(COL(I),ROW(I))+T(I)*FGRD(I)
      ENDDO
      RETURN
      END


      SUBROUTINE G2TI(T,G,NC,NR,NCH,COL,ROW)
      IMPLICIT NONE
      INTEGER NCH,NC,NR,I,COL(NCH),ROW(NCH)
      INTEGER T(NCH),G(NC,NR)

      DO I=1,NCH
       T(I)=G(COL(I),ROW(I))
      ENDDO
      RETURN
      END

      SUBROUTINE G2TR(T,G,NC,NR,NCH,COL,ROW)
      IMPLICIT NONE
      INTEGER NCH,NC,NR,I,COL(NCH),ROW(NCH)
      REAL T(NCH),G(NC,NR)

      DO I=1,NCH
       T(I)=G(COL(I),ROW(I))
      ENDDO
      RETURN
      END









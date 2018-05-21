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
!  This is a subroutine for the LDAS driver
!  i've commented out the program part of the this file.
!  because the driver only needs the subroutine.
!  -Jared
!
!
!      	This program computes 4 albedos given :
!         LAI, Greenness fraction, Cosine of the Zenith Angle,
!         Snow depth (set to 0 now since snow-related albedo 
!         calculations have been disabled in this subroutine and 
!         are taken care of in the main code of Mosaic), Veg type,
!         IRUN (set to 1) and Canopy Temp. (set to 0 b/c this 
!         is only needed for snow related calculations that have
!         been disabled in this subroutine
!
!     	Current version of mosaic uses avisdr,anirdr,avisdf,anirdf 
!         that are all set to the albedo taken from EDAS data
!       Using this subroutine, these four parameters can be
!         individually computed and passed back to the main model.
!
!	This subroutine is a modified version of the one written
!	  by Randy Koster / GSFC---it was modified by 
!         Brian Cosgrove / GSFC


!	REAL AVISDR,ANIRDR,AVISDF,ANIRDF,VLAI
!        REAL VGRN,ZTH,SNW,tc
!
!        INTEGER ITYP,irun
!
!
!
!  	   TC=0
!        ityp=4
!        irun=1
!        snw=0.0
!        ZTH=.707
!        VLAI=1.0
!        VGRN=0.6
!
!        call SIBALB(AVISDR, ANIRDR, AVISDF, ANIRDF,
!     &              VLAI, VGRN, ZTH, SNW, ITYP, IRUN,tc)
!
!        print *,'AVISDR, ANIRDR, AVISDF, ANIRDF'
!        print *,'VLAI, VGRN, ZTH, SNW, ITYP, IRUN'
!
!        print *,AVISDR, ANIRDR, AVISDF, ANIRDF
!        print *,VLAI, VGRN, ZTH, SNW, ITYP, IRUN
!        end
! 
!
!
	SUBROUTINE oldSIBALB ( & 
     			 AVISDR, ANIRDR, AVISDF, ANIRDF,& 
     			 VLAI, VGRN, ZTH, SNW, ITYP, IRUN, TC) 
            

	IMPLICIT NONE

	REAL ALVDRS, ALIDRS
	REAL ALVDRD, ALIDRD
	REAL ALVDRI, ALIDRI


	INTEGER IRUN
      INTEGER ITYP (IRUN)

	PARAMETER (ALVDRS = 0.100, ALIDRS = 0.200, & 
             ALVDRD = 0.300,ALIDRD = 0.350, & 
             ALVDRI = 0.700,ALIDRI = 0.700)

! ALVDRS:  Albedo of soil for visible   direct  solar radiation.
! ALIDRS:  Albedo of soil for infra-red direct  solar radiation.
! ALVDFS:  Albedo of soil for visible   diffuse solar radiation.
! ALIDFS:  Albedo of soil for infra-red diffuse solar radiation.


	INTEGER NLAI

	PARAMETER (NLAI = 14 )


	REAL EPSLN, BLAI, DLAI

      PARAMETER (EPSLN = 1.E-6)
      PARAMETER (BLAI = 0.5)
      PARAMETER (DLAI = 0.5)

        REAL ZERO, ONE

      PARAMETER (ZERO=0., ONE=1.0)

	REAL ALATRM

      PARAMETER (ALATRM = BLAI + (NLAI - 1) * DLAI - EPSLN)

        INTEGER NTYPS,UMDNTYPS

      PARAMETER (NTYPS=9)
      PARAMETER (UMDNTYPS=13)


	REAL AVISDR (IRUN), ANIRDR (IRUN), AVISDF (IRUN), ANIRDF (IRUN), & 
     	     VLAI   (IRUN),   VGRN (IRUN),   ZTH  (IRUN),    SNW (IRUN), & 
              TC   (IRUN)

! OUTPUTS:

! AVISDR:   visible, direct albedo.
! ANIRDR:   near infra-red, direct albedo.
! AVISDF:   visible, diffuse albedo.
! ANIRDF:   near infra-red, diffuse albedo.

! INPUTS:

! VLAI:     the leaf area index.
! VGRN:     the greenness index.
! ZTH:      The cosine of the solar zenith angle.
! SNW:      Snow cover in meters water equivalent.


! MOSAIC ITYP: Vegetation type as follows:
!                  1:  BROADLEAF EVERGREEN TREES
!                  2:  BROADLEAF DECIDUOUS TREES
!                  3:  NEEDLELEAF TREES
!                  4:  GROUND COVER
!                  5:  BROADLEAF SHRUBS
!                  6:  DWARF TREES (TUNDRA)
!                  7:  BARE SOIL
!  IRUN: Chip index
!

! UMD ITYP: Vegetation type as follows:
!       1.  Evergreen Needleleaf Forest
!       2.  Evergreen Broadleaf Forest
!       3.  Deciduous Needleleaf Forest
!       4.  Deciduous Broadleaf Forest
!       5.  Mixed Cover
!       6.  Woodland
!       7.  Wooded Grassland
!       8.  Closed Shrubland
!       9.  Open Shrubland
!       10.  Grassland
!       11.  Cropland
!       12.  Bare Ground
!       13.  Urban and Built-Up


! [ Definition of Variables: ]
!
	INTEGER I, J, LAI


	REAL FAC,               GAMMA,          BETA,          ALPHA, & 
     	     DX,                DY,             ALA,           GRN (2), & 
     	     SNWALBold (4, NTYPS), SNWMIDold (NTYPS),  WRMFAC, TICE, & 
            SNWALB (4,UMDNTYPS),SNWMID(UMDNTYPS)



! [ Definition of Functions: ]
!
	REAL oldCOEFF

! --------------------------------------------------



!   Constants used in albedo calculations:

      REAL ALVDRold (NLAI, 2, NTYPS)
      REAL BTVDRold (NLAI, 2, NTYPS)
      REAL GMVDRold (NLAI, 2, NTYPS)
      REAL ALIDRold (NLAI, 2, NTYPS)
      REAL BTIDRold (NLAI, 2, NTYPS)
      REAL GMIDRold (NLAI, 2, NTYPS)

      REAL ALVDR (NLAI, 2, UMDNTYPS)
      REAL BTVDR (NLAI, 2, UMDNTYPS)
      REAL GMVDR (NLAI, 2, UMDNTYPS)
      REAL ALIDR (NLAI, 2, UMDNTYPS)
      REAL BTIDR (NLAI, 2, UMDNTYPS)
      REAL GMIDR (NLAI, 2, UMDNTYPS)

      
!  (Data statements for ALVDR described in full; data statements for
!   other constants follow same framework.)

!    BROADLEAF EVERGREEN (ITYP=1); GREEN=0.33; LAI: .5-7
	DATA (ALVDRold (I, 1, 1), I = 1, 14) & 
     	  /0.0808, 0.0796, 0.0792, 0.0790, 10*0.0789/

!    BROADLEAF EVERGREEN (ITYP=1); GREEN=0.67; LAI: .5-7
	DATA (ALVDRold (I, 2, 1), I = 1, 14) & 
     	  /0.0788, 0.0775, 0.0771, 0.0769, 10*0.0768/

!    BROADLEAF DECIDUOUS (ITYP=2); GREEN=0.33; LAI: .5-7
	DATA (ALVDRold (I, 1, 2), I = 1, 14) & 
     	  /0.0803, 0.0790, 0.0785, 0.0784, 3*0.0783, 7*0.0782/

!    BROADLEAF DECIDUOUS (ITYP=2); GREEN=0.67; LAI: .5-7
	DATA (ALVDRold (I, 2, 2), I = 1, 14) & 
     	  /0.0782, 0.0770, 0.0765, 0.0763, 10*0.0762/

!    NEEDLELEAF (ITYP=3); GREEN=0.33; LAI=.5-7
	DATA (ALVDRold (I, 1, 3), I = 1, 14) & 
     	  /0.0758, 0.0746, 0.0742, 0.0740, 10*0.0739/

!    NEEDLELEAF (ITYP=3); GREEN=0.67; LAI=.5-7
	DATA (ALVDRold (I, 2, 3), I = 1, 14) &
     	  /0.0683, 0.0672, 0.0667, 2*0.0665, 9*0.0664/

!    GROUNDCOVER (ITYP=4); GREEN=0.33; LAI=.5-7
	DATA (ALVDRold (I, 1, 4), I = 1, 14) &
     	  /0.2436, 0.2470, 0.2486, 0.2494, 0.2498, 0.2500, 2*0.2501,&
     		6*0.2502 /
!    GROUNDCOVER (ITYP=4); GREEN=0.67; LAI=.5-7
	DATA (ALVDRold (I, 2, 4), I = 1, 14) /14*0.1637/

!    BROADLEAF SHRUBS (ITYP=5); GREEN=0.33,LAI=.5-7
        DATA (ALVDRold (I, 1, 5), I = 1, 14) &
         /0.0807, 0.0798, 0.0794, 0.0792, 0.0792, 9*0.0791/

!    BROADLEAF SHRUBS (ITYP=5); GREEN=0.67,LAI=.5-7
        DATA (ALVDRold (I, 2, 5), I = 1, 14) &
        /0.0787, 0.0777, 0.0772, 0.0771, 10*0.0770/

!    DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.33,LAI=.5-7
        DATA (ALVDRold (I, 1, 6), I = 1, 14) &
         /0.0802, 0.0791, 0.0787, 0.0786, 10*0.0785/

!    DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.67,LAI=.5-7
        DATA (ALVDRold (I, 2, 6), I = 1, 14) &
         /0.0781, 0.0771, 0.0767, 0.0765, 0.0765, 9*0.0764/


!    BARE SOIL
	DATA (ALVDRold (I, 1, 7), I = 1, 14) /14*ALVDRS/
	DATA (ALVDRold (I, 2, 7), I = 1, 14) /14*ALVDRS/

!    DESERT
	DATA (ALVDRold (I, 1, 8), I = 1, 14) /14*ALVDRD/
	DATA (ALVDRold (I, 2, 8), I = 1, 14) /14*ALVDRD/

!    ICE
	DATA (ALVDRold (I, 1, 9), I = 1, 14) /14*ALVDRI/
	DATA (ALVDRold (I, 2, 9), I = 1, 14) /14*ALVDRI/


        do j=1,2
                do i=1,14
                        ALVDR(i,j,1)=ALVDROLD(i,j,3)

                        ALVDR(i,j,2)=ALVDROLD(i,j,1)

                        ALVDR(i,j,3)=ALVDROLD(i,j,3)

                        ALVDR(i,j,4)=ALVDROLD(i,j,2)

                        ALVDR(i,j,5)=(.5*ALVDROLD(i,j,2) + &
                                     .5*ALVDROLD(i,j,3))

                        ALVDR(i,j,6)=(.5702*ALVDROLD(i,j,3) + & 
                                     .2200*ALVDROLD(i,j,4) + & 
                                     .2098*ALVDROLD(i,j,2))

                        ALVDR(i,j,7)=(.6372*ALVDROLD(i,j,4) + & 
                                     .1955*ALVDROLD(i,j,3) + & 
                                     .1564*ALVDROLD(i,j,2) + & 
                                     .0109*ALVDROLD(i,j,5))

                        ALVDR(i,j,8)=(.4403*ALVDROLD(i,j,8) + & 
                                     .4365*ALVDROLD(i,j,4) + & 
                                     .0743*ALVDROLD(i,j,6) + & 
                                     .0489*ALVDROLD(i,j,5))

                        ALVDR(i,j,9)=(.8506*ALVDROLD(i,j,8) + & 
                                     .0950*ALVDROLD(i,j,5) + & 
                                     .0399*ALVDROLD(i,j,4) + & 
                                     .0145*ALVDROLD(i,j,6))

                        ALVDR(i,j,10)=ALVDROLD(i,j,4)

                        ALVDR(i,j,11)=ALVDROLD(i,j,4)

                        ALVDR(i,j,12)=ALVDROLD(i,j,7)

                        ALVDR(i,j,13)=(.7114*ALVDROLD(i,j,4) + & 
                                      .1055*ALVDROLD(i,j,2) + & 
                                      .0723*ALVDROLD(i,j,3) + & 
                                      .0526*ALVDROLD(i,j,8) + & 
                                      .0178*ALVDROLD(i,j,6) + & 
                                      .0077*ALVDROLD(i,j,5) + & 
                                      .0327*ALVDROLD(i,j,7))
                enddo
        enddo

!****
!**** -------------------------------------------------
	DATA (BTVDRold (I, 1, 1), I = 1, 14) & 
     	  /0.0153, 0.0372, 0.0506, 0.0587, 0.0630, 0.0652, 0.0663, &
     		0.0668, 0.0671, 0.0672, 4*0.0673 & 
     	  /
	DATA (BTVDRold (I, 2, 1), I = 1, 14) & 
     	  /0.0135, 0.0354, 0.0487, 0.0568, 0.0611, 0.0633, 0.0644, & 
     		0.0650, 0.0652, 0.0654, 0.0654, 3*0.0655 & 
     	  /
	DATA (BTVDRold (I, 1, 2), I = 1, 14) &
     	  /0.0148, 0.0357, 0.0462, 0.0524, 0.0554, 0.0569, 0.0576, &
     		0.0579, 0.0580, 0.0581, 0.0581, 3*0.0582 &
     	  /
	DATA (BTVDRold (I, 2, 2), I = 1, 14) &
     	  /0.0131, 0.0342, 0.0446, 0.0508, 0.0539, 0.0554, 0.0560, &
     		0.0564, 0.0565, 5*0.0566 &
     	  /
	DATA (BTVDRold (I, 1, 3), I = 1, 14) &
     	  /0.0108, 0.0334, 0.0478, 0.0571, 0.0624, 0.0652, 0.0666, &
     		0.0673, 0.0677, 0.0679, 4*0.0680 &
     	  /
	DATA (BTVDRold (I, 2, 3), I = 1, 14) &
     	  /0.0034, 0.0272, 0.0408, 0.0501, 0.0554, 0.0582, 0.0597, &
     		0.0604, 0.0608, 0.0610, 4*0.0611 &
     	  /
	DATA (BTVDRold (I, 1, 4), I = 1, 14) &
     	  /0.2050, 0.2524, 0.2799, 0.2947, 0.3022, 0.3059, 0.3076, &
     		0.3085, 0.3088, 0.3090, 4*0.3091 &
     	  /
	DATA (BTVDRold (I, 2, 4), I = 1, 14) &
     	  /0.1084, 0.1404, 0.1617, 0.1754, 0.1837, 0.1887, 0.1915, &
     		0.1931, 0.1940, 0.1946, 0.1948, 0.1950, 2*0.1951 &
     	  /
        DATA (BTVDRold (I, 1, 5), I = 1, 14) &
         /0.0203, 0.0406, 0.0548, 0.0632, 0.0679, 0.0703, 0.0716, &
          0.0722, 0.0726, 0.0727, 0.0728, 0.0728, 0.0728, 0.0729 &
     	  /

        DATA (BTVDRold (I, 2, 5), I = 1, 14) &
         /0.0184, 0.0385, 0.0526, 0.0611,  0.0658, 0.0683, 0.0696, &
          0.0702, 0.0705, 0.0707, 4*0.0708 &
     	  /

        DATA (BTVDRold (I, 1, 6), I = 1, 14) &
         /0.0199, 0.0388, 0.0494,  0.0554, 0.0584, 0.0599, 0.0606, &
          0.0609, 0.0611, 5*0.0612 &
     	  /

        DATA (BTVDRold (I, 2, 6), I = 1, 14) &
         /0.0181, 0.0371, 0.0476, 0.0537,  0.0568, 0.0583, 0.0590, &
          0.0593, 0.0595, 0.0595, 4*0.0596 &
     	  /

	DATA (BTVDRold (I, 1, 7), I = 1, 14) /14*0./
	DATA (BTVDRold (I, 2, 7), I = 1, 14) /14*0./

	DATA (BTVDRold (I, 1, 8), I = 1, 14) /14*0./
	DATA (BTVDRold (I, 2, 8), I = 1, 14) /14*0./

	DATA (BTVDRold (I, 1, 9), I = 1, 14) /14*0./
	DATA (BTVDRold (I, 2, 9), I = 1, 14) /14*0./

        do j=1,2
                do i=1,14
                        BTVDR(i,j,1)=BTVDROLD(i,j,3)

                        BTVDR(i,j,2)=BTVDROLD(i,j,1)

                        BTVDR(i,j,3)=BTVDROLD(i,j,3)

                        BTVDR(i,j,4)=BTVDROLD(i,j,2)

                        BTVDR(i,j,5)=(.5*BTVDROLD(i,j,2) + &
                                     .5*BTVDROLD(i,j,3))

                        BTVDR(i,j,6)=(.5702*BTVDROLD(i,j,3) + &
                                     .2200*BTVDROLD(i,j,4) + &
                                     .2098*BTVDROLD(i,j,2))

                        BTVDR(i,j,7)=(.6372*BTVDROLD(i,j,4) + &
                                     .1955*BTVDROLD(i,j,3) + &
                                     .1564*BTVDROLD(i,j,2) + &
                                     .0109*BTVDROLD(i,j,5))

                        BTVDR(i,j,8)=(.4403*BTVDROLD(i,j,8) + &
                                     .4365*BTVDROLD(i,j,4) + &
                                     .0743*BTVDROLD(i,j,6) + &
                                     .0489*BTVDROLD(i,j,5))

                        BTVDR(i,j,9)=(.8506*BTVDROLD(i,j,8) + &
                                     .0950*BTVDROLD(i,j,5) + & 
                                     .0399*BTVDROLD(i,j,4) + &
                                     .0145*BTVDROLD(i,j,6))

                        BTVDR(i,j,10)=BTVDROLD(i,j,4)

                        BTVDR(i,j,11)=BTVDROLD(i,j,4)

                        BTVDR(i,j,12)=BTVDROLD(i,j,7)

                        BTVDR(i,j,13)=(.7114*BTVDROLD(i,j,4) + &
                                      .1055*BTVDROLD(i,j,2) + &
                                      .0723*BTVDROLD(i,j,3) + &
                                      .0526*BTVDROLD(i,j,8) + &
                                      .0178*BTVDROLD(i,j,6) + &
                                      .0077*BTVDROLD(i,j,5) + &
                                      .0327*BTVDROLD(i,j,7))
                enddo
        enddo


!**** -----------------------------------------------------------
	DATA (GMVDRold (I, 1, 1), I = 1, 14) &
     	  /0.0814, 0.1361, 0.2078, 0.2650, 0.2986, 0.3169,  0.3265, &
      	   0.3313, 0.3337, 0.3348, 0.3354, 0.3357, 2*0.3358 &
     	  /
	DATA (GMVDRold (I, 2, 1), I = 1, 14) &
     	  /0.0760, 0.1336, 0.2034, 0.2622, 0.2969, 0.3159,  0.3259, &
     	   0.3309, 0.3333, 0.3346, 0.3352, 0.3354, 2*0.3356 &
     	  /
	DATA (GMVDRold (I, 1, 2), I = 1, 14) &
     	  /0.0834, 0.1252, 0.1558, 0.1927, 0.2131,   0.2237, 0.2290, &
     	   0.2315, 0.2327, 0.2332, 0.2335, 2*0.2336, 0.2337 &
     	  /
	DATA (GMVDRold (I, 2, 2), I = 1, 14) &
     	  /0.0789, 0.1235, 0.1531, 0.1912, 0.2122, 0.2232,  0.2286, &
     	   0.2312, 0.2324, 0.2330, 0.2333, 0.2334, 2*0.2335 &
     	  /
	DATA (GMVDRold (I, 1, 3), I = 1, 14) &
     	  /0.0647, 0.1342, 0.2215, 0.2968, 0.3432, 0.3696, 0.3838, &
     	   0.3912, 0.3950, 0.3968, 0.3978, 0.3982, 0.3984, 0.3985 &
     	  /
	DATA (GMVDRold (I, 2, 3), I = 1, 14) &
     	  /0.0258, 0.1227, 0.1999, 0.2825, 0.3339, 0.3634, 0.3794, &
     	   0.3877, 0.3919, 0.3940, 0.3950, 0.3956, 0.3958, 0.3959 &
     	  /
	DATA (GMVDRold (I, 1, 4), I = 1, 14) &
     	  /0.3371, 0.5762, 0.7159, 0.7927, 0.8324, 0.8526,  0.8624, &
     	   0.8671, 0.8693, 0.8704, 0.8709, 0.8710, 2*0.8712 &
     	  /
	DATA (GMVDRold (I, 2, 4), I = 1, 14) &
     	  /0.2634, 0.4375, 0.5532, 0.6291, 0.6763, 0.7048, 0.7213, &
     	   0.7310, 0.7363, 0.7395, 0.7411, 0.7420, 0.7426, 0.7428 &
     	  /
        DATA (GMVDRold (I, 1, 5), I = 1, 14) &
         /0.0971, 0.1544, 0.2511, 0.3157, 0.3548, 0.3768, 0.3886, &
          0.3948, 0.3978, 0.3994, 0.4001, 0.4006, 0.4007, 0.4008 &
     	  /

        DATA (GMVDRold (I, 2, 5), I = 1, 14) &
         /0.0924, 0.1470, 0.2458, 0.3123, 0.3527, 0.3756, 0.3877, &
          0.3942, 0.3974, 0.3990, 0.3998, 0.4002, 0.4004, 0.4005 &
     	  /

        DATA (GMVDRold (I, 1, 6), I = 1, 14) &
         /0.0970, 0.1355, 0.1841, 0.2230, 0.2447,  0.2561, 0.2617, &
          0.2645, 0.2658, 0.2664, 0.2667, 3*0.2669 &
     	  /

        DATA (GMVDRold (I, 2, 6), I = 1, 14) &
         /0.0934, 0.1337, 0.1812, 0.2213, 0.2437, 0.2554, 0.2613, &
          0.2642, 0.2656, 0.2662, 0.2665, 0.2667, 0.2667, 0.2668 &
     	  /

	DATA (GMVDRold (I, 1, 7), I = 1, 14) /14*1./
	DATA (GMVDRold (I, 2, 7), I = 1, 14) /14*1./

	DATA (GMVDRold (I, 1, 8), I = 1, 14) /14*1./
	DATA (GMVDRold (I, 2, 8), I = 1, 14) /14*1./

	DATA (GMVDRold (I, 1, 9), I = 1, 14) /14*1./
	DATA (GMVDRold (I, 2, 9), I = 1, 14) /14*1./

        do j=1,2
                do i=1,14
                        GMVDR(i,j,1)=GMVDROLD(i,j,3)

                        GMVDR(i,j,2)=GMVDROLD(i,j,1)

                        GMVDR(i,j,3)=GMVDROLD(i,j,3)

                        GMVDR(i,j,4)=GMVDROLD(i,j,2)

                        GMVDR(i,j,5)=(.5*GMVDROLD(i,j,2) + &
                                     .5*GMVDROLD(i,j,3))

                        GMVDR(i,j,6)=(.5702*GMVDROLD(i,j,3) + &
                                     .2200*GMVDROLD(i,j,4) + &
                                     .2098*GMVDROLD(i,j,2))

                        GMVDR(i,j,7)=(.6372*GMVDROLD(i,j,4) + &
                                     .1955*GMVDROLD(i,j,3) + &
                                     .1564*GMVDROLD(i,j,2) + &
                                     .0109*GMVDROLD(i,j,5))

                        GMVDR(i,j,8)=(.4403*GMVDROLD(i,j,8) + &
                                     .4365*GMVDROLD(i,j,4) + & 
                                     .0743*GMVDROLD(i,j,6) + &
                                     .0489*GMVDROLD(i,j,5))

                        GMVDR(i,j,9)=(.8506*GMVDROLD(i,j,8) + &
                                     .0950*GMVDROLD(i,j,5) + &
                                     .0399*GMVDROLD(i,j,4) + &
                                     .0145*GMVDROLD(i,j,6))

                        GMVDR(i,j,10)=GMVDROLD(i,j,4)

                        GMVDR(i,j,11)=GMVDROLD(i,j,4)

                        GMVDR(i,j,12)=GMVDROLD(i,j,7)

                        GMVDR(i,j,13)=(.7114*GMVDROLD(i,j,4) + &
                                      .1055*GMVDROLD(i,j,2) + &
                                      .0723*GMVDROLD(i,j,3) + &
                                      .0526*GMVDROLD(i,j,8) + &
                                      .0178*GMVDROLD(i,j,6) + &
                                      .0077*GMVDROLD(i,j,5) + &
                                      .0327*GMVDROLD(i,j,7))
                enddo
        enddo



!****  -----------------------------------------------------------

	DATA (ALIDRold (I, 1, 1), I = 1, 14) &
     	  /0.2867,  0.2840, 0.2828, 0.2822, 0.2819, 0.2818, 2*0.2817, &
     	   6*0.2816 &
     	  /
	DATA (ALIDRold (I, 2, 1), I = 1, 14) &
     	  /0.3564, 0.3573, 0.3577, 0.3580, 2*0.3581, 8*0.3582/ 
	DATA (ALIDRold (I, 1, 2), I = 1, 14) &
     	  /0.2848, 0.2819, 0.2804, 0.2798, 0.2795, 2*0.2793, 7*0.2792/
	DATA (ALIDRold (I, 2, 2), I = 1, 14) &
     	  /0.3544, 0.3550, 0.3553, 2*0.3555, 9*0.3556/
	DATA (ALIDRold (I, 1, 3), I = 1, 14) &
     	  /0.2350, 0.2311, 0.2293, 0.2285, 0.2281, 0.2280, 8*0.2279/
	DATA (ALIDRold (I, 2, 3), I = 1, 14) &
     	  /0.2474, 0.2436, 0.2418, 0.2410, 0.2406, 0.2405, 3*0.2404, &
     	   5*0.2403 &
         /
	DATA (ALIDRold (I, 1, 4), I = 1, 14) &
     	  /0.5816, 0.6157, 0.6391, 0.6556, 0.6673, 0.6758, 0.6820, &
     	   0.6866, 0.6899, 0.6924, 0.6943, 0.6956, 0.6966, 0.6974 &
     	  /
	DATA (ALIDRold (I, 2, 4), I = 1, 14) &
     	  /0.5489, 0.5770, 0.5955, 0.6079, 0.6163, 0.6221, 0.6261, &
     	   0.6288, 0.6308, 0.6321, 0.6330, 0.6337, 0.6341, 0.6344 &
     	  /
        DATA (ALIDRold (I, 1, 5), I = 1, 14) &
         /0.2845, 0.2837, 0.2832, 0.2831, 0.2830, 9*0.2829/
        DATA (ALIDRold (I, 2, 5), I = 1, 14) &
         /0.3532, 0.3562, 0.3578,  0.3586, 0.3590, 0.3592, 0.3594, &
          0.3594, 0.3594, 5*0.3595 &
     	  /
        DATA (ALIDRold (I, 1, 6), I = 1, 14) &
         /0.2825, 0.2812, 0.2806, 0.2803, 0.2802, 9*0.2801/
        DATA (ALIDRold (I, 2, 6), I = 1, 14) &
         /0.3512, 0.3538,  0.3552, 0.3559, 0.3562, 0.3564, 0.3565, &
          0.3565, 6*0.3566 &
     	  /

	DATA (ALIDRold (I, 1, 7), I = 1, 14) /14*ALIDRS/
	DATA (ALIDRold (I, 2, 7), I = 1, 14) /14*ALIDRS/

	DATA (ALIDRold (I, 1, 8), I = 1, 14) /14*ALIDRD/
	DATA (ALIDRold (I, 2, 8), I = 1, 14) /14*ALIDRD/

	DATA (ALIDRold (I, 1, 9), I = 1, 14) /14*ALIDRI/
	DATA (ALIDRold (I, 2, 9), I = 1, 14) /14*ALIDRI/

        do j=1,2
                do i=1,14
                        ALIDR(i,j,1)=ALIDROLD(i,j,3)

                        ALIDR(i,j,2)=ALIDROLD(i,j,1)

                        ALIDR(i,j,3)=ALIDROLD(i,j,3)

                        ALIDR(i,j,4)=ALIDROLD(i,j,2)

                        ALIDR(i,j,5)=(.5*ALIDROLD(i,j,2) + &
                                     .5*ALIDROLD(i,j,3))

                        ALIDR(i,j,6)=(.5702*ALIDROLD(i,j,3) + &
                                     .2200*ALIDROLD(i,j,4) + &
                                     .2098*ALIDROLD(i,j,2))

                        ALIDR(i,j,7)=(.6372*ALIDROLD(i,j,4) + &
                                     .1955*ALIDROLD(i,j,3) + &
                                     .1564*ALIDROLD(i,j,2) + &
                                     .0109*ALIDROLD(i,j,5))

                        ALIDR(i,j,8)=(.4403*ALIDROLD(i,j,8) + &
                                     .4365*ALIDROLD(i,j,4) + &
                                     .0743*ALIDROLD(i,j,6) + &
                                     .0489*ALIDROLD(i,j,5))

                        ALIDR(i,j,9)=(.8506*ALIDROLD(i,j,8) + &
                                     .0950*ALIDROLD(i,j,5) + &
                                     .0399*ALIDROLD(i,j,4) + &
                                     .0145*ALIDROLD(i,j,6))

                        ALIDR(i,j,10)=ALIDROLD(i,j,4)

                        ALIDR(i,j,11)=ALIDROLD(i,j,4)

                        ALIDR(i,j,12)=ALIDROLD(i,j,7)

                        ALIDR(i,j,13)=(.7114*ALIDROLD(i,j,4) + &
                                      .1055*ALIDROLD(i,j,2) + &
                                      .0723*ALIDROLD(i,j,3) + &
                                      .0526*ALIDROLD(i,j,8) + &
                                      .0178*ALIDROLD(i,j,6) + &
                                      .0077*ALIDROLD(i,j,5) + &
                                      .0327*ALIDROLD(i,j,7))
                enddo
        enddo



!**** -----------------------------------------------------------
	DATA (BTIDRold (I, 1, 1), I = 1, 14) &
     	  /0.1291, 0.1707, 0.1969, 0.2125, 0.2216,   0.2267, 0.2295, &
     	   0.2311, 0.2319, 0.2323, 0.2326, 2*0.2327, 0.2328 &
     	  /
	DATA (BTIDRold (I, 2, 1), I = 1, 14) &
     	  /0.1939, 0.2357, 0.2598, 0.2735, 0.2810,  0.2851, 0.2874, &
     	   0.2885, 0.2892, 0.2895, 0.2897, 3*0.2898 &
     	  /
	DATA (BTIDRold (I, 1, 2), I = 1, 14) &
     	  /0.1217, 0.1522, 0.1713, 0.1820,   0.1879,  0.1910, 0.1926, &
     	   0.1935, 0.1939, 0.1942, 2*0.1943, 2*0.1944 &
     	  /
	DATA (BTIDRold (I, 2, 2), I = 1, 14) &
     	  /0.1781, 0.2067, 0.2221, 0.2301,   0.2342,  0.2363, 0.2374, &
     	   0.2379, 0.2382, 0.2383, 2*0.2384, 2*0.2385 &
     	  /
	DATA (BTIDRold (I, 1, 3), I = 1, 14) &
     	  /0.0846, 0.1299, 0.1614, 0.1814, 0.1935,   0.2004, 0.2043, &
     	   0.2064, 0.2076, 0.2082, 0.2085, 2*0.2087, 0.2088 &
     	  /
	DATA (BTIDRold (I, 2, 3), I = 1, 14) &
     	  /0.0950, 0.1410, 0.1722, 0.1921, 0.2042, 0.2111,  0.2151, &
     	   0.2172, 0.2184, 0.2191, 0.2194, 0.2196, 2*0.2197 &
     	  /
	DATA (BTIDRold (I, 1, 4), I = 1, 14) &
     	  /0.5256, 0.7444, 0.9908, 1.2700, 1.5680, 1.8505, 2.0767, &
     	   2.2211, 2.2808, 2.2774, 2.2362, 2.1779, 2.1160, 2.0564 &
     	  /
	DATA (BTIDRold (I, 2, 4), I = 1, 14) &
     	  /0.4843, 0.6714, 0.8577, 1.0335, 1.1812, 1.2858, 1.3458, &
     	   1.3688, 1.3685, 1.3546, 1.3360, 1.3168, 1.2989, 1.2838 &
     	  / 
	DATA (BTIDRold (I, 1, 5), I = 1, 14) &
         /0.1498, 0.1930, 0.2201, 0.2364, 0.2460, 0.2514, 0.2544, &
          0.2560, 0.2569, 0.2574, 0.2577, 0.2578, 0.2579, 0.2579 &
     	  /

        DATA (BTIDRold (I, 2, 5), I = 1, 14) &
         /0.2184, 0.2656, 0.2927, 0.3078, 0.3159,  0.3202, 0.3224, &
          0.3235, 0.3241, 0.3244, 0.3245, 3*0.3246 &
     	  /

        DATA (BTIDRold (I, 1, 6), I = 1, 14) &
         /0.1369, 0.1681, 0.1860, 0.1958, 0.2010,  0.2038, 0.2053, &
          0.2060, 0.2064, 0.2066, 0.2067, 3*0.2068 &
     	  /

        DATA (BTIDRold (I, 2, 6), I = 1, 14) &
         /0.1969, 0.2268, 0.2416,  0.2488, 0.2521, 0.2537, 0.2544, &
          0.2547, 0.2548, 5*0.2549 &
     	  / 


	DATA (BTIDRold (I, 1, 7), I = 1, 14) /14*0./
	DATA (BTIDRold (I, 2, 7), I = 1, 14) /14*0./

	DATA (BTIDRold (I, 1, 8), I = 1, 14) /14*0./
	DATA (BTIDRold (I, 2, 8), I = 1, 14) /14*0./

	DATA (BTIDRold (I, 1, 9), I = 1, 14) /14*0./
	DATA (BTIDRold (I, 2, 9), I = 1, 14) /14*0./
        do j=1,2
                do i=1,14
                        BTIDR(i,j,1)=BTIDROLD(i,j,3)

                        BTIDR(i,j,2)=BTIDROLD(i,j,1)

                        BTIDR(i,j,3)=BTIDROLD(i,j,3)

                        BTIDR(i,j,4)=BTIDROLD(i,j,2)

                        BTIDR(i,j,5)=(.5*BTIDROLD(i,j,2) + &
                                     .5*BTIDROLD(i,j,3))

                        BTIDR(i,j,6)=(.5702*BTIDROLD(i,j,3) + &
                                     .2200*BTIDROLD(i,j,4) + &
                                     .2098*BTIDROLD(i,j,2))

                        BTIDR(i,j,7)=(.6372*BTIDROLD(i,j,4) + &
                                     .1955*BTIDROLD(i,j,3) + &
                                     .1564*BTIDROLD(i,j,2) + &
                                     .0109*BTIDROLD(i,j,5))

                        BTIDR(i,j,8)=(.4403*BTIDROLD(i,j,8) + &
                                     .4365*BTIDROLD(i,j,4) + &
                                     .0743*BTIDROLD(i,j,6) + &
                                     .0489*BTIDROLD(i,j,5))

                        BTIDR(i,j,9)=(.8506*BTIDROLD(i,j,8) + &
                                     .0950*BTIDROLD(i,j,5) + &
                                     .0399*BTIDROLD(i,j,4) + &
                                     .0145*BTIDROLD(i,j,6))

                        BTIDR(i,j,10)=BTIDROLD(i,j,4)

                        BTIDR(i,j,11)=BTIDROLD(i,j,4)

                        BTIDR(i,j,12)=BTIDROLD(i,j,7)

                        BTIDR(i,j,13)=(.7114*BTIDROLD(i,j,4) + &
                                      .1055*BTIDROLD(i,j,2) + &
                                      .0723*BTIDROLD(i,j,3) + &
                                      .0526*BTIDROLD(i,j,8) + &
                                      .0178*BTIDROLD(i,j,6) + &
                                      .0077*BTIDROLD(i,j,5) + &
                                      .0327*BTIDROLD(i,j,7))
                enddo
        enddo


!**** --------------------------------------------------------------
	DATA (GMIDRold (I, 1, 1), I = 1, 14) &
     	  /0.1582, 0.2581, 0.3227, 0.3635, 0.3882, 0.4026, 0.4108, &
     	   0.4154, 0.4179, 0.4193, 0.4200, 0.4204, 0.4206, 0.4207 &
     	  /
	DATA (GMIDRold (I, 2, 1), I = 1, 14) &
     	  /0.1934, 0.3141, 0.3818, 0.4200, 0.4415, 0.4533, 0.4598, &
     	   0.4633, 0.4651, 0.4662, 0.4667, 0.4671, 2*0.4672 &
     	  /
	DATA (GMIDRold (I, 1, 2), I = 1, 14) &
     	  /0.1347, 0.1871, 0.2277, 0.2515, 0.2651, 0.2727, 0.2768, &
     	   0.2790, 0.2801, 0.2808, 0.2811, 0.2812, 0.2813, 0.2814 &
     	  /
	DATA (GMIDRold (I, 2, 2), I = 1, 14) &
     	  /0.1440, 0.2217, 0.2629, 0.2839, 0.2947, 0.3003, 0.3031, &
     	   0.3046, 0.3054, 0.3058, 0.3060, 2*0.3061, 0.3062 &
     	  /
	DATA (GMIDRold (I, 1, 3), I = 1, 14) &
     	  /0.1372, 0.2368, 0.3235, 0.3839, 0.4229, 0.4465, 0.4602, &
     	   0.4679, 0.4722, 0.4745, 0.4758, 0.4764, 0.4768, 0.4770 &
     	  /
	DATA (GMIDRold (I, 2, 3), I = 1, 14) &
     	  /0.1435, 0.2524, 0.3370, 0.3955, 0.4332, 0.4563, 0.4697, &
     	   0.4773, 0.4815, 0.4839, 0.4851, 0.4858, 0.4861, 0.4863 &
     	  /
	DATA (GMIDRold (I, 1, 4), I = 1, 14) &
     	  /0.4298, 0.9651, 1.6189, 2.4084, 3.2992, 4.1928, 4.9611, &
     	   5.5095, 5.8085, 5.9069, 5.8726, 5.7674, 5.6346, 5.4944 &
     	  /
	DATA (GMIDRold (I, 2, 4), I = 1, 14) &
     	  /0.4167, 0.8974, 1.4160, 1.9414, 2.4147, 2.7803, 3.0202, &
     	   3.1468, 3.1954, 3.1932, 3.1676, 3.1328, 3.0958, 3.0625 &
     	  /
        DATA (GMIDRold (I, 1, 5), I = 1, 14) &
         /0.1959, 0.3203, 0.3985, 0.4472, 0.4766, 0.4937, 0.5034, &
          0.5088, 0.5117, 0.5134, 0.5143, 0.5147, 0.5150, 0.5152 &
     	  /
 
        DATA (GMIDRold (I, 2, 5), I = 1, 14) &
         /0.2328, 0.3859, 0.4734, 0.5227, 0.5498, 0.5644, 0.5720, &
          0.5761, 0.5781, 0.5792, 0.5797, 0.5800, 0.5802, 0.5802 &
     	  /

        DATA (GMIDRold (I, 1, 6), I = 1, 14) &
         /0.1447, 0.2244, 0.2698, 0.2953, 0.3094, 0.3170, 0.3211, &
          0.3233, 0.3244, 0.3250, 0.3253, 0.3255, 0.3256, 0.3256 &
     	  /

        DATA (GMIDRold (I, 2, 6), I = 1, 14) &
         /0.1643, 0.2624, 0.3110, 0.3347, 0.3461, 0.3517, 0.3543, &
          0.3556, 0.3562, 0.3564, 0.3565, 0.3566, 0.3566, 0.3566 &
     	  /

	DATA (GMIDRold (I, 1, 7), I = 1, 14) /14*1./
	DATA (GMIDRold (I, 2, 7), I = 1, 14) /14*1./

	DATA (GMIDRold (I, 1, 8), I = 1, 14) /14*1./
	DATA (GMIDRold (I, 2, 8), I = 1, 14) /14*1./

	DATA (GMIDRold (I, 1, 9), I = 1, 14) /14*1./
	DATA (GMIDRold (I, 2, 9), I = 1, 14) /14*1./

        do j=1,2
                do i=1,14
                        GMIDR(i,j,1)=GMIDROLD(i,j,3)

                        GMIDR(i,j,2)=GMIDROLD(i,j,1)

                        GMIDR(i,j,3)=GMIDROLD(i,j,3)

                        GMIDR(i,j,4)=GMIDROLD(i,j,2)

                        GMIDR(i,j,5)=(.5*GMIDROLD(i,j,2) + &
                                     .5*GMIDROLD(i,j,3)) 

                        GMIDR(i,j,6)=(.5702*GMIDROLD(i,j,3) + &
                                     .2200*GMIDROLD(i,j,4) + &
                                     .2098*GMIDROLD(i,j,2))

                        GMIDR(i,j,7)=(.6372*GMIDROLD(i,j,4) + &
                                     .1955*GMIDROLD(i,j,3) + &
                                     .1564*GMIDROLD(i,j,2) + &
                                     .0109*GMIDROLD(i,j,5))

                        GMIDR(i,j,8)=(.4403*GMIDROLD(i,j,8) + &
                                     .4365*GMIDROLD(i,j,4) + &
                                     .0743*GMIDROLD(i,j,6) + &
                                     .0489*GMIDROLD(i,j,5))

                        GMIDR(i,j,9)=(.8506*GMIDROLD(i,j,8) + &
                                     .0950*GMIDROLD(i,j,5) + &
                                     .0399*GMIDROLD(i,j,4) + &
                                     .0145*GMIDROLD(i,j,6))

                        GMIDR(i,j,10)=GMIDROLD(i,j,4)

                        GMIDR(i,j,11)=GMIDROLD(i,j,4)

                        GMIDR(i,j,12)=GMIDROLD(i,j,7)

                        GMIDR(i,j,13)=(.7114*GMIDROLD(i,j,4) + &
                                      .1055*GMIDROLD(i,j,2) + &
                                      .0723*GMIDROLD(i,j,3) + &
                                      .0526*GMIDROLD(i,j,8) + &
                                      .0178*GMIDROLD(i,j,6) + &
                                      .0077*GMIDROLD(i,j,5) + &
                                      .0327*GMIDROLD(i,j,7))
                enddo
        enddo

!**** -----------------------------------------------------------

      DATA GRN /0.33, 0.67/

      DATA SNWALBold /.85, .50, .85, .50, &
                  .85, .50, .85, .50, &
                  .85, .50, .85, .50, &
                  .85, .50, .85, .50, &
                  .85, .50, .85, .50, &
                  .85, .50, .85, .50, &
                  .85, .50, .85, .50, &
                  .85, .50, .85, .50, &
                  .85, .50, .85, .50  &
     		  /

      DATA SNWMIDold /50.,50.,50.,2.,50.,2.,2.,2.,2./
      DATA TICE/273.16/


!FPP$ EXPAND (oldCOEFF)

      DO 100 I=1,IRUN
	  ALA = AMIN1 (AMAX1 (ZERO, VLAI(I)), ALATRM)
	  LAI = 1 + MAX(0, INT((ALA-BLAI)/DLAI) )
	  DX = (ALA - (BLAI+(LAI-1)*DLAI)) * (ONE/DLAI)
	  DY = (VGRN(I)- GRN(1)) * (ONE/(GRN(2) - GRN(1)))

	  ALPHA = oldCOEFF (ALVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
	  BETA  = oldCOEFF (BTVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
	  GAMMA = oldCOEFF (GMVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

	  AVISDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
	  AVISDF(I) = ALPHA-BETA &
               + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))

	  ALPHA = oldCOEFF (ALIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
	  BETA  = oldCOEFF (BTIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
	  GAMMA = oldCOEFF (GMIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

	  ANIRDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
	  ANIRDF(I) = ALPHA-BETA &
               + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))


!----------------------------------------------------
!	This section commented out b/c snow albedo already accounted
!	for in main Mosaic code.
!
!	If this part of the code will be used, then the SNWALBold
!	and SNWMIDold need to be mapped to UMD values, and the 
!	section of code below needs to be un-commented out.

!      EXP039: ALLOW REDUCTION IN ALBEDO WHEN SNOW IS CLOSE TO MELTING.
!      (FROM SIB). WE USE 1K (INSTEAD OF .1) BECAUSE RADIATION
!      IS CALLED EVERY 3 HOURS.

!	  IF (SNW (I) .GT. ZERO) THEN
!	   FAC = SNW(I) / (SNW(I) + SNWMID(ITYP(I)))

!           WRMFAC=1.0
!           IF(TC(I) .GT. TICE-1.0) WRMFAC=0.6

!	   AVISDR(I) = AVISDR(I) +
!     &  (SNWALB(1,ITYP(I))*WRMFAC - AVISDR(I)) * FAC
!	   ANIRDR(I) = ANIRDR(I) +
!     &  (SNWALB(2,ITYP(I))*WRMFAC - ANIRDR(I)) * FAC
!	   AVISDF(I) = AVISDF(I) +
!     &  (SNWALB(3,ITYP(I))*WRMFAC - AVISDF(I)) * FAC
!	   ANIRDF(I) = ANIRDF(I) +
!     &  (SNWALB(4,ITYP(I))*WRMFAC - ANIRDF(I)) * FAC
!	  ENDIF

 100  CONTINUE

      RETURN
      END


      FUNCTION oldCOEFF(TABLE, NTABL, LAI ,DX, DY)

      INTEGER NTABL, LAI

      REAL TABLE (NTABL, 2), DX, DY

      oldCOEFF = (TABLE(LAI,  1) & 
           + (TABLE(LAI  ,2) - TABLE(LAI  ,1)) * DY ) * (1.0-DX) &
           + (TABLE(LAI+1,1) &
           + (TABLE(LAI+1,2) - TABLE(LAI+1,1)) * DY ) * DX

      RETURN
      END

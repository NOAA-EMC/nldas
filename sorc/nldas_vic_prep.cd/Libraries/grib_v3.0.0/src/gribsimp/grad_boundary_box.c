/***********************************************************************
File:  	 grad_boundary_box.c
Date:  	 May 1996
Author:  Alice Nakajima, SAIC
Purpose: fill the Pdef, Xdef, Ydef strings needed for GrADS .ctl file; 

Gribsimp used to use 'geom.h' which just has defn for Reg_Geom only;
now use isdg.h;
06/26/97 atn:  +prototypes;
07/02/98 atn:  *latlon proj's longitude computation; * x/ydef printing;
***********************************************************************/
#include <stdio.h>
#include <math.h>
#include "isdb.h"  	/* for reg_geom */
#include "dprints.h" 	/* for dprints */
#include "gribfuncs.h"  /* function prototypes */

#define  min(a, b)  ( (a)<(b) ? (a) : (b) )
#define  max(a, b)  ( (a)>(b) ? (a) : (b) )
enum corners { TL, TR, BL, BR };	/* bottom left,rt, top lt,rt*/


#if PROTOTYPE_NEEDED
extern void xy_latlon (REG_GEOM *geom, double *x, double *y, double *lat, 
		double *lon, int *Status);
extern void latlon_xy (REG_GEOM *geom, double *lat, double *lon, double *x, 
		double *y, int *Status);
#else
extern void xy_latlon ();
extern void latlon_xy ();
#endif

/***********************************************************************
* A.   FUNCTION:  grad_boundary_box
*      PURPOSE:   
*        to find the bounding box coordinates of a particular grid type.
*        Only Lat/Lon (Spherical) projection and Lambert Projection are
*        currently supported.
*
*      ARGUMENTS:
*        char    **Options;     Outp=  charptr, to hold 'Options' line in .ctl 
*        char    **Pdefn; 	Outp=  charptr, to hold 'pdef' line in .ctl 
*        char    **Xdefn; 	Outp=  charptr, to hold 'xdef' line in .ctl
*        char    **Ydefn; 	Outp=  charptr, to hold 'ydef' line in .ctl
*        struct grid_desc_sec *gds; Input= holds Grid Defn Sect info for Msg
*        char   *errmsg;       Outp=  remains empty unless error occurred
*
*      RETURN CODE:
*        0:  no error & fill 'xdef' & "ydef" lines to go into the GrRads 
*            Data control file (.ctl);
*        1:  error occured;
***********************************************************************/

#if PROTOTYPE_NEEDED
int	grad_boundary_box (char **Options, char **Pdefn, char **Xdefn, 
			   char **Ydefn, grid_desc_sec *gds, char *errmsg)
#else
int	grad_boundary_box (Options, Pdefn, Xdefn, Ydefn, gds, errmsg)

char    **Options; 	/* Outp=  charptr, to hold 'Options' line in .ctl */
char    **Pdefn; 	/* Outp=  charptr, to hold 'pdef' line in .ctl */
char    **Xdefn; 	/* Outp=  charptr, to hold 'xdef' line in .ctl */
char    **Ydefn; 	/* Outp=  charptr, to hold 'ydef' line in .ctl */
struct grid_desc_sec *gds; /* Input= holds Grid Defn Sect info for Msg */
char   *errmsg;		/* Outp=  remains empty unless error occurred */
#endif
{
/*
* A.1       DEFAULT to error Status 
*/
   REG_GEOM 	geom;
   char         *func= "grad_boundary_box";
   char   	 *xdef_str, *ydef_str, *pdef_str, *option_str;
   int	 	i,st1;
   double     	dLat, dLon; 
   double	numpts_x, numpts_y;
   double 	x1=-999., y1=-999., x[4], y[4], lat[4], lon[4];
   double   	lat1, lon1, maxlat, maxlon, minlat, minlon, delx, dely;

   DPRINT1 ("Entering %s\n", func);
/*
*
* A.2       ALLOCATE storage for working arrays
*           IF (fails) THEN
*               RETURN error  !msg in Error array
*           ELSE
*               INIT them all to zero
*           ENDIF
*/
  if (! (xdef_str= (char*)malloc(100))|| ! (ydef_str=(char*)malloc(100))||
      ! (option_str= (char*)malloc(100))|| ! (pdef_str= (char*)malloc(150)))
	{ sprintf(errmsg,
	  "%s: Malloc ydef/xdef/option/pdef string failed\n", func);
	  goto BAIL;
	}
  else  { option_str[0]='\0'; pdef_str[0]='\0';
	  xdef_str[0]='\0'; ydef_str[0]='\0';
	}

/*
*
* A.3       SWITCH (type of Projection)
*/
  switch (gds->head.usData_type) 
  {
/*
* A.3.1      Latlon Grid:
*                !calculate Lat, Lon
*                !assign option string to 'option yrev'
*                !initializes xdef, ydef strings
*/
    case LATLON_PRJ:    /* Latitude/Longitude Grid (Equidistant Cylindrical
                  or Plate Carree )         */
	DPRINT2 ("Lat/Lon Projection=%d, scancode=%d\n", 
	gds->head.usData_type, gds->llg.usScan_mode);
	dLat= min (gds->llg.lLat1, gds->llg.lLat2) / 1000.0;

	/* Select Western bounding coordinate ( based on scan mode)  */
	switch (gds->llg.usScan_mode) {
	  case   0: 
	    /* GRIB's 1st pt is top Left  (+i, -j, Horiz first) */
	  case  0x20:  /* dec. 32 */
	    /* GRIB's 1st pt is top Left (+i, -j, Vertical first) */
	  case  0x40:  /* dec. 64 */
	    /* GRIB's 1st pt is Bottom Left (+i, +j, Horiz first) */
	  case  0x60: /* dec. 96 */
	    /* GRIB's 1st pt is Bottom Left (+i, +j, Vertical first) */
		dLon = gds->llg.lLon1 / 1000.0; 
		break;
		
	  default:
	  /*case 0x80:  dec. 128 */
	    /* GRIB's 1st pt is Top Right (-i, -j, Horiz first) */
	  /*case 0xA0: dec. 192 */
	    /* GRIB's 1st pt is Top Right (-i, -j, Vertical first) */
	  /*case 0xC0: dec. 160 */
	    /* GRIB's 1st pt is Bottom Right (-i, +j, Horiz first) */
	  /*case 0xE0: dec. 224 */
	    /* GRIB's 1st pt is Bottom Right (-i, +j, Vertical first) */

		dLon = gds->llg.lLon2 / 1000.0; 
		break;
	}

	if ((gds->llg.usScan_mode&0x40)== 0x00) 
	strcpy(option_str,"options yrev");
	sprintf(xdef_str, "xdef %5d LINEAR %6.3f %.3f",
       	gds->llg.usNi, dLon, (float)gds->llg.iDi/1000.0);
	sprintf(ydef_str, "ydef %5d LINEAR %6.3f %.3f",
       	gds->llg.usNj, dLat, (float)gds->llg.iDj/1000.0);
	break;
	
/*
* A.3.1      Lambert Conformal:
*                !Need to calculate what the Origin is, since GRIB's origin
*                !can be anywhere but Database requires it to be at top Left
*                !and data to flow left to right, top to bottom;
*
*                !Fill GEOM_IN struct with  the internal GDS struct's content
*                !Hard-code:  stor_dsc to "+x in -y"
*/
    case LAMB_PRJ:    /* Lambert conformal */
	DPRINT2("Lambert Conformal Projection=%d, scancode=%d\n", 
	gds->head.usData_type, gds->lam.usScan_mode);
 	strcpy (geom.prjn_name, "lambert");
	geom.nx = gds->lam.iNx;
	geom.ny = gds->lam.iNy;
	geom.lat = (double)gds->lam.lLat1/1000.;
	geom.lon = (double)gds->lam.lLon1/1000.;
	if (geom.lon < -180.) geom.lon+= 360.;
  	else if (geom.lon > 180) geom.lon -= 360.;
	geom.x_int_dis = (double)gds->lam.ulDx /1000.;
	geom.y_int_dis = (double)gds->lam.ulDy /1000.;
	geom.parm_1 = (double)gds->lam.lLat_cut1/1000.;
	geom.parm_2 = (double)gds->lam.lLat_cut2/1000.;
	geom.parm_3 = (double)gds->lam.lLon_orient/1000.;
	if (geom.parm_3 < -180.) geom.parm_3+= 360.;
  	else if (geom.parm_3 > 180) geom.parm_3 -= 360.;
        strcpy (geom.stor_dsc, "+x in -y");


/* 
*                CONVERT GRIB's origin (any corner) to DATABASE's origin
*                (always left top corner)=
*                !Depending on the Scan Mode, find Orig_ix & Orig_iy =
*                !  if GRIB scans in +I dir, set Db's orig_ix= col 1;
*                !  if GRIB scans in -I dir, set Db's orig_ix= #cols;
*                !  if GRIB scans in -J dir, set Db's orig_iy= row 1;
*                !  if GRIB scans in +J dir, set Db's orig_iy= #rows, opts=YREV;
*                !FUNCTION xy_latlon  !convert X/Y coords to Lat/Lon
*                !initializes xdef, ydef, options strings
*/
        DPRINT1("LAM.usScan_mode= %d\n", gds->lam.usScan_mode);
	switch (gds->lam.usScan_mode) {
	  case   0: 
	    /* GRIB's 1st pt is top Left  (+i, -j, Horiz first) */
	    geom.orig_ix= (double) 1; geom.orig_iy= (float) 1; break;

	  case  0x20:  /* dec. 32 */
	    /* GRIB's 1st pt is top Left (+i, -j, Vertical first) */
	    geom.orig_ix= (double) 1; geom.orig_iy= (float) 1; break;

	  case  0x40:  /* dec. 64 */
	    /* GRIB's 1st pt is Bottom Left (+i, +j, Horiz first) */
	    geom.orig_ix= (double) 1; 
	    geom.orig_iy= (float) gds->lam.iNy; 
	    strcpy(option_str,"options yrev"); break;

	  case  0x60: /* dec. 96 */
	    /* GRIB's 1st pt is Bottom Left (+i, +j, Vertical first) */
	    geom.orig_ix= (double) 1; 
	    geom.orig_iy= (float) gds->lam.iNy; 
	    strcpy(option_str,"options yrev"); break;

	  case 0x80: /* dec. 128 */
	    /* GRIB's 1st pt is Top Right (-i, -j, Horiz first) */
	    geom.orig_ix= (double) gds->lam.iNx; 
	    geom.orig_iy= (float) 1; break;

	  case 0xA0: /* dec. 192 */
	    /* GRIB's 1st pt is Top Right (-i, -j, Vertical first) */
	    geom.orig_ix= (double) gds->lam.iNx; 
	    geom.orig_iy= (float) 1; break;

	  case 0xC0: /* dec. 160 */
	    /* GRIB's 1st pt is Bottom Right (-i, +j, Horiz first) */
	    geom.orig_ix= (double) gds->lam.iNx; 
	    geom.orig_iy= (float) gds->lam.iNy;
	    strcpy(option_str,"options yrev"); break;

	  case 0xE0: /* dec. 224 */
	    /* GRIB's 1st pt is Bottom Right (-i, +j, Vertical first) */
	    geom.orig_ix= (double) gds->lam.iNx; 
	    geom.orig_iy= (float) gds->lam.iNy;
	    strcpy(option_str,"options yrev"); break;

	  default : 
	    sprintf(errmsg,"%s: Invalid Scanmode=%d\n", 
	    func, gds->lam.usScan_mode); goto BAIL;
	}

	DPRINT8 ( "\nGEOM: '%s', nx=%d, ny=%d, lat=%.2f, lon=%.2f, "
	"x_int_dis=%.2f,\ny_int_dis=%.2f, parm1=%.2f",
	geom.prjn_name, geom.nx, geom.ny, geom.lat, geom.lon,
	geom.x_int_dis, geom.y_int_dis, geom.parm_1);

	DPRINT4 ( 
	", parm2=%.2f, parm3=%.2f\norig_ix= (double)%f orig_iy= (float)%f\n",
	geom.parm_2, geom.parm_3, geom.orig_ix, geom.orig_iy);

/*
*                !FOR each of the 4 corners:
*                !  FUNCTION xy_latlon     !convert its X/Y to Lat/Lon
*/
       /*   NEONS origin is (1,1),  TOP LEFT CORNER!!
 	    Latitude goes up/down,  Longitude goes left/right

                                PARM3, x=constant 
				  |        
				  |         '*' is highest pt if Parm3 line
			     _  - *-  -   _     is not at the Edge of grid;
			   -      |          -   _
			 -        |	           x (cc,1)TR
			-         |               /
	      TL (1,1) x..........|............../..
			\         |(x,y)	/
			 \        | 	       /
			  \   _ - |-  -  _    /
			   -      |        - _
		      BL(1,rc)	            (cc,rc)BR
	*/
	x[TL]= 1.0; x[TR]= (double) gds->lam.iNx;
	x[BL]= 1.0; x[BR]= (double) gds->lam.iNx;
	y[TL]= 1.0; y[TR]= 1.0;
	y[BL]= (double) gds->lam.iNy; y[BR]= (double) gds->lam.iNy;
	for (i=0; i < 4; i++) 		/* fill Lat,Lon of the 4 corners */
	  {
	    xy_latlon (&geom, &x[i], &y[i], &lat[i], &lon[i], &st1);
	    if (st1 != 0) {  		/* ck Status */
		DPRINT1 ("XY_LATLON has error %d\n",st1);
		sprintf("%s:  xy_latlon has error %d\n", func,st1);
		goto BAIL;
		}
	    DPRINT4("Corner (%.f, %.f) is at Lon/Lat= (%.3f, %.3f)\n",
	    x[i], y[i], lon[i], lat[i]);
	  }

	DPRINT8 ("  X [BL, BR, TL, TR]= %.2lf %.2lf %.2lf %.2lf\n" \
	 	 "  Y [BL, BR, TL, TR]= %.2lf %.2lf %.2lf %.2lf\n" ,
		x[0], x[1], x[2], x[3], y[0], y[1], y[2], y[3]);

	DPRINT8 ("LAT [BL, BR, TL, TR]= %.2lf %.2lf %.2lf %.2lf\n" \
		 "LON [BL, BR, TL, TR]= %.2lf %.2lf %.2lf %.2lf\n",
		 lat[0], lat[1], lat[2], lat[3],
		 lon[0], lon[1], lon[2], lon[3]);

/*
*                FIND the MINIMUM  Latitude (lowest point of grid)
*                !which is the smaller lat of Bottom Left & Bottom Right 
*                !corners
*/
	minlat= min (lat[BL], lat[BR]);
	DPRINT1 ("**MINLAT= min(lat[bottom left & bott rt]= %.3f\n",minlat);

/*
*                    
*                FIND the MAXIMUM Latitude  (hightest point of grid)
*                !IF (Parm3 is on the edge of the grid)
*                !THEN  
*                !   set maxlat to lat at top of grid at parm3
*/
	if (geom.lon == geom.parm_3)  { 
		/* Parm 3 vertical line is on the Edge of the grid,
		   means that the highest row is the Top Left corner
		   of the grid;
		*/
		maxlat=geom.lat;
		DPRINT1 ("/* Parm3  is on edge of grid* /, maxlat=%.3f\n",
		maxlat);
		}

/*
*                !ELSE
*                !   Parm3 is in not on the edge of the grid,
*                !   FUNCTION latlon_xy    !find XY where Lon=Parm3
*                !                         !and Lat=top left corner's Lat
*                !   IF (error) Return error
*                !   Now use x = X1 from above, with y = 1
*                !   FUNCTION xy_latlon    !compute max (Lat,Lon)
*                !   IF (error) Return error
*                !   IF (Parm3 is on grid)
*                !        use max lat from above
*                !   ELSE
*                !        use max lat of TL vs. TR
*                !   ENDIF
*                !   Redefine origin such that ref lon = parm3, ref lat = maxlat
*
*                !ENDIF
*/
	else {
		/* find XY of point where Lon=Parm3, lat=lat of
		   top left corner of grid;
		  */
		latlon_xy (&geom, &lat[TL], &geom.parm_3, &x1, &y1, &st1);
		if (st1!=0) {
		   sprintf(errmsg, "%s:  latlon_xy call has error %d\n",
		   func, st1); goto BAIL; }

		DPRINT4 ("/* PARm3 is middle of grid * /: \n"
			"conv lat/lon (%.3f, Parm3=%.3f) to X1=%.2f y=%.2f\n",
			  lat[TL], geom.parm_3, x1, y1);

		/* Hold X1 constant, compute (lat,lon) at top edge */
		y1= (double) 1;
		xy_latlon (&geom, &x1, &y1, &lat1,&lon1,&st1);
		if (st1!=0) {
		   sprintf(errmsg, "%s:  xy_latlon call has error %d\n",
		   func, st1); goto BAIL; }

		/* if Parm3 line is inside the Grid,
		   then use calculated Max lat from above;

		   if Parm3 line is Outside the Grid, then Max lat is the
		   maximum lat of top left and top right corners;
		*/
                if (lon[TL]<geom.parm_3 && geom.parm_3<lon[TR]) {
			maxlat = lat1;
			DPRINT4("Use (X1 %.2f) and (geom.ny %.2f) conv to:\n"
			"(**MAXLAT=%.3f**, templon=%.3f)\n", 
			x1, y1, maxlat, lon1);
		}
		else
	    	{  maxlat = max(lat[TL], lat[TR]);
	           DPRINT1("/* PArm3 outside of grid */, maxlat=%.3f\n",maxlat);
		}

                /* Re-define origin of grid */
                geom.lat = lat1;
                geom.lon = lon1;
                geom.orig_ix = x1;
                geom.orig_iy = y[BL];
            }
		
/*
*                GET the MIN Longitude (Lon of Top Left corner)
*                GET the MAX Longitude (Lon of Top Right corner)
*/
	minlon= lon[TL]; maxlon= lon[TR];
	DPRINT2("**MINLON=%.3f, MAXLON=%.3f\n", minlon,maxlon);

/*
*                !CONVERT min lat and lon to DB's unit (km to degrees)
*                !PREPARE Pdef, xdef, ydef strings
*/
	/* make min lat/lon a 0.5 deg increment then
	   compute the bounding box parameters;
	   div. by 111 to convert from (km) to (degrees)
	*/
        minlat = floor(minlat *2.)/2.; /* 29.784 to 29.50  */
        minlon = floor(minlon *2.)/2.; /* -123.26 to -123.50 */
	maxlat = ceil (maxlat *2.)/2.; /* 39.977 to 40.0 */
	maxlon = ceil (maxlon *2.)/2.; /* -111.185 to -111.00 */
	delx   = floor((double)geom.x_int_dis/ 111.*1000.)/1000.;
	dely   = floor((double)geom.y_int_dis/ 111.*1000.)/1000.;
	numpts_x = (ceil(maxlon- minlon) / delx) +1;
	numpts_y = (ceil(maxlat- minlat) / dely) +1;

	DPRINT7 ( "ROUND OFF, new minlat= %.2f, minlon=%.2f, "\
	"maxlat=%.2f, maxlon=%.2f\n"\
	"delx= floor(geom.xintdis:%f)/111.*1000/1000= %.3f;  dely=%.3f\n",
	minlat,minlon,maxlat, maxlon, geom.x_int_dis, delx, dely);

	DPRINT8 ("numptsx=(maxlon:%.2f - minlon:%.2f)/(delx:%.2f)+1= %.2lf\n"\
	"numptsY=(maxlat:%.2f - minlat:%.2f)/(dely:%.2f)= %.2f\n",
	maxlon, minlon, delx, numpts_x, maxlat, minlat, dely,numpts_y);

	sprintf(pdef_str, 
	"pdef %d %d lcc %.2f %.2f %.0f %.0f %.0f %.0f %.0f %ld %ld",
		geom.nx, geom.ny, geom.lat, geom.lon, 
		geom.orig_ix, geom.orig_iy,
		geom.parm_1, geom.parm_2, geom.parm_3, 
		gds->lam.ulDx, gds->lam.ulDy);	
  	sprintf(xdef_str, "xdef %.0f linear %6.3f %.3f", 
		numpts_x,minlon,delx); 
  	sprintf(ydef_str, "ydef %.0f linear %6.3f %.3f", 
		numpts_y,minlat,dely); 
	/*
	"Mel2> pdef 79 79 lcc 31 -120  22 10 60 30 -120 15000 15000\n"
	"Mel2> xdef 117 linear -124 .125\n"
	"Mel2> ydef 100 linear 29 .125 \n\n");
	*/
	break;

    default:  sprintf(errmsg,"%s: invalid  gds->head.usData_type=  %d\n",
	func, gds->head.usData_type); 
	goto BAIL;

   /* not supported */
    /* case 4:    /- Gaussian Latitude/Longitude grid */
    /* case 90:    /- Space View perspective or Orthographic Grid    */
    /* case 10:   /- Rotated Lat/Lon */
    /* case 14:   /- Rotated Gaussian */
    /* case 20:   /- Stretched Lat/Lon */
    /* case 24:   /- Stretched Gaussian */
    /* case 30:   /- Stretched and Rotated Lat/Lon */
    /* case 34:   /- Stretched and Rotated Gaussian */
    /* case 8:    /- Albers equal-area, secant or tangent, conical or bipolar */
    /* case 13:   /- Oblique Lambert conformal */
    /* case 1:    /- Mercator Projection Grid               */
    /* case 5:    /- Polar Stereographic Projection Grid    */
/*
* A.3       ENDSWITCH     !type of Projection
*/
  }
  
/*
*
* A.4       ASSIGN the Xdef, Ydef, Pdef, Options string to caller's ptrs;
*
* A.5       RETURN with no errors
*/ 
  *Xdefn= xdef_str;
  *Ydefn= ydef_str;
  *Pdefn= pdef_str;
  *Options= option_str;
  DPRINT1 ("Leaving %s, no errors\n", func);
  return (0);

BAIL:
/*
*
* A.6       IF (error) THEN
*               SET caller's pointers to Null 
*               FREE all the working arrays
*           ENDIF
*
* A.7       RETURN with error Stat
*/ 
  *Xdefn= 0;
  *Ydefn= 0;
  *Pdefn= 0;
  *Options= 0;
  if (xdef_str) free (xdef_str);
  if (ydef_str) free (ydef_str);
  if (pdef_str) free (pdef_str);
  if (option_str) free (option_str);
  DPRINT1 ("Leaving %s, w/ Error Status\n", func);
  return (1);
}

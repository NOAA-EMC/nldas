/*
06/26/97 atn:  +prototypes;
*/
#include <stdio.h>
#include "isdb.h"  /* for struct DATE */
#include "gribfuncs.h"	/* prototypes */

#define ISLEAP(year)    (!((year) % 4) && ((year) % 100) || !((year) % 400))
static int days_sum[12] = {0,31,59,90,120,151,181,212,243,273,304,334};

/***************************************************************/
/*
 * compute day-of-year from year, month and day
 * works for leap years
 */
#if PROTOTYPE_NEEDED
void md_doy(DATE *date, int *doy)
#else
void md_doy(date, doy)
DATE    *date;                          /* DATE structure from isdb.h */
int     *doy;                           /* day of year */
#endif
{
/*      
 *      Input:
 *              date.year  (int)  (year since 0 BC)
 *              date.month (int)  (month in year)
 *              date.day   (int)  (day in month)
 *      Output: doy        (int)  (day-of-year)
 */
        int leap_add = 0;               /* add 1 day if leap year */
        int month;                      /* month */

        month   = date->month;
        if (month < 1 || month > 12) {
                fprintf(stderr, "Md_doy error; month: %d\n", month);
                month = 1;      
        }

        if (ISLEAP(date->year) && month > 2) leap_add = 1;
        *doy     = days_sum[month-1] + date->day + leap_add;
        return;
}
/***************************************************************/
/*
 * convert from human time to epochal time (hours since 00 jan 1, 1970)
 */
#if PROTOTYPE_NEEDED
void h_e_time(DATE *date, double *hour, double *etime)
#else
void h_e_time(date, hour, etime)
DATE    *date;                          /* DATE structure from isdb.h */
double  *hour;                          /* hour in day */
double  *etime;                         /* epochal time */
#endif
{
        int     ytemp;                  /* temporary year holder */
        int     day_cnt;                /* count of days */
        int     doy;                    /* day of year */
        extern void    md_doy();

        md_doy(date, &doy);
        day_cnt = 0;
        if (date->year > 1970) {
                for (ytemp = date->year - 1; ytemp >= 1970; ytemp--) {
                        day_cnt += ISLEAP(ytemp) ? 366 : 365;
                }
        } else if (date->year < 1970) {
                for (ytemp = date->year; ytemp < 1970; ytemp++) {
                        day_cnt -= ISLEAP(ytemp) ? 366 : 365;
                }
        }       
        *etime  = (double) (day_cnt + doy - 1) * 24. + *hour;
        return;
}       

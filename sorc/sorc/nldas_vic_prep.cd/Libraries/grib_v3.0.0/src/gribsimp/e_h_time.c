#include <math.h> 
#include <isdb.h> 

#define ISLEAP(year)    (!((year) % 4) && ((year) % 100) || !((year) % 400))

static int mon_day_cnt[12] = {31,28,31,30,31,30,31,31,30,31,30,31};

/***************************************************************/
/*
 * compute month and day from year and day-of-year
 * works for leap years
 */
void doy_md(doy, date)
int     *doy;                           /* day of year */
DATE    *date;                          /* date structure from isdb.h */
{
/*
 *      Input:
 *              doy         (int)  (day-of-year)
 *              date.year   (int)  (year since 0 BC)
 *      Output: 
 *              date.month  (int)  (month in year) 
 *              date.day    (int)  (day in month)
 */
        int i;                          /* month counter */
        int idoy;                       /* day of year counter */

        if ((idoy = *doy) < 1) {
                date->month = 0;
                date->day   = 0;
                return;
        }
        if (ISLEAP(date->year)) {
                mon_day_cnt[1] = 29;
        } else {
                mon_day_cnt[1] = 28;
        }
        date->month     = 0;
        for (i = 0; i < 12; i++) {
                (date->month)++;
                date->day       = idoy;
                if ((idoy -= mon_day_cnt[date->month-1]) <= 0) {
                        return;
                }
        }
        return;
}

/***************************************************************/
/*
 * convert from epochal time (hours since 00 jan 1, 1970) to human time
 */
void e_h_time(etime, date, hour)
double  *etime;                         /* epochal time */
DATE    *date;                          /* date structure from isdb.h */
double  *hour;                          /* hour in day */
{
        int     ytemp;                  /* temporary year holder */
        int     yr_day_cnt;             /* count of days in year */
        int     doy;                    /* day of year */
        extern void    doy_md();

        doy     = (int) floor(*etime / 24.) + 1;
        *hour   = *etime - (double) (doy - 1) * 24.;
        if (doy > 0) {
                for (ytemp = 1970; ; ytemp++) {
                        yr_day_cnt = ISLEAP(ytemp) ? 366 : 365;
                        if (doy <= yr_day_cnt) break;
                        doy -= yr_day_cnt;
                }
        } else {
                for (ytemp = 1969; ; ytemp--) {
                        yr_day_cnt = ISLEAP(ytemp) ? 366 : 365;
                        doy += yr_day_cnt;
                        if (doy > 0) break;
                }
        }
        date->year = ytemp;
        doy_md(&doy, date);
        return;
}

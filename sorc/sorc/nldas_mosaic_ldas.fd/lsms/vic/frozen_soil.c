//-------------------------------------------------------------------------
// NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
// Released October 2005
//
// See SOFTWARE DISTRIBUTION POLICY for software distribution policies
//
// The LIS source code and documentation are in the public domain,
// available without fee for educational, research, non-commercial and
// commercial purposes.  Users may distribute the binary or source
// code to third parties provided this statement appears on all copies and
// that no charge is made for such copies.
//
// NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
// SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
// IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
// LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
//
// See COPYRIGHT.TXT for copyright details.
//
//-------------------------------------------------------------------------
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "vicNl.h"

#define MAXIT 1000

void setup_frozen_soil(soil_con_struct   *soil_con,
		       layer_data_struct *layer_wet,
		       layer_data_struct *layer_dry,
		       layer_data_struct *layer,
		       energy_bal_struct  energy,
		       int                rec,
		       int                veg,
		       int                Nnodes,
		       float             mu,
		       float            *kappa,
		       float            *Cs,
		       float            *moist) {
/**********************************************************************
  setup_frozen_soil	Keith Cherkauer		July 27, 1998

  This subroutine prepares the data arrays needed to solve the soil
  thermal fluxes through all soil thermal nodes.

  soil_con_struct    soil_con   soil parameter structure
  layer_data_struct *layer_wet  soil variables for wet fraction
  layer_data_struct *layer_dry  soil variables for dry fraction
  layer_data_struct *layer      average soil variables
  energy_bal_struct  energy     energy balance structure
  int                rec        record number
  int                veg        vegetation type number
  int                Nnodes     number of soil thermal nodes 
  float             mu         fraction of grid cell that received precip
  float            *kappa      soil layer thermal conductivity (W/m/K)
  float            *Cs         soil layer heat capacity (J/m^3/K)
  float            *moist      soil layer moisture (mm)
  
**********************************************************************/
  int nidx;

  for(nidx=0;nidx<Nnodes;nidx++) {
    moist[nidx] = energy.moist[nidx];
    kappa[nidx] = energy.kappa_node[nidx];
    Cs[nidx]    = energy.Cs_node[nidx];
  }

}

void finish_frozen_soil_calcs(energy_bal_struct *energy,
			      layer_data_struct *layer_wet,
			      layer_data_struct *layer_dry,
			      layer_data_struct *layer,
			      soil_con_struct   *soil_con,
			      int                Nnodes,
			      int                veg,
			      float             mu,
			      float            *T,
			      float            *kappa,
			      float            *Cs,
			      float            *moist,
			      int               Nlayer,
			      int               frozen_soil_flag) {
/******************************************************************
  finish_frozen_soil_calcs      Keith Cherkauer      July 27, 1998

  This subroutine redistributes soil properties based on the 
  thermal solutions found for the current time step.

******************************************************************/

  int     i;

  find_0_degree_fronts(energy, soil_con->dz_node, T, Nnodes);

  /** Store Layer Temperature Values **/
  for(i=0;i<Nnodes;i++) energy->T[i] = T[i];

  if(energy->Nfrost>0) energy->frozen = TRUE;
  else energy->frozen = FALSE;

  /** Redistribute Soil Properties for New Frozen Soil Layer Size **/
  if(soil_con->FS_ACTIVE && frozen_soil_flag)
    estimate_layer_ice_content(layer_wet, soil_con->dz_node, energy->T,
			     soil_con->max_moist_node, 
			     soil_con->expt_node, soil_con->bubble_node, 
			     soil_con->depth, soil_con->max_moist, 
			     soil_con->expt, soil_con->bubble, 
			     soil_con->bulk_density,
			     soil_con->soil_density, soil_con->quartz,
			     soil_con->layer_node_fract, Nnodes, 
			     Nlayer, soil_con->FS_ACTIVE, frozen_soil_flag);
  /**  if(options.DIST_PRCP)
    estimate_layer_ice_content(layer_dry, soil_con->dz_node, energy->T,
			       soil_con->max_moist_node, 
#if QUICK_FS
			       soil_con->ufwc_table_node,
#else
			       soil_con->expt_node, soil_con->bubble_node, 
#endif
			       soil_con->depth, soil_con->max_moist, 
#if QUICK_FS
			       soil_con->ufwc_table_layer,
#else
			       soil_con->expt, soil_con->bubble, 
#endif
			       soil_con->bulk_density, soil_con->soil_density, 
			       soil_con->quartz, soil_con->layer_node_fract, 
			       Nnodes, options.Nlayer, soil_con->FS_ACTIVE);
  */
}

void solve_T_profile(float *T,
		     float *T0,
		     float *dz,
		     float *kappa,
		     float *Cs,
		     float *moist,
		     float  deltat,
		     float *max_moist,
		     float *bubble,
		     float *expt,
		     float *ice,
		     float *alpha,
		     float *beta,
		     float *gamma,
		     int     Nnodes,
		     char   *FIRST_SOLN,
		     char    FIRST_TIME, 
		     int     FS_ACTIVE, 
		     int    frozen_soil_flag) {
/**********************************************************************
  This subroutine was written to iteratively solve the soil temperature
  profile using a numerical difference equation.  The solution equation
  is second order in space, and first order in time.
**********************************************************************/

  static float A[MAX_NODES];
  static float B[MAX_NODES];
  static float C[MAX_NODES];
  static float D[MAX_NODES];
  static float E[MAX_NODES];

  float *aa, *bb, *cc, *dd, *ee;

  char   Error;
  int    j;
  char NOFLUX = FALSE; 

  if(FIRST_SOLN[0]) {
    FIRST_SOLN[0] = FALSE;
    for(j=1;j<Nnodes-1;j++) {
      A[j] = beta[j-1]*deltat*(kappa[j+1]-kappa[j-1]);
      B[j] = 2.*alpha[j-1]*alpha[j-1]*deltat*kappa[j];
      C[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*Cs[j]*T0[j];
      D[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*ice_density*Lf;
      E[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*Cs[j] 
	+ 4.*kappa[j]*alpha[j-1]*alpha[j-1]*deltat;
    }
    if(NOFLUX) {
      j = Nnodes-1;
      A[j] = beta[j-1]*deltat*(kappa[j]-kappa[j-1]);
      B[j] = 2.*alpha[j-1]*alpha[j-1]*deltat*kappa[j];
      C[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*Cs[j]*T0[j];
      D[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*ice_density*Lf;
      E[j] = alpha[j-1]*alpha[j-1]*beta[j-1]*Cs[j] 
	+ 4.*kappa[j]*alpha[j-1]*alpha[j-1]*deltat;
    }
  }
    
  aa = &A[0];
  bb = &B[0];
  cc = &C[0];
  dd = &D[0];
  ee = &E[0];

  for(j=0;j<Nnodes;j++) T[j]=T0[j];

  Error = calc_soil_thermal_fluxes(Nnodes, T, T0, moist, max_moist, ice, 
				   bubble, expt, alpha, gamma, aa, bb, cc, 
				   dd, ee, FS_ACTIVE, frozen_soil_flag);
  
}
  

int calc_soil_thermal_fluxes(int     Nnodes,
			     float *T,
			     float *T0,
			     float *moist,
			     float *max_moist,
			     float *ice,
			     float *bubble,
			     float *expt,
			     float *alpha,
			     float *gamma,
			     float *A, 
			     float *B, 
			     float *C, 
			     float *D, 
			     float *E,
			     char    FS_ACTIVE,
			     int frozen_soil_flag) {
  
  /** Eventually the nodal ice contents will also have to be updated **/


  int    Error;
  char   Done;
  int    j;
  int    ItCount;
  float threshold = 1.e-2;	/* temperature profile iteration threshold */
  float maxdiff;
  float diff;
  float oldT;
  float fprime;
  char NOFLUX = FALSE;
  Error = 0;
  Done = FALSE;
  ItCount = 0;
  
  while(!Done && Error==0 && ItCount<MAXIT) {
    ItCount++;
    maxdiff=threshold;
    for(j=1;j<Nnodes-1;j++) {
      oldT=T[j];
      
      /**	2nd order variable kappa equation **/
      fprime = (T[j+1]-T[j-1])/alpha[j-1];
      
      if(T[j] >= 0 || !FS_ACTIVE || !frozen_soil_flag) {
	T[j] = (A[j]*(T[j+1]-T[j-1])
		+ B[j]*(T[j+1]+T[j-1]-gamma[j-1]*fprime)
		+ C[j] + D[j]*(0.-ice[j])) / (E[j]);
      }
      else {

	T[j] = root_brent(T0[j]-(SOIL_DT),
			  T0[j]+(SOIL_DT),soil_thermal_eqn, 
			  T[j+1], T[j-1], T0[j], moist[j], max_moist[j], 
			  bubble[j], expt[j], ice[j], gamma[j-1], fprime, 
			  A[j], B[j], C[j], D[j], E[j]);
	
	if(T[j] <= MISSING+1) 
	  error_solve_T_profile(T[j], T[j+1], T[j-1], T0[j], moist[j], 
				max_moist[j], bubble[j], expt[j], ice[j], 
				gamma[j-1], fprime, A[j], B[j], C[j], D[j], 
				E[j]);
      }
      
      diff=fabs(oldT-T[j]);
      if(diff > maxdiff) maxdiff=diff;
    }
    
    if(NOFLUX) { 
      /** Solve for bottom temperature if using no flux lower boundary **/
      oldT=T[Nnodes-1];
      
      fprime = (T[Nnodes-1]-T[Nnodes-2])/alpha[Nnodes-2];
      
      j = Nnodes-1;
      
      if(T[j] >= 0 || !FS_ACTIVE || !frozen_soil_flag) {
	T[j] = (A[j]*(T[j]-T[j-1]) + B[j]*(T[j] + T[j-1] - gamma[j-1]*fprime)
		+ C[j] + D[j]*(0.-ice[j])) / E[j];
      }
      else {

	T[Nnodes-1] = root_brent(T0[Nnodes-1]-SOIL_DT,T0[Nnodes-1]
				 +SOIL_DT,soil_thermal_eqn, T[Nnodes-1],
				 T[Nnodes-2], T0[Nnodes-1], 
				 moist[Nnodes-1], max_moist[Nnodes-1], 
				 bubble[j], expt[Nnodes-1], ice[Nnodes-1], 
				 gamma[Nnodes-2], fprime, 
				 A[j], B[j], C[j], D[j], E[j]);
	
	if(T[j] <= MISSING+1) 
	  error_solve_T_profile(T[Nnodes-1], T[Nnodes-1],
				T[Nnodes-2], T0[Nnodes-1], 
				moist[Nnodes-1], max_moist[Nnodes-1], 
				bubble[Nnodes-1], 
				expt[Nnodes-1], ice[Nnodes-1], 
				gamma[Nnodes-2], fprime, 
				A[j], B[j], C[j], D[j], E[j]);
      }
      
      diff=fabs(oldT-T[Nnodes-1]);
      if(diff>maxdiff) maxdiff=diff;
    }
    
    if(maxdiff <= threshold) Done=TRUE;
    
  }
  
  if(!Done && !Error) {
    fprintf(stderr,"ERROR: Temperature Profile Unable to Converge!!!\n");
    fprintf(stderr,"Dumping Profile Temperatures (last, new).\n");
    for(j=0;j<Nnodes;j++) fprintf(stderr,"%f\t%f\n",T0[j],T[j]);
    //    vicerror("ERROR: Cannot solve temperature profile:\n\tToo Many Iterations in solve_T_profile");
  }

  return (Error);

}

float error_solve_T_profile (float Tj, ...) {

  va_list ap;

  float error;

  va_start(ap,Tj);
  error = error_print_solve_T_profile(Tj, ap);
  va_end(ap);

  return error;

}

float error_print_solve_T_profile(float T, va_list ap) {

  float TL;
  float TU;
  float T0;
  float moist;
  float max_moist;
  float bubble;
  float expt;
  float ice0;
  float gamma;
  float fprime;
  float A;
  float B;
  float C;
  float D;
  float E;

  int CELL_NUM = -1;
  char lismsg[MAXSTRING];

  TL        = (float) va_arg(ap, double);
  TU        = (float) va_arg(ap, double);
  T0        = (float) va_arg(ap, double);
  moist     = (float) va_arg(ap, double);
  max_moist = (float) va_arg(ap, double);
  bubble    = (float) va_arg(ap, double);
  expt      = (float) va_arg(ap, double);
  ice0      = (float) va_arg(ap, double);
  gamma     = (float) va_arg(ap, double);
  fprime    = (float) va_arg(ap, double);
  A         = (float) va_arg(ap, double);
  B         = (float) va_arg(ap, double);
  C         = (float) va_arg(ap, double);
  D         = (float) va_arg(ap, double);
  E         = (float) va_arg(ap, double);
  
  sprintf(lismsg,"DBG: error_print -- TL = %f [ %d ]",TL,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- TU = %f [ %d ]",TU,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- T0 = %f [ %d ]",T0,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- moist = %f [ %d ]",moist,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- max_moist = %f [ %d ]",max_moist,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- bubble = %f [ %d ]",bubble,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- expt = %f [ %d ]",expt,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- ice0 = %f [ %d ]",ice0,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- gamma = %f [ %d ]",gamma,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- fprime = %f [ %d ]",fprime,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- A = %f [ %d ]",A,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- B = %f [ %d ]",B,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- C = %f [ %d ]",C,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- D = %f [ %d ]",D,CELL_NUM);
  lis_log_msgC(lismsg);
  sprintf(lismsg,"DBG: error_print -- E = %f [ %d ]",E,CELL_NUM);
  lis_log_msgC(lismsg);

  //  vicerror("Finished dumping values for solve_T_profile.\nTry increasing SOIL_DT to get model to complete cell.\nThen check output for instabilities.");
  
  return(0.0);

}

#undef MAXIT

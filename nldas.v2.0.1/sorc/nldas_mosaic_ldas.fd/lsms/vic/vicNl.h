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
#include <stdio.h>
#include <math.h>
#include "vicNl_def.h"
#include "ftn.h"
#include "lis_headers.h"

float arno_evap(layer_data_struct *, layer_data_struct *, float, float, 
                float, float, float, float, float, float, float, 
                float, float, float, float, float, float, float);

void CalcAerodynamic(char, int, int, float, float, float, float, float *, 
                     float *, float *, float, float *, float *);

float calc_rainonly(float,float,float);

float CalcSnowPackEnergyBalance(float Tsurf, ...);

void calc_root_fractions_(int *,int *);

int calc_soil_thermal_fluxes(int, float *, float *, float *, float *, 
                             float *, float *, float *, float *, 
                             float *, float *, float *, float *, 
                             float *, float *, char, int);

float calc_surf_energy_bal(int, int, int, int, int, int, int, int,
                           float, float, float, float, float, 
                           float, float, float, float, float, 
                           float, float, float, float, float, 
                           float, float, float *, float *, float *,
                           float *, float *, atmos_data_struct *, 
                           veg_var_struct *, veg_var_struct *, 
                           energy_bal_struct *, snow_data_struct *, 
                           layer_data_struct *, layer_data_struct *, 
                           soil_con_struct *, dmy_struct *, 
                           int, int, int, int);

float calc_veg_displacement(float);

float calc_veg_height(float);

float calc_veg_roughness(float);

float canopy_evap(layer_data_struct *, layer_data_struct *, veg_var_struct *, 
                  veg_var_struct *, char, int, int, float, float *, float, 
                  float, float, float, float, float, float, float, float, 
                  float, float *, float *, float *, float *, float *, int);

void compute_dz(float *, float *, int , float); 

//void compute_dz(double *, float *, int , double); 

void compute_penman_constants(float, float, float, float, float, 
                                float, float, float, float, float *, 
                                float *, float *, float *, float *);

void compute_soil_layer_thermal_properties(layer_data_struct *, float *,
                                           float *, float *, float *, 
                                           int);

void distribute_node_moisture_properties(float *, float *, float *, float *, 
                                        float *, float *, float *, float *, 
                                        float *, float *, float *, float *, 
                                        float *, float *, int, int, char, int);

float ErrorPrintSnowPackEnergyBalance(float, va_list);

float ErrorSnowPackEnergyBalance(float Tsurf, ...);

float error_calc_surf_energy_bal(float Tsurf, ...);

float error_print_solve_T_profile(float, va_list);

float error_solve_T_profile(float Tsurf, ...);

void estimate_layer_ice_content(layer_data_struct *, float *, float *,
                              float *, float *, float *, float *,
                              float *, float *, float *, float *,
                              float *, float *, float [MAX_LAYERS+1][MAX_NODES],
                              int, int, char, int);

float estimate_T1(float, float, float, float, float, 
                  float, float, float, float, float, float);

float exp_interp(float,float,float,float,float);

void find_0_degree_fronts(energy_bal_struct *, float *, float *, int);

layer_data_struct find_average_layer(layer_data_struct *, layer_data_struct *,
                                     float, float);

void finish_frozen_soil_calcs(energy_bal_struct *, layer_data_struct *,
                              layer_data_struct *, layer_data_struct *,
                              soil_con_struct *, int, int, float, 
                              float *, float *, float *, float *, int , int);

void full_energy(int, atmos_data_struct *, soil_con_struct *,
                 veg_con_struct *, dist_prcp_struct *,
                 dmy_struct, int, int, int, int, int, int, int);

float func_surf_energy_bal(float, va_list);

void FTN(initialize_model_state)(int *, int *, int *, int *, int *, int *, float *);
void FTN(initialize_model_state2)(int *, int *, int *, int *, int *, int *);

void initialize_snow(snow_data_struct **, int, int, int);

void initialize_soil(cell_data_struct **, soil_con_struct *, int, int, int);

void initialize_veg(veg_var_struct **, int, int);

float linear_interp(float, float, float, float, float);

cell_data_struct **make_cell_data(int, int);

void make_dist_prcp(int *);

energy_bal_struct **make_energy_bal(int, int);

snow_data_struct **make_snow_data(int, int);

veg_var_struct **make_veg_var(int, int);

float massrelease(float *, float *, float *, float *);

float maximum_unfrozen_water(float, float, float, float);

float modify_Ksat(float, int);

float new_snow_density(float);

void nrerror(char *);

FILE  *open_file(char string[], char type[]);

void prepare_full_energy(int, int, int, dist_prcp_struct *, soil_con_struct *, 
                         float *, float *, int, int, int);

float quick_penman(float, float, float, float, float, float, float, float);

void FTN(read_initial_model_state)(char *, int *, int *, int *, int *);

//void  FTN(read_parammap)(char *, int *,
//                         char *, int *,
//                         char *, int *,
//                         char *, int *,
//                         char *, int *,
//                         char *, int *,
//                         char *, int *,
//                         int *, int *, int *, 
//                         int *, int *, int *);

void FTN(read_snowband)(int *, int *);

//void FTN(read_soilparam)(char *, int *, int *, int *, int *, char *, int *, 
//                         char *, int *, int *, int *, int *, int *);

void FTN(read_veglib)(char *, int *, int *, int *);

void FTN(read_vegparam)(int *);

float root_brent(float, float, float (*Function)(float, va_list), ...);

void runoff(layer_data_struct *, layer_data_struct *, energy_bal_struct *, 
            soil_con_struct *, float *, float *, float *, float *, float *, 
            float, int, int, int, int, int, int, int);

void setup_frozen_soil(soil_con_struct *, layer_data_struct *, 
                       layer_data_struct *, layer_data_struct *, 
                       energy_bal_struct , int, int, int, float, 
                       float *, float *, float *);

void set_node_parameters(float *, float *, float *, float *, float *, float *, 
                         float *, float *, float *, float *, float *, float *, 
                         float [MAX_LAYERS+1][MAX_NODES], int, int, char);

float SnowPackEnergyBalance(float, va_list);

float snow_albedo(float, float, float, float, int, char);

float snow_density(int, float, float, float, float, float, float, float);

void snow_intercept(float, float, float, float, float, float, float, 
                    float, float, float, float, float, float, float, 
                    float *, float *, float *, float *, float *, float *, 
                    float *, float *, int, int);

void snow_melt(soil_con_struct *, int, float, float, float, snow_data_struct *,
               float, float, float, float, float, float, float, float, float, 
               float, float, float, float, float, float, float *, float *, 
               float *, float *, float *, float *, float *, float *);

float soil_conductivity(float, float, float, float, float);

float soil_thermal_eqn(float, va_list);

float solve_snow(int, snow_data_struct *, layer_data_struct *, 
                 layer_data_struct *, veg_var_struct *, veg_var_struct *, 
                 int, int, energy_bal_struct *, soil_con_struct *, 
                 char, int, int, int, int, int, int, int, 
                 float, float, float, float, float, float, float, float, 
                 float, float, float, float, float, float, float, float, float,
                 float *, float *, float *, float *, float *, float *, float *,
                 float *, float *, float *, float *, float *, float *, 
                 float *, float *, int, int, int);

float solve_surf_energy_bal(float Tsurf, ...);

void  solve_T_profile(float *, float *, float *, float *,float *,
                      float *, float, float *, float *, float *,
                      float *, float *, float *, float *, int, char *, 
                      char, int, int);

float StabilityCorrection(float, float, float, float, float, float);

float svp(float);

float svp_slope(float);

void surface_fluxes(int, char, int, int, int, int, int, int, int, int, int, 
                    float, float, float, float, float, float,
                    float, float, float *, float *, float *, float *, 
                    float *, float *, float *, float *, float *, float *,
                    float *, float *, float *, float *, float *, 
                    atmos_data_struct *, soil_con_struct *, dmy_struct, 
                    energy_bal_struct *, snow_data_struct *, 
                    layer_data_struct *, layer_data_struct *, 
                    veg_var_struct *, veg_var_struct *, int, int, int, int);

void transpiration(layer_data_struct *, int, int, float, float, float, 
                   float, float, float, float, float, float, float, 
                   float *, float *, float *, float *, float *, 
                   float *, float *, int);

void FTN(vic_allocate)(int *, int *, int *, int *);

void FTN(vic_almaout)(int *, int *,int *, int*, int *, int *, int *, int *);

void FTN(vic_singleout)(int *, int *, int *, int *, int *, int *, int *, int *,
                        int *, int *, int *, int *, int *, char *);

void FTN(vic_bcast)(int *);

void FTN(vic_f2t)(int *, float *, int *);

void FTN(vic_gather)();

void FTN(vic_singlegather)(int *, float *);

float vic_penman(float, float, float, float, float, float, float, 
                 float, float, float, float);

void FTN(vic_run)(int *, int *, int *, int *, int *, int *, 
                  int *, int *, int *, int *, int *, int *, int *);

void FTN(vic_scatter)();

void vic_delcalcs();

void FTN(vic_soil_paramalloc)(int *, int *, int *);

/* void vic_stats(float *, float , int ,float *, float *, float *, float *); */

void FTN(vic_totinit)();

float volumetric_heat_capacity(float, float, float);

void FTN(write_model_state)(char *, int *, int *, int *, int *, int *);

void MassRelease(float *, float *, float *, float *); 

void make_dist_prcp(int *);


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
//BOP
//
//  !ROUTINE: LIS_soils_FTable
//  
//
// !DESCRIPTION:
//   Function table implementation for different soils options
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>

#include "ftn_drv.h"
typedef struct
{ 
  void (*func)(float*);
} SAND_TABLE; 
SAND_TABLE sand_table[6];

typedef struct
{ 
  void (*func)(float*);
} CLAY_TABLE; 
CLAY_TABLE clay_table[6];

typedef struct
{ 
  void (*func)(float*);
} SILT_TABLE; 
SILT_TABLE silt_table[2];

typedef struct
{ 
  void (*func)(float*);
} W_SAT_TABLE; 
W_SAT_TABLE w_sat_table[3];

typedef struct
{ 
  void (*func)(float*);
} W_SAT_MATP_TABLE; 
W_SAT_MATP_TABLE w_sat_matp_table[3];

typedef struct
{ 
  void (*func)(float*);
} W_SAT_HYDC_TABLE; 
W_SAT_HYDC_TABLE w_sat_hydc_table[3];

typedef struct
{ 
  void (*func)(float*);
} W_BPOWER_TABLE; 
W_BPOWER_TABLE w_bpower_table[3];

typedef struct
{ 
  void (*func)(float*);
} W_WILT_TABLE; 
W_WILT_TABLE w_wilt_table[3];

typedef struct
{ 
  void (*func)(float*);
} SOIL_CLASS_TABLE; 
SOIL_CLASS_TABLE soil_class_table[6];

//BOP
// !ROUTINE: registerreadsand
//  
// !DESCRIPTION: Registers the routines to open and 
// read sand data
// 
// !INTERFACE:
void FTN(registerreadsand)(int *i,void (*func)())
  //EOP
{ 
  sand_table[*i].func = func; 
}

//BOP
// !ROUTINE: readsand
//  
// !DESCRIPTION: Delegates the routines for 
// reading sand files
// 
// !INTERFACE:
void FTN(readsand)(int *i,float *array)
//EOP
{ 
  sand_table[*i].func(array); 
}
//BOP
// !ROUTINE: registerreadclay
//  
// !DESCRIPTION: Registers the routines to open and 
// read clay data
// 
// !INTERFACE:
void FTN(registerreadclay)(int *i,void (*func)())
  //EOP
{ 
  clay_table[*i].func = func; 
}
//BOP
// !ROUTINE: readsand
//  
// !DESCRIPTION: Delegates the routines for 
// reading sand files
// 
// !INTERFACE:
void FTN(readclay)(int *i,float *array)
//EOP
{ 
  clay_table[*i].func(array); 
}

//BOP
// !ROUTINE: registerreadsilt
//  
// !DESCRIPTION: Registers the routines to open and 
// read silt data
// 
// !INTERFACE:
void FTN(registerreadsilt)(int *i,void (*func)())
  //EOP
{ 
  silt_table[*i].func = func; 
}
//BOP
// !ROUTINE: readsand
//  
// !DESCRIPTION: Delegates the routines for 
// reading sand files
// 
// !INTERFACE:
void FTN(readsilt)(int *i,float *array)
//EOP
{ 
  silt_table[*i].func(array); 
}

//BOP
// !ROUTINE: registerreadwsat
//  
// !DESCRIPTION: Registers the routines to open and 
// read maximum soil moisture content data.
// 
// !INTERFACE:
void FTN(registerreadwsat)(int *i,void (*func)())
  //EOP
{ 
  w_sat_table[*i].func = func; 
}

//BOP
// !ROUTINE: readwsat
//  
// !DESCRIPTION: Delegates the routines for 
// reading maximum soil moisture content files.
// 
// !INTERFACE:
void FTN(readwsat)(int *i,float *array)
//EOP
{ 
  w_sat_table[*i].func(array); 
}

//BOP
// !ROUTINE: registerreadwsatmatp
//  
// !DESCRIPTION: Registers the routines to open and 
// read saturated soil potential data.
// 
// !INTERFACE:
void FTN(registerreadwsatmatp)(int *i,void (*func)())
  //EOP
{ 
  w_sat_matp_table[*i].func = func; 
}

//BOP
// !ROUTINE: readwsatmatp
//  
// !DESCRIPTION: Delegates the routines for 
// reading saturated soil potential files.
// 
// !INTERFACE:
void FTN(readwsatmatp)(int *i,float *array)
//EOP
{ 
  w_sat_matp_table[*i].func(array); 
}

//BOP
// !ROUTINE: registerreadwsathydc
//  
// !DESCRIPTION: Registers the routines to open and 
// read saturated soil hydraulic conductivity data.
// 
// !INTERFACE:
void FTN(registerreadwsathydc)(int *i,void (*func)())
  //EOP
{ 
  w_sat_hydc_table[*i].func = func; 
}

//BOP
// !ROUTINE: readwsathydc
//  
// !DESCRIPTION: Delegates the routines for 
// reading saturated soil hydraulic conductivity files.
// 
// !INTERFACE:
void FTN(readwsathydc)(int *i,float *array)
//EOP
{ 
  w_sat_hydc_table[*i].func(array); 
}

//BOP
// !ROUTINE: registerreadwbpower
//  
// !DESCRIPTION: Registers the routines to open and 
// read b parameter data.
// 
// !INTERFACE:
void FTN(registerreadwbpower)(int *i,void (*func)())
  //EOP
{ 
  w_bpower_table[*i].func = func; 
}

//BOP
// !ROUTINE: readwbpower
//  
// !DESCRIPTION: Delegates the routines for 
// reading b parameter files.
// 
// !INTERFACE:
void FTN(readwbpower)(int *i,float *array)
//EOP
{ 
  w_bpower_table[*i].func(array); 
}

//BOP
// !ROUTINE: registerreadwwilt
//  
// !DESCRIPTION: Registers the routines to open and 
// read wilting point parameter data.
// 
// !INTERFACE:
void FTN(registerreadwwilt)(int *i,void (*func)())
  //EOP
{ 
  w_wilt_table[*i].func = func; 
}

//BOP
// !ROUTINE: readwwilt
//  
// !DESCRIPTION: Delegates the routines for 
// reading wilting point parameter files.
// 
// !INTERFACE:
void FTN(readwwilt)(int *i,float *array)
//EOP
{ 
  w_wilt_table[*i].func(array); 
}

//BOP
// !ROUTINE: registerreadsoilclass
//  
// !DESCRIPTION: Registers the routines to open and 
// read soil classification data.
// 
// !INTERFACE:
void FTN(registerreadsoilclass)(int *i,void (*func)())
  //EOP
{ 
  soil_class_table[*i].func = func; 
}

//BOP
// !ROUTINE: readsoilclass
//  
// !DESCRIPTION: Delegates the routines for 
// reading soil classification files.
// 
// !INTERFACE:
void FTN(readsoilclass)(int *i,float *array)
//EOP
{ 
  soil_class_table[*i].func(array); 
}


#ifndef INTEGRATE_H_
#define INTEGRATE_H_

/* standard C stuff */
#include <stdlib.h>
/* serial N_Vector types, fcts., and macros */
#include <nvector/nvector_serial.h>
/* gsl random number generation */
#include <gsl/gsl_rng.h>
/* model specific header */
#include "model.h"


/* stochastic differential equation? */
#ifndef __SDE__
#define __SDE__ 0
#endif
/* default MaxNumSteps */
#ifndef __CVODE_MAX_NUM_STEPS__
#define __CVODE_MAX_NUM_STEPS__ 1000
#endif
/* default MaxStep */
#ifndef __CVODE_MAX_STEP__
#define __CVODE_MAX_STEP__ 0
#endif
/* default sparse solver */
#ifndef __CVODE_SPARSE__
#define __CVODE_SPARSE__ 0
#endif
/* default absolute tolerance (scalar) */
#ifndef __CVODE_ATOL__
#define __CVODE_ATOL__ 1.0e-8
#endif
/* default relative tolerance */
#ifndef __CVODE_RTOL__
#define __CVODE_RTOL__ 1.0e-8
#endif




/* function declarations */
int  integrate   ( size_t n_sim_pts, double * t_sim, double * x_sim0, double * x_sim, struct function_data * data, double seed );
int  add_state_noise ( N_Vector x, N_Vector noise );
void print_value ( char * name, realtype value );
int  check_flag  ( void * flagvalue, char * funcname, int opt );

#endif

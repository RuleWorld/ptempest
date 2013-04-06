#ifndef MODEL_H_
#define MODEL_H_

/* standard C stuff */
#include <stdlib.h>
/* serial N_Vector types, fcts., and macros */
#include <nvector/nvector_serial.h>
/* gsl random number generation */
#include <gsl/gsl_rng.h>


/* include header with model specific defintions */
#ifndef CUSTOM_H_
#define CUSTOM_H_
#include "rhs.h"
#endif



/* struct containing disturbances */
struct  disturbance
{
    size_t   n_pts;
    double * time;
    double * x[__N_DIST__];
};

/* struct containing controls */
struct  control
{
    size_t   n_pts;
    double * time;
    double * x[__N_CTRL__];
};

/* struct combining params, disturbances and controls */
struct  function_data
{
    struct model_parameters  * par;
    struct disturbance       * dist;
    struct control           * ctrl;
};


/* function declarations */
int  rhs_func ( realtype t, N_Vector x, N_Vector xdot, void *f_data );
int  generate_state_noise ( realtype t0, realtype t1, N_Vector x, N_Vector noise, gsl_rng *rand_gen, struct function_data * data );

#endif

/*
**   file: rhs.h
**
**   Example header file for RHS function: Lorenz attractor.
**   --J.Hogg, 19Apr2012
**   
**/

/* Problem Dimensions */
#define __N_EQNS__        3    /* number of equations */
#define __N_PARAMS__      3    /* number of total parameters */
#define __N_CTRL__        1    /* number of control variables */
#define __N_DIST__        1    /* number of distrurbance variables */

/* stochastic differential equation? */
#define __SDE__  1

/* CVode integrator configuration */
#define __CVODE_MAX_NUM_STEPS__  1000
#define __CVODE_MAX_STEP__       0
#define __CVODE_SPARSE__         0
#define __CVODE_ATOL__           1.0e-7
#define __CVODE_RTOL__           1.0e-9

/* process noise parameter */
#define  __PROCESS_NOISE_SIGMA__  0.2

/*Indices for state variables */
#define  __X_IDX__  0
#define  __Y_IDX__  1
#define  __Z_IDX__  2

/* Friendly names for state variables */
#define  __X__  NV_Ith_S(x,__X_IDX__)
#define  __Y__  NV_Ith_S(x,__Y_IDX__)
#define  __Z__  NV_Ith_S(x,__Z_IDX__)

/* Friendly names for rates */
#define  __DOT_X__  NV_Ith_S(xdot,__X_IDX__)
#define  __DOT_Y__  NV_Ith_S(xdot,__Y_IDX__)
#define  __DOT_Z__  NV_Ith_S(xdot,__Z_IDX__)


/* struct containing all model parameters */
struct model_parameters
{
    double sigma;
    double rho;
    double beta;
};


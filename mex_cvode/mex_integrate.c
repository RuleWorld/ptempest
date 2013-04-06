/*   
**   FILE: mex_integrate.c 
**
**   MEX-wrapper for CVode integrator
**
**   AUTHOR:  Justin Hogg
**   EMAIL:   justin.s.hogg@gmail.com
**   UPDATE:  19 April 2012
**
**   rhs.c contains model specific functions:
**     rhs_func             : right hand side (derivative) of the ODE
**     generate_state_noise : noise generating function for SDE models
**
**   rhs.h contains model specific defintions:
**     N_EQNS           : number of model equations
**     N_PARAMS         : number of model parameters
**     N_CTRL           : number of control variables
**     N_DIST           : number of disturbance variables
**     model_parameters : parameter names
**
**   Requires the CVode (v2.6) and GSL libraries.
**   To compile in Matlab:
**      mex -lgsl -lgslcblas -lsundials_nvecserial -lsundials_cvode
**          mex_integrate.c integrate.c functions.c rhs.c -o MODEL
**
**   To run from the Matlab console: 
**      [err_flag, x_sim] = MODEL( t_sim, x_sim0, params, ...
**                                 t_dist, x_dist, t_ctrl, x_ctrl, seed );
**
**      t_sim   -- column vector of timepoints to return from the simulation
**      x_sim0  -- row vector of initial values
**      params  -- row vector of parameters
**      t_dist  -- column vector of disturbance times
**      x_dist  -- array (or column vector) of disturbance values
**      t_ctrl  -- column vector of control times
**      x_ctrl  -- array (or column vector) of control values
**      seed    -- a seed for generating random process noise
**                   (no process noise if seed==0)
**
**   (note: rows always correspond to timepoints!)
*/




/**
***   COMPILER DIRECTIVES, ETC
**/

/* Standard C stuff */
#include <stdlib.h>
#include <math.h>
/* matlab MEX headers */
#include "mex.h"
#include "matrix.h" 
/* headers for CVode integrator */
#include "integrate.h"


/* function declarations */
void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );


/*
 *  main MEX: process input, call the integrator, return output
 */
 
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    /* output variables */
    double * err_flag;
    double * x_sim;
    /* input variables */
    double * t_sim;
    double * x_sim0;
    double * params;
    double * t_dist;
    double * x_dist;
    double * t_ctrl;
    double * x_ctrl;
    double   seed;
    
    /* iterator */
    size_t  ii;
    /* number of time points */
    size_t  n_sim_pts;
    size_t  n_dist_pts;
    size_t  n_ctrl_pts;

    /* structures for function data */
    struct disturbance       dist;
    struct control           ctrl;
    struct function_data     data;


    /* Start handling input/output arguments from Matlab */
    /* check number of input/output arguments */
    if (nlhs != 2)
    {   mexErrMsgTxt
        (   "invalid syntax:  [err_flag, x_sim] = MODEL( t_sim, x_sim0, params, t_dist, x_dist, t_ctrl, x_ctrl, seed )"   );
    }
    if (nrhs != 8)
    {   mexErrMsgTxt
        (   "invalid syntax:  [err_flag, x_sim] = MODEL( t_sim, x_sim0, params, t_dist, x_dist, t_ctrl, x_ctrl, seed )"   );
    }    
    

    /* check for correct dimensions */
    /* t_sim */
    if ( mxGetN(prhs[0]) != 1 )
    {   mexErrMsgTxt("invalid argument: TIME must be a column vector");   }

    /* x_sim0 */
    if ( (mxGetM(prhs[1]) != 1)  ||  (mxGetN(prhs[1]) != __N_EQNS__) )
    {   mexErrMsgTxt("invalid argument: X_SIM0 must be a row vector with N_EQNS elements");   } 

    /* params */
    if ( (mxGetM(prhs[2]) != 1)  ||  (mxGetN(prhs[2]) != __N_PARAMS__) )
    {   mexErrMsgTxt("invalid argument: PARAMS must be a row vector with N_PARAMS elements");   }

    /* t_dist */
    if ( mxGetN(prhs[3]) > 1 )
    {   mexErrMsgTxt("invalid argument: T_DIST must be a column vector");   }

    /* u_dist */
    if ( mxGetM(prhs[4])*mxGetN(prhs[4]) != __N_DIST__*mxGetM(prhs[3]) )
    {   mexErrMsgTxt("invalid argument: X_DIST must be a column vector with 0*|T_DIST| elements");   }

    /* t_ctrl */
    if ( mxGetN(prhs[5]) > 1 )
    {   mexErrMsgTxt("invalid argument: T_CTRL must be a column vector");   }   

    /* u_ctrl */
    if ( mxGetM(prhs[6])*mxGetN(prhs[6]) != __N_CTRL__*mxGetM(prhs[5]) )
    {   mexErrMsgTxt("invalid argument: X_CTRL must be a column vector with 1*|T_CTRL| elements");   }

    /* seed */
    if ( (mxGetM(prhs[7]) != 1)  ||  (mxGetN(prhs[7]) != 1) )
    {   mexErrMsgTxt("invalid argument: SEED must be a scalar value");   }

   
    /* get pointers to input arrays */
    t_sim   = mxGetPr(prhs[0]);
    x_sim0  = mxGetPr(prhs[1]);
    params  = mxGetPr(prhs[2]);
    t_dist  = mxGetPr(prhs[3]);
    x_dist  = mxGetPr(prhs[4]);
    t_ctrl  = mxGetPr(prhs[5]);
    x_ctrl  = mxGetPr(prhs[6]);
    /* get random number seed */
    seed    = *mxGetPr(prhs[7]);
   
    /* get number of timepoints */
    n_sim_pts  = mxGetM(prhs[0]);
    n_dist_pts = mxGetM(prhs[3]);
    n_ctrl_pts = mxGetM(prhs[5]);

    /* Create an mxArray for output trajectories */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL );
    plhs[1] = mxCreateDoubleMatrix(n_sim_pts, __N_EQNS__, mxREAL);

    /* get pointers to output arrays */
    err_flag = mxGetPr(plhs[0]);
    x_sim    = mxGetPr(plhs[1]);


    /* set up disturbance structure */
    dist.n_pts  = n_dist_pts;
    dist.time   = t_dist;
    for ( ii = 0; ii < __N_DIST__; ++ii )
    {
        dist.x[ii] = x_dist + ii*__N_DIST__;
    }

    /* setup control structure */
    ctrl.n_pts = n_ctrl_pts;
    ctrl.time  = t_ctrl;
    for ( ii = 0; ii < __N_CTRL__; ++ii )
    {
        ctrl.x[ii] = x_ctrl + ii*__N_DIST__;
    }

    /* setup data  */
    data.par  = (struct model_structure *)params;  /* TODO: suppress warning? */
    data.dist = &dist;
    data.ctrl = &ctrl;

    /* integrate system */
    err_flag[0] = integrate( n_sim_pts, t_sim, x_sim0, x_sim, &data, seed );
    return;

}















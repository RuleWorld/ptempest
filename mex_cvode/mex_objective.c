/*   
**   FILE: mex_objective.c 
**
**   MEX-wrapper for calculating objective function
**
**   AUTHOR:  Justin Hogg
**   EMAIL:   justin.s.hogg@gmail.com
**   UPDATE:  19 Apr 2012
**
**   rhs.c contains model specific functions:
**     rhs_func             : right hand side (derivative) of the ODE
**     generate_state_noise : noise generating function 
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
**          mex_objective.c integrate.c functions.c rhs.c -o OBJECTIVE
**
**   To run from the Matlab console: 
**      [val] = OBJECTIVE( t_sim, x_sim0, params, t_dist, x_dist, ...
**                         t_ctrl, x_ctrl, x_ref, w_sim, w_ctrl, w_dctrl );
**
**      t_sim   -- column vector of timepoints for objective evaluation
**      x_sim0  -- row vector of initial values
**      params  -- row vector of parameters
**      t_dist  -- column vector of disturbance times
**      x_dist  -- array (or column vector) of disturbance values
**      t_ctrl  -- column vector of control times
**      x_ctrl  -- array (or column vector) of control values
**      x_ref   -- reference trajectory
**      w_sim   -- objective weights for deviation from reference (vector)
**      w_ctrl  -- objective weights for controls (vector)
**      w_dctrl -- objective weights for delta controls (vector)
**
**   (note: rows always correspon to timepoints!)
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
 *  main MEX: process input, call the integrator, calculate objective, return output
 */
 
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    /* output variables */
    double * objective;
    /* input variables */
    double * t_sim;
    double * x_sim0;
    double * params;
    double * t_dist;
    double * x_dist;
    double * t_ctrl;
    double * x_ctrl;
    double * x_ref;
    double * w_sim;
    double * w_ctrl;
    double * w_dctrl;

    /* An array for simulation trajectories. NOTE: ideally this would be a 2-D
       array, but integrate() expects a double pointer, so we stick with this. */
    double * x_sim;
    /* holds return flag from integrate() */
    int      err_flag;
    /* random number seed. seed==0 turns off random noise */
    double   seed = 0;
    
    /* number of time points */
    size_t  n_sim_pts;
    size_t  n_dist_pts;
    size_t  n_ctrl_pts;
    /* iterators */
    size_t  ii, jj;
    

    /* objective components */
    double  objective_ctrl;
    double  objective_dctrl;
    double  objective_sim;

    /* structures for function data */
    struct disturbance    dist;
    struct control        ctrl;
    struct function_data  data;
  
    
    /* Start handling input/output arguments from Matlab */
    /* check number of input/output arguments */
    if (nlhs != 1)
    {   mexErrMsgTxt
        (   "syntax: [val] = OBJECTIVE( t_sim, x_sim0, params, t_dist, x_dist, t_ctrl, x_ctrl, x_ref, w_sim, w_ctrl, w_dctrl )"   );
    }
    if (nrhs != 11 )
    {   mexErrMsgTxt
        (   "syntax: [val] = OBJECTIVE( t_sim, x_sim0, params, t_dist, x_dist, t_ctrl, x_ctrl, x_ref, w_sim, w_ctrl, w_dctrl )"   );
    }    


    /* check for correct dimensions */
    /* t_sim */
    if ( mxGetN(prhs[0]) != 1 )
    {   mexErrMsgTxt("invalid argument: T_SIM must be a column vector");   }

    /* x_sim0 */
    if ( (mxGetM(prhs[1]) != 1)  ||  (mxGetN(prhs[1]) != __N_EQNS__) )
    {   mexErrMsgTxt("invalid argument: X_SIM0 must be a row vector with N_EQNS elements");   } 

    /* params */
    if ( (mxGetM(prhs[2]) != 1)  ||  (mxGetN(prhs[2]) != __N_PARAMS__) )
    {   mexErrMsgTxt("invalid argument: PARAMS must be a row vector with N_PARAMS elements");   }

    /* t_dist */
    if ( mxGetN(prhs[3]) > 1 )
    {   mexErrMsgTxt("invalid argument: T_DIST must be a column vector");   }

    /* x_dist */
    if ( mxGetM(prhs[4])*mxGetN(prhs[4]) != mxGetM(prhs[3]) )
    {   mexErrMsgTxt("invalid argument: X_DIST must be a column vector with |T_DIST| elements");   }

    /* t_ctrl */
    if ( mxGetN(prhs[5]) > 1 )
    {   mexErrMsgTxt("invalid argument: T_CTRL must be a column vector");   }   

    /* x_ctrl */
    if ( mxGetM(prhs[6])*mxGetN(prhs[6]) != __N_CTRL__*mxGetM(prhs[5]) )
    {   mexErrMsgTxt("invalid argument: X_CTRL must be a column vector with 4*|T_CTRL| elements");   }

    /* x_ref */
    if ( (mxGetM(prhs[7]) != mxGetM(prhs[0]) )  ||  (mxGetN(prhs[7]) != __N_EQNS__) )
    {   mexErrMsgTxt("invalid argument: X_REF is an array with size |T_SIM| x N_EQNS");   }
    
    /* w_sim */    
    if ( (mxGetM(prhs[8]) != 1)  ||  (mxGetN(prhs[8]) != __N_EQNS__) )
    {   mexErrMsgTxt("invalid argument: W_SIM is a row vector with length N_EQNS");   }

    /* w_ctrl */    
    if ( (mxGetM(prhs[9]) != 1)  ||  (mxGetN(prhs[9]) != __N_CTRL__) )
    {   mexErrMsgTxt("invalid argument: W_CTRL is a row vector with length = number of control channels");   }

    /* w_dctrl */    
    if ( (mxGetM(prhs[10]) != 1)  ||  (mxGetN(prhs[10]) != __N_CTRL__) )
    {   mexErrMsgTxt("invalid argument: W_DCTRL is aa row vector with length = number of control channels");   }

   
    /* get pointers to input arrays */
    t_sim   = mxGetPr(prhs[0]);
    x_sim0  = mxGetPr(prhs[1]);
    params  = mxGetPr(prhs[2]);
    t_dist  = mxGetPr(prhs[3]);
    x_dist  = mxGetPr(prhs[4]);
    t_ctrl  = mxGetPr(prhs[5]);
    x_ctrl  = mxGetPr(prhs[6]);
    x_ref   = mxGetPr(prhs[7]);
    w_sim   = mxGetPr(prhs[8]);
    w_ctrl  = mxGetPr(prhs[9]);
    w_dctrl = mxGetPr(prhs[10]);
   
    /* get number of timepoints */
    n_sim_pts  = mxGetM(prhs[0]);
    n_dist_pts = mxGetM(prhs[3]);
    n_ctrl_pts = mxGetM(prhs[5]);

    /* Create an mxArray for output objective */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL );
    /* get pointer to output arrays */
    objective = mxGetPr(plhs[0]);


    /* Allocate memory for x_sim */
    x_sim = (double *) malloc( n_sim_pts*__N_EQNS__ * sizeof(double) );
    if (x_sim == NULL)
    {   /* Memory could not be allocated, abort! */
        mexErrMsgTxt("Unable to allocate dynamic memory for simulation trajectory.");
    }


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


    /* initialize objective */
    *objective = 0.0;
        
    /* calculate control contribution to the objective! */
    objective_ctrl = 0.0;
    /* loop over ctrl time points, starting from the second point
         NOTE: the first point is ignored, since it's not part of the optimization */
    for ( ii = 1;  ii < n_ctrl_pts;  ++ii )
    {
        /* loop over control channels */
        for ( jj = 0;  jj < __N_CTRL__;  ++jj )
        {
            objective_ctrl += w_ctrl[jj] * pow( x_ctrl[ jj*n_ctrl_pts + ii ], 2.0);
        }
    }
    *objective += objective_ctrl;
    
    
    /* calculate delta control contribution to the objective! */
    objective_dctrl = 0.0;
    /* loop over ctrl time points
        NOTE: first point is only needed to calculate the delta step of the second control point */
    for ( ii = 1;  ii < n_ctrl_pts;  ++ii )
    {
        /* loop over control channels */    
        for ( jj = 0;  jj < __N_CTRL__;  ++jj )
        {
            objective_dctrl += w_dctrl[jj] * pow( (x_ctrl[ jj*n_ctrl_pts + ii ] - x_ctrl[ jj*n_ctrl_pts + ii-1 ] ), 2.0);
        }
    }
    *objective += objective_dctrl;



    /* Integrate System */
    err_flag = integrate( n_sim_pts, t_sim, x_sim0, x_sim, &data, seed );

    if (err_flag == 0)
    {   /* calculate trajectory error contribution to objective (initial timepoint is ignored) */
        objective_sim = 0.0;
        for ( ii = 1;  ii < n_sim_pts;  ++ii )
        {
            for ( jj = 0;  jj < __N_EQNS__;  ++jj )
            {
                objective_sim += w_sim[jj] * pow( x_sim[jj*n_sim_pts + ii] - x_ref[jj*n_sim_pts + ii], 2.0);
            }
        }

        /* finalize objective */
        *objective += objective_sim;
    }
    else
    {   /* error! */
        mexErrMsgTxt("Some problem integrating model equations.");
    }
    
    /* free up memory allocatd for x_sim */
    free(x_sim);
    x_sim = NULL;

    return;
}















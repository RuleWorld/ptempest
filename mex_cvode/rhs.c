/*
**   FILE: rhs.c
**
**   Example source file for RHS function: Lorenz attractor
**   with control channel on X, disturbance on Y, and state noise.
**   --J.Hogg, 19Apr2012
**
**   See mex_integrate.c for compilation instructions.  
**
**   dX/dt = sigma*(Y - X) + U
**   dY/dt = X*(rho - Z) - Y + D
**   dZ/dt = X*Y - beta*Z
**
**   where: U=control channel
**          D=disturbance
**
**   Trial simulation in MATLAB:
**   >> t = [0:.01:40]';
**   >> x0 = [1 1 1];
**   >> pars = [10 28 8/3];
**   >> [err, x] = lorenz(t, x0, pars, [], [], [], [], 0);
**   >> plot3( x(:,1), x(:,2), x(:,3) );
**
**   Add a bit of state noise by passing a random number seed
**   >> seed = 93257;
**   >> [err, x] = lorenz(t, x0, pars, [], [], [], [], seed);
**   >> plot3( x(:,1), x(:,2), x(:,3) );
**
**   Add source term to X through the control channel at t=20
**   >> tc = [0 20]';  % control time points
**   >> uc = [0 100]'; % piecewise constant values
**   >> [err, x] = lorenz(t, x0, pars, [], [], tc, uc, seed);
**   >> plot3( x(:,1), x(:,2), x(:,3) );
**
**/

/* math */
#include "math.h"
/* model header */
#include "model.h"
/* headers for useful mathematical functions */
#include "functions.h"
/* gsl random number generation */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/* User defined RHS function */
int rhs_func (realtype t, N_Vector x, N_Vector xdot, void *f_data)
{
    /* local variables */
    struct function_data     * data;
    struct model_parameters  * par;
    struct disturbance       * dist;
    struct control           * ctrl;
    
    /* cast f_data to 'function_data' type */
    data = (struct function_data *)f_data;

    /* get pointer to parameters, disturbance and control structures */
    par  = data->par;
    dist = data->dist;
    ctrl = data->ctrl;


    /*** USER DEFINED CONTENT STARTS HERE ***/

    /* calculate ctrl value at this time point */
    size_t ii;  
    double ctrl_X = 0.0;
    if ( ctrl->n_pts > 0 )
    {   /*
         * control trajectories are piecewise constant:
         *  use a heaviside approximation to avoid discontinuity
         */
        ctrl_X = ctrl->x[0][0] * heaviside_approx(t - ctrl->time[0]);
        for ( ii = 1;  ii < ctrl->n_pts;  ++ii )
        {   ctrl_X += (ctrl->x[0][ii] - ctrl->x[0][ii-1]) * heaviside_approx(t - ctrl->time[ii]);   }    
    }

    /* calculate disturbance values at this time point */ 
    double dist_Y = 0.0;
    if ( dist->n_pts > 0 )
    {   /*
         * disturbance trajectories are piecewise constant:
         *  use a heaviside approximation to avoid discontinuity
         */
        dist_Y = dist->x[0][0] * heaviside_approx(t - dist->time[0]);
        for ( ii = 1;  ii < dist->n_pts;  ++ii )
        {   dist_Y += (dist->x[0][ii] - dist->x[0][ii-1]) * heaviside_approx(t - dist->time[ii]);   }    
    }


    /* calculate RHS function now */

    /* X rate equation */
    __DOT_X__ = par->sigma*(__Y__ - __X__)  +  ctrl_X;

    /* Y rate equation */
    __DOT_Y__ = __X__*(par->rho - __Z__) - __Y__  +  dist_Y;

    /* Z rate equation */
    __DOT_Z__ = __X__*__Y__ - par->beta*__Z__;

        
    /* debugging output */
    /*
    mexPrintf( "t = %.3e  X = %.3e  X'(t) = %.3e\n", t, __X__, __DOT_X__ );
    mexPrintf( "t = %.3e  Y = %.3e  Y'(t) = %.3e\n", t, __Y__, __DOT_Y__ );
    mexPrintf( "t = %.3e  Z = %.3e  Z'(t) = %.3e\n", t, __Z__, __DOT_Z__ );
    */    

    /*** USER DEFINED CONTENT ENDS HERE ***/
 
   
    /* All done.  RHS evaluation is loaded into vector xdot. */
    return(0);
}




/* Generate state noise (but don't apply) */
int generate_state_noise (realtype t0, realtype t1, N_Vector x, N_Vector noise, gsl_rng *rand_gen, struct function_data * data)
{
    /* local variables */
    realtype                   sig;
    struct model_parameters  * par;
    struct disturbance       * dist;
    struct control           * ctrl;

    /* get pointer to parameters, disturbance and control structures */
    par  = data->par;
    dist = data->dist;
    ctrl = data->ctrl;
    
    /* calculate mu and sig */
    sig = sqrt(t1 - t0) * __PROCESS_NOISE_SIGMA__;


    /*** USER DEFINED CONTENT STARTS HERE ***/

    /* generate and add state noise (Weiner-type noise) */
    NV_Ith_S(noise, __X_IDX__) = gsl_ran_gaussian( rand_gen, sig );
    NV_Ith_S(noise, __Y_IDX__) = gsl_ran_gaussian( rand_gen, sig );
    NV_Ith_S(noise, __Z_IDX__) = gsl_ran_gaussian( rand_gen, sig );

    /*** USER DEFINED CONTENT ENDS HERE ***/


    /* All done. State noise was added to x vector. */
    return(0);
}




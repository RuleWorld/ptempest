/*   
**   FILE: integrate.c 
**
**   CVode integrator with support for step-disturbances and control variables
**   derived examples included in the SUNDIALS distributions
**
**   AUTHOR:  Justin Hogg
**   EMAIL:   justin.s.hogg@gmail.com
**   UPDATE:  13feb2012
*/


/* headers */
#include "integrate.h"
/* CVode headers */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <cvode/cvode_spgmr.h>       /* prototype for CVSpgmr */
/* gsl random number generation */
#include <gsl/gsl_rng.h>


int integrate ( size_t n_sim_pts, double * t_sim, double * x_sim0, double * x_sim, struct function_data * data, double seed )
{
    /* CVODE specific variables */
    realtype  t;
    N_Vector  x;
    N_Vector  noise;
    void      *cvode_mem;
    int       flag;

    /* index counters */
    size_t    ii;
    size_t    jj;

    /* variables for GSL random number generation */
    const gsl_rng_type * T;
    gsl_rng * rand_gen;

    /* initialize GSL random number generator */
    gsl_rng_env_setup();  
    T = gsl_rng_default;
    rand_gen = gsl_rng_alloc (T);
    gsl_rng_set( rand_gen, (unsigned long)seed );
    

    /* Create serial vector to hold state */
    x = NULL;
    x = N_VNew_Serial(__N_EQNS__);
    if (check_flag((void *)x, "N_VNew_Serial", 0))
    {   /* Free dynamic memory */
        gsl_rng_free( rand_gen );
        return 1;
    }

    /* Create serial vector to hold noise */
    noise = NULL;
    noise = N_VNew_Serial(__N_EQNS__);
    if (check_flag((void *)noise, "N_VNew_Serial", 0))
    {   /* Free dynamic memory */
        N_VDestroy_Serial(x);
        gsl_rng_free( rand_gen ); 
        return 1;
    }

    /* Create serial vector to hold vector atol */
    #ifdef __CVODE_VECTOR_ATOL__
    N_Vector atol = NULL;
    atol = N_VNew_Serial(__N_EQNS__);
    if (check_flag((void *)atol, "N_VNew_Serial", 0))
    {   /* Free dynamic memory */
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(noise);   
        gsl_rng_free( rand_gen );    
        return 1;
    }
    /* load vector absolute tolerances */
    double atol_values [__N_EQNS__] = __CVODE_VECTOR_ATOL__;
    for ( jj = 0; jj < __N_EQNS__; jj++ )
    {
        NV_Ith_S(atol,jj) = atol_values[jj];
    }
    #endif

    /* copy initial conditions into state vector and the output array*/
    for ( jj = 0; jj < __N_EQNS__; jj++ )
    {
        x_sim[jj*n_sim_pts] = x_sim0[jj];
        NV_Ith_S(x,jj)      = x_sim0[jj];
        /* set noise to zero */
        NV_Ith_S(noise,jj)  = 0.0;
    }

    
    /*   Call CVodeCreate to create the solver memory:    
     *   CV_BDF     specifies the Backward Differentiation Formula
     *   CV_NEWTON  specifies a Newton iteration
     *   A pointer to the integrator problem memory is returned and stored in cvode_mem.
     */
    cvode_mem = NULL;     
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0))
    {   /* Free dynamic memory */
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(noise);
        #ifdef __CVODE_VECTOR_ATOL__
        N_VDestroy_Serial(atol);
        #endif
        gsl_rng_free( rand_gen );
        return 1;
    }
  
    /*  Call CVodeInit to initialize the integrator memory */
    flag = CVodeInit(cvode_mem, rhs_func, t_sim[0], x);
    if (check_flag(&flag, "CVodeInit", 1))
    {   /* Free dynamic memory */
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(noise);
        #ifdef __CVODE_VECTOR_ATOL__
        N_VDestroy_Serial(atol);
        #endif
        CVodeFree( &cvode_mem );
        gsl_rng_free( rand_gen );
        return 1;
    }
    
    #ifdef __CVODE_VECTOR_ATOL__
    flag = CVodeSVtolerances(cvode_mem, __CVODE_RTOL__, atol);
    if (check_flag(&flag, "CVodeSVtolerances", 1))
    {   /* Free dynamic memory */
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(noise);
        N_VDestroy_Serial(atol);
        CVodeFree( &cvode_mem );
        gsl_rng_free( rand_gen );
        return 1;
    }
    #else
    /* Set scalar relative and absolute tolerances */
    flag = CVodeSStolerances(cvode_mem, __CVODE_RTOL__, __CVODE_ATOL__);
    if (check_flag(&flag, "CVodeSStolerances", 1))
    {   /* Free dynamic memory */
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(noise);
        CVodeFree( &cvode_mem );
        gsl_rng_free( rand_gen );
        return 1;
    }
    #endif

   
    /* pass params to rhs_func */
    flag = CVodeSetUserData(cvode_mem, data);
    if (check_flag(&flag, "CVodeSetFdata", 1))
    {   /* Free dynamic memory */
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(noise);
        #ifdef __CVODE_VECTOR_ATOL__
        N_VDestroy_Serial(atol);
        #endif
        CVodeFree( &cvode_mem );
        gsl_rng_free( rand_gen );
        return 1;
    }

    /* Choose Sparse or dense method for Jacobian approximation */
    if ( __CVODE_SPARSE__ )
    {
        flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
        if (check_flag(&flag, "CVSpgmr", 1))
        {   /* Free dynamic memory */
            N_VDestroy_Serial(x);
            N_VDestroy_Serial(noise);
            #ifdef __CVODE_VECTOR_ATOL__
            N_VDestroy_Serial(atol);
            #endif
            CVodeFree( &cvode_mem );
            gsl_rng_free( rand_gen );
            return 1;
        }    
    }
    else
    {   /* dense */
        flag = CVDense(cvode_mem, __N_EQNS__); 
        if (check_flag(&flag, "CVDense", 1))
        {   /* Free dynamic memory */
            N_VDestroy_Serial(x);
            N_VDestroy_Serial(noise);
            #ifdef __CVODE_VECTOR_ATOL__
            N_VDestroy_Serial(atol);
            #endif
            CVodeFree( &cvode_mem );
            gsl_rng_free( rand_gen );
            return 1;
        }   
    }
        

    flag = CVodeSetMaxNumSteps(cvode_mem, __CVODE_MAX_NUM_STEPS__);
    if (check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    {   /* Free dynamic memory */
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(noise);
        #ifdef __CVODE_VECTOR_ATOL__
        N_VDestroy_Serial(atol);
        #endif
        CVodeFree( &cvode_mem );
        gsl_rng_free( rand_gen );  
        return 1;
    }


    flag = CVodeSetMaxStep(cvode_mem, __CVODE_MAX_STEP__);
    if (check_flag(&flag, "CVodeSetMaxStep", 1))
    {   /* Free dynamic memory */
        N_VDestroy_Serial(x);
        N_VDestroy_Serial(noise);
        #ifdef __CVODE_VECTOR_ATOL__
        N_VDestroy_Serial(atol);
        #endif
        CVodeFree( &cvode_mem );
        gsl_rng_free( rand_gen );  
        return 1;
    }
    

    /* iterate over output timepoints */
    for ( ii=1;  ii < n_sim_pts;  ii++ )
    {
        #if __SDE__
        if ( seed != 0 )
        {   /* generate state noise */
            flag = generate_state_noise( t_sim[ii-1], t_sim[ii], x, noise, rand_gen, data );
            if ( flag < 0 )
            {   /* Free dynamic memory */
                N_VDestroy_Serial(x);
                N_VDestroy_Serial(noise);
                #ifdef __CVODE_VECTOR_ATOL__
                N_VDestroy_Serial(atol);
                #endif
                CVodeFree( &cvode_mem );
                gsl_rng_free( rand_gen );
                return 1;
            }
        }
        #endif
        
        /* integrate to next timepoint */
        flag = CVode(cvode_mem, t_sim[ii], x, &t, CV_NORMAL);    
        if (check_flag(&flag, "CVode", 1))
        {   /* Free dynamic memory */
            N_VDestroy_Serial(x);
            N_VDestroy_Serial(noise);
            #ifdef __CVODE_VECTOR_ATOL__
            N_VDestroy_Serial(atol);
            #endif
            CVodeFree( &cvode_mem );
            gsl_rng_free( rand_gen );     
            return 1;
        }
       
        #if __SDE__
        if ( seed != 0 )
        {            
            /* add state noise */
            add_state_noise( x, noise );      

            /* reset integrator for the next time point */
            flag = CVodeReInit(cvode_mem, t_sim[ii], x);
            if (check_flag(&flag, "CVodeReInit", 1))
            {   /* Free dynamic memory */
                N_VDestroy_Serial(x);
                N_VDestroy_Serial(noise);
                #ifdef __CVODE_VECTOR_ATOL__
                N_VDestroy_Serial(atol);
                #endif
                CVodeFree( &cvode_mem );
                gsl_rng_free( rand_gen );
                return 1;
            }            
        }
        #endif

        /* copy output from nvector to matlab array */
        for ( jj = 0; jj < __N_EQNS__; ++jj )
        {   x_sim[jj*n_sim_pts + ii] = NV_Ith_S(x,jj);   }
    
    }

    /* Free dynamic memory */
    N_VDestroy_Serial(x);
    N_VDestroy_Serial(noise);
    #ifdef __CVODE_VECTOR_ATOL__
    N_VDestroy_Serial(atol);
    #endif
    CVodeFree( &cvode_mem );
    gsl_rng_free( rand_gen );

    return 0;   
}




/* Add state noise */
int  add_state_noise ( N_Vector x, N_Vector noise )
{
    size_t ii;
    for ( ii = 0; ii < __N_EQNS__; ++ii )
    {
        NV_Ith_S(x, ii) += NV_Ith_S(noise, ii);
    }
    return(0);
}



/* output function for debugging purposes */
void print_value( char * name, realtype value )
{
    mexPrintf( "%s = %lf\n", name, value );
}




/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */
int check_flag( void * flagvalue, char * funcname, int opt )
{
    int *err_flag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL)
    {
        mexPrintf( "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n", funcname );       
        return(1);
    }

    /* Check if flag < 0 */
    else if (opt == 1)
    {
        err_flag = (int *) flagvalue;
        if (*err_flag < 0)
        {
            mexPrintf( "\nSUNDIALS_ERROR: %s() failed with flag = %d\n", funcname, *err_flag );
            return(1);
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL)
    {
        mexPrintf( "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n", funcname );
        return(1);
    }

    return(0);
}



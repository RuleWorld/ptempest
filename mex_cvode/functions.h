#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_


/* heaviside step function approximation constant */
#define __HEAVISIDE_APPROX_CONST__  200

/* heaviside step function ( we use the convention heaviside(0) = 1 ) */
double heaviside        ( double x );
double heaviside_approx ( double x );

/* ramp(x) = x if x >= 0 and 0 otherwise */
double ramp ( double x );

/* hill functions */
double hill     ( double x, double n, double k );
double hilln    ( double x, double n, double k );
double hillnlog ( double x, double max_n, double n, double k );

/* continuous or function */
double continuous_or ( double x, double y );

/* saturation function */
double satur ( double x, double lb, double ub );


#endif

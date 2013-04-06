#include "functions.h"
#include <math.h>


/* user-function definintions */
double heaviside ( double x )
{   return ( x >= 0.0 ? 1.0 : 0.0 );   }

double heaviside_approx ( double x )
{   return (  1.0 / (1.0 + exp( -__HEAVISIDE_APPROX_CONST__*x )) );   }

double ramp ( double x )
{   return ( heaviside(x)*x );   }

double hill ( double x, double n, double k )
{   return ( pow(x,n)/(pow(x,n) + pow(k,n)) );   }

double hilln ( double x, double n, double k )
{   return ( pow(x,n)/(pow(x,n) + pow(k,n))*(1.0 + pow(k,n)) );   }

double hillnlog ( double x, double x_max, double n, double k )
{   return ( x >= 1.0 ? pow(log(x)/log(x_max), n)/(pow(log(x)/log(x_max), n) + pow(k, n)) * (1.0 + pow(k, n)) : 0.0 );   }

double continuous_or ( double x, double y )
{   return ( x*(1.0 - y) + (1.0 - x)*y + x*y );   }

double satur ( double x, double lb, double ub )
{   return ( (x < lb) ? lb : ( (x >= ub) ? ub : x ) );   }


/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk) and Anne Cori (a.cori@imperial.ac.uk), 2012.
  Licence: GPL >=2.
*/


/* 
   This files contains includes for relevant libraries and some global variables and variable types.
*/

#ifndef __COMMON_H
#define __COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* Calls to GNU Scientific Library */
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_pow_int.h>

/* #include <gsl/gsl_blas.h> */
/* #include <gsl/gsl_eigen.h> */
/* #include <gsl/gsl_linalg.h> */
/* #include <gsl/gsl_math.h> */
/* #include <gsl/gsl_matrix.h> */
/* #include <gsl/gsl_permutation.h> */
/* #include <gsl/gsl_sf.h> */
/* #include <gsl/gsl_sf_gamma.h> */
/* #include <gsl/gsl_sort.h> */

/* R headers of R API */
#include <R.h>

#define NEARZERO 0.00000000000000000001 /* 1e-20 */
#define NEARPLUSINF 100000000000000000000.0 /* 10e20 */
#define NEARMINUSINF -100000000000000000000.0 /* -10e20 */
#define TRUE 1
#define FALSE 0


typedef short int bool;




#endif

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector_double.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>
 
#include "mydefn.h"
#include "matrix.h"

double den_nor ( double x, double m, double s2 );
double den_norm_prec ( double x, double m, double prec );
double log_nor (double x, double mu, double s2);
double den_st_prec ( double x, double m, double inv_s2, double df);
double den_st (double x, double mu, double s2, double df);
double log_mvt_den(const gsl_vector *x, const gsl_vector *M, double df, const gsl_matrix *V);
double log_mnv_den(const gsl_vector *x, const gsl_vector *M, const gsl_matrix *V);
double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mn, const gsl_matrix *vr);
double gen_Student_t(gsl_vector *v_t, gsl_matrix *m_Q, int i_n);
double dens_inv_gamma(const double d_x, const double d_nu, const double d_s);
double ln_dens_inv_gamma(const double d_x, const double d_nu, const double d_s);

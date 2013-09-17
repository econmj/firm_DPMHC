#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics_double.h>

#include "matrix.h"
#include "ran.h"
#include "distr.h"
#include "mydefn.h"

void olsg(gsl_vector *y, gsl_matrix *X, gsl_vector *XTy, gsl_matrix *XTX);
void linear_gibbs_beta(const gsl_matrix *XTX, const gsl_vector *XTy,
		       const gsl_vector *b0, const double s2,
		       const gsl_matrix *T0, gsl_vector *draw);
double linear_gibbs_sigma2(const gsl_vector *y, const gsl_matrix *X,
		       const gsl_vector *beta, const double v0,
			   const double s0);
double sample_gibbs_mean(gsl_vector *y, gsl_vector *x,
			 double m, double v2, double s2);
double sample_gibbs_var(gsl_vector *y, gsl_vector *x,
			double beta, double v0, double s0);
gsl_vector* linear_ols_beta(gsl_vector *v_y, gsl_matrix *m_X);


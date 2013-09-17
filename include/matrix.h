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

#include "mydefn.h"


void vec_copy(gsl_vector *x , int s, int n, gsl_vector *y, int t);
double sum(const gsl_vector *y);
double sumr(const gsl_vector *y,const size_t st,const size_t ed);
double mean(const gsl_vector *y);
void mean_matrix_col(const gsl_matrix *A , gsl_vector *sample_mean);
double var(const gsl_vector *y);
void acf(gsl_vector *y, long n, gsl_vector *rho);
void cov_matrix_col(const gsl_matrix *A, gsl_matrix *cov);
void den_interval(const gsl_vector *y,const double a,double * q_a, const double b,double * q_b);
void mat_xtx(const gsl_matrix *x, gsl_matrix *xtx);
void pmat(const gsl_matrix *x);
void pmatr(const gsl_matrix *x, size_t a1,size_t a2,size_t b1, size_t b2);
void pvec(const gsl_vector *x);
void pvecr(const gsl_vector *x,size_t st,size_t ed);
void pvec_int(const gsl_vector_int *x);
void mypause(void);
double dens_est(double x, gsl_vector *data, double h);
double epan( double u);
void inv_pd_sym( gsl_matrix *A );
void cholesky_partial_solve(const gsl_matrix *chol, gsl_vector *y);

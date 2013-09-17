#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "matrix.h"
#include "mydefn.h"

double ran_gamma(gsl_rng *r, double v, double s);
double ran_invgamma(gsl_rng *r,double v,double s);
void normal_gamma_univ(gsl_vector *x,double v, double s,double b_0, double tau);
void normal_invgamma_univ(gsl_vector *x,double v, double s,double b_0, double tau);
size_t ran_multinomial(const gsl_vector *p,const size_t n);
//int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mn,
//	    const gsl_matrix *v, gsl_vector *result);
int ran_mvn(const gsl_vector *mu, const gsl_matrix *V, gsl_vector *draw);
int ran_mvn_alt(const gsl_vector *mu, const gsl_matrix *V, gsl_vector *draw);
int ran_mvt(const gsl_vector *mu, const gsl_matrix *V, double df, gsl_vector *draw);
int ran_wish(const gsl_matrix *m_S0, const int i_v, gsl_matrix *m_Draw);
int ran_inv_wish(const gsl_matrix *m_S0, const int i_v, gsl_matrix *m_Draw);

#ifndef _JENSEN_UTIL_H_
#define _JENSEN_UTIL_H_

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>

//gsl_rng *rng; /* needed to produce global seed for RNG */



void gsl_gnuplot_vector(gsl_vector *v_X);
void dgmpnt_gsl_vec(gsl_vector *v_A,FILE *outputfile);

void dgmpnt_gsl_matrix(gsl_matrix *m_A,FILE *outputfile);

void Quasi_Newton_Hessian(gsl_matrix *m_QuasiHess, 
			  gsl_vector *v_grad_new, gsl_vector *v_x_new,
			  gsl_vector *v_grad_old, gsl_vector *v_x_old);

int hessian( const gsl_vector* x, double (*f)( const gsl_vector*, void* ),
	     void* params, double h, gsl_matrix* H );

gsl_vector *jensen_ran_multv_tdist(const gsl_rng *gsl_r, gsl_vector *v_mu,
				   gsl_matrix *m_Sigma, double d_nu);

double Student_t_kernel(gsl_vector *v_x, gsl_vector *v_mu, 
			gsl_matrix *m_Sigma, double d_nu);

double find_machine_precision(void);

void pnt_gsl_matrix(gsl_matrix *m_A,FILE *outputfile);

#endif  /* _JENSEN_UTIL_H_ */

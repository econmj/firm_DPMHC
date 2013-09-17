#ifndef PS_H_INCLUDED
#define PS_H_INCLUDED

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

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

void posterior_summary(const gsl_matrix *theta, FILE* ofile, long M);
void den_est_file (gsl_vector *y , double Lower , double Upper, int Points , FILE *fp, double h);

void write_posterior_summary(const gsl_matrix *theta, char *s_dirname, long M);

#endif // PS_H_INCLUDED

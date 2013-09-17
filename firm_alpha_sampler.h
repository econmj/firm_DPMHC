#ifndef FIRM_ALPHA_SAMPLER_H_INCLUDED
#define FIRM_ALPHA_SAMPLER_H_INCLUDED

#include "read_mf_data.h"
#include "jensen_util.h"
#include "matrix.h"
#include "DPMN_polya_smplr.h"
#include "DPMHC_smplr.h"
#include "kernel_histogram.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_randist.h>

extern gsl_rng *rng;


int firm_alpha_sampler(struct str_firm_data *a_firm_data, int i_n_firms, struct str_DPMN *ptr_DPMN_alpha);
int firm_alpha_DPMHC_sampler(struct str_firm_data *a_firm_data, int i_n_firms, struct str_DPMHC *ptr_DPMHC_alpha);
int firm_RndEffect_alpha_sampler(struct str_firm_data *a_firm_data, int i_n_firms, double d_mu_a, double d_s2_a);

int firm_alpha_shrinkage(struct str_firm_data *a_firm_data, int i_n_firms, struct str_DPMN *ptr_DPMN_alpha);

#endif // FIRM_ALPHA_SAMPLER_H_INCLUDED

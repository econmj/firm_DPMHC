#ifndef READ_MF_DATA_INCLUDED
#define READ_MF_DATA_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "linear_model.h"
#include "ps.h"
#include "mydefn.h"

struct str_firm_data{
    int i_id;
    int i_ni;
    gsl_vector *v_ret;
    gsl_matrix  *m_factors;
    gsl_matrix *m_XTX;
    gsl_vector *v_XTy;
    double d_s2;
    double d_alpha;
    gsl_vector *v_beta;
    gsl_vector *v_alpha; // alpha_t = mu (1-rho) + rho alpha{t-1} + w e_t
    double d_mu;
    double d_rho;
    double d_w2;
    int i_si;
    gsl_matrix *m_mcmc_draws;
};

 struct str_firm_data* read_mf_data(char *c_filename, int *p_i_n_firms);

 void alloc_mcmc_matrix(struct str_firm_data *a_firm_data, int i_n_firms, int i_draws, int i_factors);

 int write_mf_draws(struct str_firm_data *a_firm_data, int i_n_firms);
 int write_mfs_alpha_draws(struct str_firm_data *a_firm_data, int i_n_firms);


#endif // READ_MF_DATA_INCLUDED

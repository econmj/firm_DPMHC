#ifndef FIRM_LINEAR_MODEL_H_INCLUDED
#define FIRM_LINEAR_MODEL_H_INCLUDED

#include "read_mf_data.h"
#include "linear_model.h"
#include "jensen_util.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>



int firm_linear_beta_s2(struct str_firm_data *a_firm_data, int i_J, gsl_vector *v_beta0, gsl_matrix *m_Lambda0, double d_nu_y0, double d_s2_y0);
int firm_put_MCMC(struct str_firm_data *a_firm_data, int i_J, int i_draw);
int firm_AR_beta_s2(struct str_firm_data *a_firm_data, int i_J, gsl_vector *v_beta0, gsl_matrix *m_Lambda0, double d_nu_y0, double d_s2_y0);
int firm_persist_beta_s2(struct str_firm_data *a_firm_data, int i_J, gsl_vector *v_beta0, gsl_matrix *m_Lambda0, double d_nu_y0, double d_s2_y0);
int firm_linear_alpha_beta_s2(struct str_firm_data *a_firm_data, int i_J, gsl_vector *v_beta0, gsl_matrix *m_Lambda0, double d_nu_y0, double d_s2_y0);

int firm_mu_var_alpha(struct str_firm_data *a_firm_data, int i_J, double *d_mu_a, double *d_s2_a,
                      double d_m0, double d_kappa0, double d_v0, double d_s0);

#endif // FIRM_LINEAR_MODEL_H_INCLUDED

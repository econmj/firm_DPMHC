#ifndef DPMHC_SMPLR_H_INCLUDED
#define DPMHC_SMPLR_H_INCLUDED

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include "ran.h"
#include "distr.h"
#include "matrix.h"
#include "mydefn.h"
#include "read_mf_data.h"

extern gsl_rng *rng;

struct str_DPMHC{
    gsl_vector *v_y;
    gsl_matrix *m_DPtheta;  // K X 3 { mu*_j, xi^*_j, tau^*_j }
    gsl_vector_int *vi_S;
    gsl_vector_int *vi_n;
    int i_K;                // Truncation of mixture
    int i_m;                // # of alive clusters i_m ,= i_K
    double d_DPalpha;
    gsl_vector *v_v;
    gsl_vector *v_w;
    gsl_vector *v_u;
    gsl_vector *v_e;
    int i_ctr;
    gsl_matrix *m_DPmcmc;  // will contain draws of i_K and d_DPalpha
    gsl_matrix *m_pred_den;
    double d_m0;
    double d_s2m;
    double d_A;
    double d_a;
    double d_b;
};

int DPMHC_alloc(struct str_DPMHC *ptr_DPMHC_data, size_t i_T, size_t i_maxK, double d_m0, double d_s2m,
    double d_A, double d_a, double d_b);
int DPMHC_init(struct str_DPMHC *ptr_DPMHC_data, int i_draws);
int DPMHC_v_smplr(struct str_DPMHC *ptr_DPMHC_data);
int DPMHC_u_smplr(struct str_DPMHC *ptr_DPMHC_data);
int DPMHC_K(struct str_DPMHC *ptr_DPMHC_data);
int DPMHC_S_smplr(struct str_DPMHC *ptr_DPMHC_data);
int DPMHC_mu_smplr(struct str_DPMHC *ptr_DPMHC_data);
int DPMHC_tau_smplr(struct str_DPMHC *ptr_DPMHC_data);
int DPMHC_xi_smplr(struct str_DPMHC *ptr_DPMHC_data, int i_J, struct str_firm_data *a_firm_data);

int DPMHC_smplr(struct str_DPMHC *ptr_DPMHC_data,int i_J, struct str_firm_data *a_firm_data);

#endif // DPMHC_SMPLR_H_INCLUDED

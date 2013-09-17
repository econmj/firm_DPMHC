#ifndef DPMN_POLYA_SMPLR_H_INCLUDED
#define DPMN_POLYA_SMPLR_H_INCLUDED

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include "ran.h"
#include "distr.h"
#include "matrix.h"
#include "mydefn.h"

extern gsl_rng *rng;

struct str_DPMN{
    gsl_vector *v_y;
    gsl_matrix *m_DPtheta;
    gsl_vector_int *vi_S;
    gsl_vector_int *vi_n;
    int i_K;
    double d_DPalpha;
    int i_ctr;
    gsl_matrix *m_DPmcmc;  // will contain draws of i_K and d_DPalpha
    gsl_matrix *m_pred_den;
    double d_m0;
    double d_kappa0;
    double d_nu0;
    double d_s20;
    double d_a;
    double d_b;
};

int DPMN_alloc(struct str_DPMN *ptr_DPMN_data, int i_T, int i_maxK, double d_m0, double d_kappa0,
    double d_nu0, double d_s20, double d_a, double d_b);

int DPMN_initialize(struct str_DPMN *ptr_DPMN_data, int i_draws);

int DPMN_initialize_pred_dens(struct str_DPMN *ptr_DPMN_data);

int DPMN_polya_smplr(struct str_DPMN *ptr_DPMN_data);

int DPMN_vS_K_polya_smplr(struct str_DPMN *ptr_DPMN_data);

int DPMN_theta_polya_smplr(struct str_DPMN *ptr_DPMN_data);

int DPMN_DPalpha_smplr(struct str_DPMN *ptr_DPMN_data);

double DPMN_z_unif_prior(double d_alpha, double d_b);

double DPMN_DPalpha_ln_p_k(int i_K, double d_DPalpha, int i_n);

int DPMN_compute_pred_dens(struct str_DPMN *ptr_DPMN_data);

int DPMN_store_draw(struct str_DPMN *ptr_DPMN_data);

int DPMN_write_DPMN_draw(FILE *f_outfile, struct str_DPMN DPMN_data);

#endif // DPMN_POLYA_SMPLR_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

#include "read_mf_data.h"
#include "firm_linear_model.h"
#include "firm_alpha_sampler.h"
#include "firm_DPMN_alpha.h"
#include "DPMN_polya_smplr.h"
#include "DPMHC_smplr.h"
#include "firm_DPMHC.h"
#include "jensen_util.h"
#include "ps.h"

gsl_rng *rng;

int main()
{

    // inputs
  //    char c_filename[1024] = "/home/jensen/fisher/MutualFunds/Data/Jones_Shanken_data.csv";  // data file
char c_filename[1024] = "/home/jensen/fisher/MutualFunds/Data/Jones_Shanken_data.csv"; // data file
    // number of factors
    int i_factors = 4;
    // MCMC inputs
    int l,j,i_cntr = 0;
    int i_burn = 500;
    int i_draws = 100;
    int i_thin = 1;
    int i_ftr = 1;

    // hyperprior parameters to beta_i|sigma^2_i ~ N(v_beta0, sigma^2_i m_Lambda0^-1)
    // and sigma^2_i ~ IG(d_nu_y0/2, d_s2_y0 * d_nu_y0/2) for all i
    double d_nu_y0 = 0;
    double d_s2_y0 = 0;
    gsl_vector *v_beta0 = gsl_vector_alloc(i_factors);
    gsl_vector_set_all(v_beta0,0.0);
    gsl_matrix *m_Lambda0 = gsl_matrix_alloc(i_factors,i_factors);
    gsl_matrix_set_zero(m_Lambda0);

    // DPM storage
    int i_maxK;

    struct str_DPMHC s_DPMHC_data;

    // Base disrtribution hyperparameters eta|lambda^2 ~ N(d_m0, lambda^2/d_kappa0)
    // lambda^2 ~ IG(d_nu0/2, d_s20 * d_nu0/2), and DPMalpha ~ Gamma(d_a,d_b)
    double d_m0 = 0.;
    double d_s2m = 100;
    double d_nu0 = 0.01;
    double d_s20 = 0.01;
    double d_A = 25.0;
    double d_a = 2.0;
    double d_b = 30.0;

    // initialize RNG
    gsl_rng_env_setup();
    rng = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng,(unsigned long int) 874526599);
    fprintf (stdout,"generator type: %s\n", gsl_rng_name (rng));

    // Read in data
    int i_n_firms;
    struct str_firm_data *a_firm_data = read_mf_data(c_filename,&i_n_firms);
    //i_n_firms = 1;

    i_maxK = i_n_firms;

    DPMHC_alloc(&s_DPMHC_data, (size_t) i_n_firms, (size_t) i_maxK, d_m0, d_s2m, d_A,
                d_a, d_b);
    DPMHC_init(&s_DPMHC_data, i_draws);

    FILE *f_DPMN = fopen("DPMN_draws.txt","w");


    alloc_mcmc_matrix(a_firm_data, i_n_firms, i_draws, i_factors);

    for(l=0;l<i_burn+i_draws*i_thin;l++){

       if(l % i_ftr == 0) {
            fprintf(stdout,"\n Iter = %d\n", l);
            fprintf(stdout, "K = %d, DPalpha = %10.8lg\n",
                    s_DPMHC_data.i_K,
                    s_DPMHC_data.d_DPalpha);
            for(j=0;j<s_DPMHC_data.i_K;j++){
                fprintf(stdout, "   mu*(%d)     xi*(%d)          tau*(%d):\n",j,j,j);
                fprintf(stdout, "%12.8lg %12.8lg %12.8g\n",mget(s_DPMHC_data.m_DPtheta,j,0),
                                                mget(s_DPMHC_data.m_DPtheta,j,1),
                                                mget(s_DPMHC_data.m_DPtheta,j,2) );
            }
       }


        // draw beta_i and sigma^2_i and store in a_firm_data structure
        firm_linear_beta_s2(a_firm_data, i_n_firms, v_beta0, m_Lambda0, d_nu_y0, d_s2_y0);
        firm_alpha_DPMHC_sampler(a_firm_data, i_n_firms, &s_DPMHC_data);
        firm_DPMHC(a_firm_data, i_n_firms, &s_DPMHC_data);

        if (l == i_burn){
   //         DPMN_initialize_pred_dens(&s_DPMN_data);

        }

        if ( (l >= i_burn) && (l % i_thin == 0) ){

           firm_put_MCMC(a_firm_data, i_n_firms, i_cntr);
           i_cntr++;
/*
            DPMN_store_draw(&s_DPMN_data); // stores draw of DPalpha and i_K
            DPMN_write_DPMN_draw(f_DPMN,s_DPMN_data);
            DPMN_compute_pred_dens(&s_DPMN_data);
            */


       }

    }

    gsl_vector_view vv_pred_den;
  //  vv_pred_den = gsl_matrix_column( s_DPMN_data.m_pred_den, 1);
   // gsl_vector_scale( &vv_pred_den.vector, 1./i_draws);
/*
            FILE *f_test = fopen("pred_den.txt","w");
            int ii;
            for (ii=0;ii<(s_DPMN_data.m_pred_den)->size1;ii++){
                fprintf(f_test,"%14.8e %14.8e\n",mget(s_DPMN_data.m_pred_den,ii,0),  mget(s_DPMN_data.m_pred_den,ii,1));
            }
            fclose(f_test);
*/
 //   firm_alpha_shrinkage(a_firm_data, i_n_firms, &s_DPMN_data);

    //write_mf_draws(a_firm_data, i_n_firms);
    write_mfs_alpha_draws(a_firm_data,i_n_firms);
 //   posterior_summary(s_DPMN_data.m_DPmcmc,stdout,-1);

    printf("Hello world!\n");
    return 0;
}

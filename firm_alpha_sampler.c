#include "firm_alpha_sampler.h"

int firm_alpha_sampler(struct str_firm_data *a_firm_data, int i_n_firms, struct str_DPMN *ptr_DPMN_alpha){

    int i;
    int i_J = i_n_firms;
    int i_ni;
    int i_factors = (a_firm_data[0].v_beta)->size;
    int i_si;
    double d_s2;
    double d_sum;
    double d_s2_si, d_mu_si;
    double d_mean,d_var;
    double d_alphai;
    gsl_matrix *m_DPtheta = ptr_DPMN_alpha->m_DPtheta;

    gsl_vector *v_ret;
    gsl_matrix *m_factors;
    gsl_vector *v_betai;
    gsl_vector *v_yi;
    gsl_matrix *m_Xi;
    gsl_matrix_view mv_Xi;

    for(i=0;i<i_J;i++){

        i_ni = a_firm_data[i].i_ni;
        //i_si = a_firm_data[i].i_si;
        i_si = vget_int(ptr_DPMN_alpha->vi_S,i);
        d_mu_si = mget(m_DPtheta,i_si,0);
        d_s2_si = mget(m_DPtheta,i_si,1);
        v_betai = a_firm_data[i].v_beta;
        d_s2 = a_firm_data[i].d_s2;
        m_factors = a_firm_data[i].m_factors;


        m_Xi = gsl_matrix_alloc(i_ni,i_factors);
        mv_Xi = gsl_matrix_submatrix(m_factors,0,0,i_ni,i_factors);
        gsl_matrix_memcpy(m_Xi,&mv_Xi.matrix);


        // compute y_i = ret_i - beta_i * X_i
        v_ret = a_firm_data[i].v_ret;
        v_yi = gsl_vector_alloc(i_ni);
        gsl_vector_memcpy(v_yi,v_ret);
        gsl_blas_dgemv(CblasNoTrans, -1.0, m_Xi, v_betai, 1.0, v_yi);

        //posterior draw of alpha_i | y_i, mu_si, s2_si

        d_var = 1./(1./d_s2_si+ i_ni/d_s2);
        //d_var = 1./(0 + i_ni/d_s2);
        d_sum = i_ni*mean(v_yi);
        d_mean = (d_mu_si/d_s2_si+ d_sum/d_s2)*d_var;
        //d_mean = (0 + d_sum/d_s2)*d_var;

        d_alphai = d_mean + gsl_ran_gaussian_ziggurat(rng, sqrt(d_var));
        a_firm_data[i].d_alpha = d_alphai;

        gsl_vector_free(v_yi);
        gsl_matrix_free(m_Xi);
    }

    return 0;
}


int firm_alpha_DPMHC_sampler(struct str_firm_data *a_firm_data, int i_n_firms, struct str_DPMHC *ptr_DPMHC_alpha){

    int i;
    int i_J = i_n_firms;
    int i_ni;
    int i_factors = (a_firm_data[0].v_beta)->size;
    int i_si;
    double d_s2;
    double d_sum;
    double d_mu_si,d_xi_si,d_tau_si;
    double d_mean,d_var;
    double d_ei, d_alphai;
    gsl_matrix *m_DPtheta = ptr_DPMHC_alpha->m_DPtheta;

    gsl_vector *v_ret;
    gsl_matrix *m_factors;
    gsl_vector *v_betai;
    gsl_vector *v_yi;
    gsl_matrix *m_Xi;
    gsl_matrix_view mv_Xi;

    for(i=0;i<i_J;i++){

        i_ni = a_firm_data[i].i_ni;
        //i_si = a_firm_data[i].i_si;
        i_si = vget_int(ptr_DPMHC_alpha->vi_S,i);
        d_mu_si = mget(m_DPtheta,i_si,0);
        d_xi_si = mget(m_DPtheta,i_si,1);
        d_tau_si = mget(m_DPtheta,i_si,2);
        v_betai = a_firm_data[i].v_beta;
        d_s2 = a_firm_data[i].d_s2;
        m_factors = a_firm_data[i].m_factors;


        m_Xi = gsl_matrix_alloc(i_ni,i_factors);
        mv_Xi = gsl_matrix_submatrix(m_factors,0,0,i_ni,i_factors);
        gsl_matrix_memcpy(m_Xi,&mv_Xi.matrix);


        // compute y_i = (ret_i - beta_i * X_i - mu*_si)/(|xi*_si|/tau_si^{1/2})
        v_ret = a_firm_data[i].v_ret;
        v_yi = gsl_vector_alloc(i_ni);
        gsl_vector_memcpy(v_yi,v_ret);
        gsl_blas_dgemv(CblasNoTrans, -1.0, m_Xi, v_betai, 1.0, v_yi);

        gsl_vector_add_constant(v_yi, -d_mu_si);
        gsl_vector_scale(v_yi, 1./(d_xi_si/sqrt(d_tau_si)) ); // ????



        //posterior draw of alpha_i | y_i, mu_si, s2_si

        d_var = 1./(1. + i_ni/( d_s2 * d_tau_si/(d_xi_si * d_xi_si) ) );

        d_sum = i_ni*mean(v_yi);
        d_mean =  d_sum * d_var;

        d_ei = d_mean + gsl_ran_gaussian_ziggurat(rng, sqrt(d_var));
        vset(ptr_DPMHC_alpha->v_e,i,d_ei);

        d_alphai = d_mu_si + (d_xi_si/sqrt(d_tau_si)) * d_ei; //   ????? d_xi_si or |d_xi_si|
        a_firm_data[i].d_alpha = d_alphai;

        gsl_vector_free(v_yi);
        gsl_matrix_free(m_Xi);
    }

    return 0;
}


// Draw each Mutual Fnds alpha_i | r_i, beta_i, s2_i, \propto N(alpha_i | mu_a, s2_a) L(r_i | alpha_i)
// where N(alpha_i | mu_a, s2_a ) is the Random Effects distribution which is known except for mu_a and s2_a
// which will be learned from the cross section of alpha_i, i=1,...,i_n_firms
int firm_RndEffect_alpha_sampler(struct str_firm_data *a_firm_data, int i_n_firms, double d_mu_a, double d_s2_a){

    int i;
    int i_J = i_n_firms;
    int i_ni;
    int i_factors = (a_firm_data[0].v_beta)->size;
    double d_s2;
    double d_sum;
    double d_mean,d_var;
    double d_alphai;

    gsl_vector *v_ret;
    gsl_matrix *m_factors;
    gsl_vector *v_betai;
    gsl_vector *v_yi;
    gsl_matrix *m_Xi;
    gsl_matrix_view mv_Xi;

    for(i=0;i<i_J;i++){

        i_ni = a_firm_data[i].i_ni;
        v_betai = a_firm_data[i].v_beta;
        d_s2 = a_firm_data[i].d_s2;
        m_factors = a_firm_data[i].m_factors;


        m_Xi = gsl_matrix_alloc(i_ni,i_factors);
        mv_Xi = gsl_matrix_submatrix(m_factors,0,0,i_ni,i_factors);
        gsl_matrix_memcpy(m_Xi,&mv_Xi.matrix);


        // compute y_i = ret_i - beta_i * X_i
        v_ret = a_firm_data[i].v_ret;
        v_yi = gsl_vector_alloc(i_ni);
        gsl_vector_memcpy(v_yi,v_ret);
        gsl_blas_dgemv(CblasNoTrans, -1.0, m_Xi, v_betai, 1.0, v_yi);

        //posterior draw of alpha_i | y_i, mu_a, s2_a

        d_var = 1./(1./d_s2_a + i_ni/d_s2);
        d_sum = i_ni*mean(v_yi);
        d_mean = (d_mu_a/d_s2_a + d_sum/d_s2)*d_var;


        d_alphai = d_mean + gsl_ran_gaussian_ziggurat(rng, sqrt(d_var));
        a_firm_data[i].d_alpha = d_alphai;

        gsl_vector_free(v_yi);
        gsl_matrix_free(m_Xi);
    }

    return 0;
}


int firm_alpha_shrinkage(struct str_firm_data *a_firm_data, int i_n_firms, struct str_DPMN *ptr_DPMN_alpha){

    FILE *f_shrink = fopen("alpha_shrink.txt","w");

    int i;
    int i_J = i_n_firms;
    size_t i_ni;
    size_t i_factors = (a_firm_data[0].v_beta)->size;
    double d_alphahat,d_alphabar;
    gsl_vector *v_yi;
    gsl_vector *v_betahat;
    gsl_matrix *m_factors;
    gsl_matrix *m_Xi;
    gsl_matrix_view mv_Xi;
    gsl_vector *v_ones;
    gsl_matrix *m_Z; // Z = [1 X_i]
    gsl_matrix_view mv_Z1;
    gsl_matrix_view mv_ZX;

    gsl_matrix *m_mcmc_draws;
    gsl_vector_view vv_alphai_draws;

    gsl_vector *v_alphahat = gsl_vector_alloc(i_J);


    for (i=0;i<i_J;i++){
        i_ni = (size_t) a_firm_data[i].i_ni;

        v_yi = a_firm_data[i].v_ret;
        m_factors = a_firm_data[i].m_factors;

        m_Xi = gsl_matrix_alloc(i_ni,i_factors);
        mv_Xi = gsl_matrix_submatrix(m_factors,0,0,i_ni,i_factors);
        gsl_matrix_memcpy(m_Xi,&mv_Xi.matrix);

        v_ones = gsl_vector_alloc(i_ni);
        gsl_vector_set_all(v_ones,1.);

        // for Z = [1 X]
        m_Z = gsl_matrix_alloc(i_ni,i_factors + 1);
        mv_Z1 = gsl_matrix_submatrix(m_Z, 0, 0, i_ni, 1);
        mv_ZX = gsl_matrix_submatrix(m_Z, 0, 1, i_ni, i_factors);
        gsl_matrix_set_col(&mv_Z1.matrix, 0, v_ones );
        gsl_matrix_memcpy(&mv_ZX.matrix, m_Xi );

        v_betahat = linear_ols_beta(v_yi,m_Z);
        d_alphahat = vget(v_betahat,0);

        //compute posterior mean of alpha_i draws
        m_mcmc_draws = a_firm_data[i].m_mcmc_draws;
        vv_alphai_draws = gsl_matrix_column(m_mcmc_draws,0);
        d_alphabar = mean(&vv_alphai_draws.vector);

        fprintf(f_shrink,"%14.8lg %14.8lg\n",d_alphahat,d_alphabar);

        vset(v_alphahat,i,d_alphahat);



        gsl_matrix_free(m_Z);
        gsl_vector_free(v_ones);
        gsl_matrix_free(m_Xi);
    }

    smooth_density_vector(v_alphahat);

    fclose(f_shrink);
    return 0;
}

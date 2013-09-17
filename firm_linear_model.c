#include "firm_linear_model.h"


// Draw of (beta_i,sigma^2_i) | { r_{ti} - alpha_i } where pi(beta_i) = N(beta0, sigma^2_i Lambda0^{-1})
// and sigma^2_= IG(nu_y0/2, s2_y0^2 nu_y0/2 ), i = 1,...,i_J where y_i is independent from y_i'
int firm_linear_beta_s2(struct str_firm_data *a_firm_data, int i_J, gsl_vector *v_beta0, gsl_matrix *m_Lambda0, double d_nu_y0, double d_s2_y0)
{
    size_t i_factors = v_beta0->size;
    if ( (i_factors != m_Lambda0->size1) || (i_factors != m_Lambda0->size2) ){
        fprintf(stderr,"Error in firm_linear_beta_s2(): Hyperparameter vector length and matrix row-col must be equal.\n");
        exit(1);
    }

    int i;
    int i_ni;
    gsl_vector *v_ret;
    gsl_vector *v_yi;
    gsl_matrix *m_factors;
    gsl_matrix *m_Xi;
    gsl_matrix_view mv_Xi;

    gsl_vector *v_XiTyi = gsl_vector_alloc(i_factors);
    gsl_matrix *m_XiTXi = gsl_matrix_alloc(i_factors,i_factors);


    double d_alpha;
    double d_s2;
    gsl_vector *v_beta = gsl_vector_alloc(i_factors);

    for(i=0;i<i_J;i++){

        i_ni = a_firm_data[i].i_ni;
        v_ret = a_firm_data[i].v_ret;
        m_factors = a_firm_data[i].m_factors;

        gsl_vector_memcpy(v_beta, a_firm_data[i].v_beta);
        d_alpha = a_firm_data[i].d_alpha;
        d_s2 = a_firm_data[i].d_s2;

        // form covariate matrix X from m_factors
        m_Xi = gsl_matrix_alloc(i_ni,i_factors);
        mv_Xi = gsl_matrix_submatrix(m_factors,0,0,i_ni,i_factors);
        gsl_matrix_memcpy(m_Xi,&mv_Xi.matrix);

        // compute y_i = r_i - alpha
        v_yi = gsl_vector_alloc(i_ni);
        gsl_vector_memcpy(v_yi,v_ret);
        gsl_vector_add_constant(v_yi,-d_alpha);

        //  compute X'X and X'y
        olsg(v_yi, m_Xi, v_XiTyi, m_XiTXi);

        // draw sigma^2_i | r_i - alpha, beta_i
        d_s2 = linear_gibbs_sigma2(v_yi, m_Xi, v_beta, d_nu_y0, d_s2_y0);

        // draw beta_i | r_i - alpha, sigma^2_i
        linear_gibbs_beta(m_XiTXi, v_XiTyi, v_beta0, d_s2, m_Lambda0, v_beta);

        a_firm_data[i].d_s2 = d_s2;
        gsl_vector_memcpy(a_firm_data[i].v_beta,v_beta);

        gsl_vector_free(v_yi);
        gsl_matrix_free(m_Xi);
    }


    gsl_vector_free(v_beta);
    gsl_vector_free(v_XiTyi);
    gsl_matrix_free(m_XiTXi);
    return 0;
}


int firm_put_MCMC(struct str_firm_data *a_firm_data, int i_J, int i_draw)
{
    int i,ii;
    int i_factors;
    double d_alphai,d_s2i;
    gsl_vector *v_betai;
    gsl_matrix *m_mcmc;
    gsl_vector_view vv_mcmc_draw;


    for(i=0;i<i_J;i++){
        d_alphai = a_firm_data[i].d_alpha;
        d_s2i = a_firm_data[i].d_s2;
        v_betai = a_firm_data[i].v_beta;
        i_factors = v_betai->size;
        m_mcmc = a_firm_data[i].m_mcmc_draws;

        vv_mcmc_draw = gsl_matrix_row(m_mcmc,i_draw);

        vset(&vv_mcmc_draw.vector, 0, d_alphai);
        vset(&vv_mcmc_draw.vector, i_factors+1, d_s2i);
        for(ii=0;ii<i_factors;ii++){
            vset( &vv_mcmc_draw.vector, ii+1, v_betai->data[ii]);
        }




    }

    return 0;
}

// alpha_{it} = a_i + rho_i alpha{it-1} + w_i epsilon_{it} where a_i = mu_i ( 1 - rho_i )
// Draw of ( <a_i,rho_i>, w^2_i) | alpha_i where pi(<a_i,rho_i>) = N(<a_i, rho_i>, w^2_i Lambda0^{-1}) I(|rho_i|<1)
// and w^2_i_= IG(nu_alpha/2, s2_alpha^2 nu_alpha/2 ), i = 1,...,i_J
int firm_AR_beta_s2(struct str_firm_data *a_firm_data, int i_J, gsl_vector *v_beta0, gsl_matrix *m_Lambda0, double d_nu0, double d_s2_0)
{
    size_t i_k = v_beta0->size;
    if ( (i_k != 2) || (i_k != m_Lambda0->size1) || (i_k != m_Lambda0->size2) ){
        fprintf(stderr,"Error in firm_AR_beta_s2(): Hyperparameter vector length and matrix row-col must be equal.\n");
        exit(1);
    }

    int i;
    int i_ni;
    gsl_vector *v_ai;
    gsl_vector_view vv_aiL;
    gsl_matrix *m_Xi;


    gsl_vector *v_ones;
    gsl_vector *v_XiTai = gsl_vector_alloc(i_k);
    gsl_matrix *m_XiTXi = gsl_matrix_alloc(i_k,i_k);



    double d_w2_i;
    gsl_vector *v_beta = gsl_vector_alloc(i_k);

    for(i=0;i<i_J;i++){

        i_ni = a_firm_data[i].i_ni;
        v_ai = a_firm_data[i].v_alpha;
        d_w2_i = a_firm_data[i].d_w2;


        // form covariate matrix X from ones and alpha(L)
        m_Xi = gsl_matrix_alloc(i_ni-1,i_k);

        v_ones = gsl_vector_alloc(i_ni-1);
        gsl_vector_set_all(v_ones,1.0);
        gsl_matrix_set_col(m_Xi,0,v_ones);

        vv_aiL = gsl_vector_subvector(v_ai,0,i_ni - 1);
        gsl_matrix_set_col(m_Xi,1,&vv_aiL.vector);

        // set dependent variables alpha_{t+1}
        vv_aiL = gsl_vector_subvector(v_ai,1,i_ni-1);


        //  compute X'X and X'a
        olsg(&vv_aiL.vector, m_Xi, v_XiTai, m_XiTXi);

        // draw beta_i | r_i - alpha, sigma^2_i
        sample_rho:
            linear_gibbs_beta(m_XiTXi, v_XiTai, v_beta0, d_w2_i, m_Lambda0, v_beta);
            if ( fabs(v_beta->data[1]) > 0.9999 )
                goto sample_rho;

        // draw w^2_i
        d_w2_i = linear_gibbs_sigma2(&vv_aiL.vector, m_Xi, v_beta, d_nu0, d_s2_0);


        a_firm_data[i].d_mu = v_beta->data[0]/(1. - v_beta->data[1]); // mu_i = a_i/(1 - rho)
        a_firm_data[i].d_rho = v_beta->data[1];
        a_firm_data[i].d_w2 = d_w2_i;


        gsl_matrix_free(m_Xi);
        gsl_vector_free(v_ones);
    }


    gsl_vector_free(v_beta);
    gsl_vector_free(v_XiTai);
    gsl_matrix_free(m_XiTXi);
    return 0;
}

// Draw of (beta_i,sigma^2_i) | { r_{ti} - alpha_{ti} } where pi(beta_i) = N(beta0, sigma^2_i Lambda0^{-1})
// and sigma^2_= IG(nu_y0/2, s2_y0^2 nu_y0/2 ), i = 1,...,i_J where y_i is independent from y_i'
int firm_persist_beta_s2(struct str_firm_data *a_firm_data, int i_J, gsl_vector *v_beta0, gsl_matrix *m_Lambda0, double d_nu_y0, double d_s2_y0)
{
    size_t i_factors = v_beta0->size;
    if ( (i_factors != m_Lambda0->size1) || (i_factors != m_Lambda0->size2) ){
        fprintf(stderr,"Error in firm_persist_beta_s2(): Hyperparameter vector length and matrix row-col must be equal.\n");
        exit(1);
    }

    int i;
    int i_ni;
    gsl_vector *v_ret;
    gsl_vector *v_yi;
    gsl_matrix *m_factors;
    gsl_matrix *m_Xi;
    gsl_matrix_view mv_Xi;

    gsl_vector *v_XiTyi = gsl_vector_alloc(i_factors);
    gsl_matrix *m_XiTXi = gsl_matrix_alloc(i_factors,i_factors);


    double d_s2;
    gsl_vector *v_beta = gsl_vector_alloc(i_factors);
    gsl_vector *v_alpha;

    for(i=0;i<i_J;i++){

        i_ni = a_firm_data[i].i_ni;
        v_ret = a_firm_data[i].v_ret;
        m_factors = a_firm_data[i].m_factors;

        gsl_vector_memcpy(v_beta, a_firm_data[i].v_beta);
        v_alpha = a_firm_data[i].v_alpha;
        d_s2 = a_firm_data[i].d_s2;

        // form covariate matrix X from m_factors
        m_Xi = gsl_matrix_alloc(i_ni,i_factors);
        mv_Xi = gsl_matrix_submatrix(m_factors,0,0,i_ni,i_factors);
        gsl_matrix_memcpy(m_Xi,&mv_Xi.matrix);

        // compute y_{ti} = r_{ti} - alpha_{ti}
        v_yi = gsl_vector_alloc(i_ni);
        gsl_vector_memcpy(v_yi,v_ret);
        gsl_vector_sub(v_yi,v_alpha);

        //  compute X'X and X'y

        olsg(v_yi, m_Xi, v_XiTyi, m_XiTXi);

        // draw sigma^2_i | r_i - alpha, beta_i
        d_s2 = linear_gibbs_sigma2(v_yi, m_Xi, v_beta, d_nu_y0, d_s2_y0);

        // draw beta_i | r_i - alpha, sigma^2_i
        linear_gibbs_beta(m_XiTXi, v_XiTyi, v_beta0, d_s2, m_Lambda0, v_beta);

        a_firm_data[i].d_s2 = d_s2;
        gsl_vector_memcpy(a_firm_data[i].v_beta,v_beta);

        gsl_vector_free(v_yi);
        gsl_matrix_free(m_Xi);
    }


    gsl_vector_free(v_beta);
    gsl_vector_free(v_XiTyi);
    gsl_matrix_free(m_XiTXi);
    return 0;
}

// Draw of (alpha_i, beta_i,sigma^2_i) | r_{ti} where pi(alpha_i, beta_i) = N(beta0, sigma^2_i Lambda0^{-1})
// and sigma^2_= IG(nu_y0/2, s2_y0^2 nu_y0/2 ), i = 1,...,i_J where y_i is independent from y_i'
int firm_linear_alpha_beta_s2(struct str_firm_data *a_firm_data, int i_J, gsl_vector *v_beta0, gsl_matrix *m_Lambda0, double d_nu_y0, double d_s2_y0)
{
    size_t i_k = v_beta0->size; //i_k = i_factors + 1
    if ( (i_k != m_Lambda0->size1) || (i_k != m_Lambda0->size2) ){
        fprintf(stderr,"Error in firm_linear_alpha_beta_s2(): Hyperparameter vector length and matrix row-col must be equal.\n");
        exit(1);
    }

    int i;
    int i_ni;
    gsl_vector *v_ret;
    gsl_matrix *m_factors;
    gsl_matrix *m_Xi;
    gsl_matrix_view mv_Xi;
    gsl_matrix_view mv_Fi;
    gsl_vector *v_ones;

    gsl_vector *v_XiTyi = gsl_vector_alloc(i_k);
    gsl_matrix *m_XiTXi = gsl_matrix_alloc(i_k,i_k);


    double d_s2;
    gsl_vector *v_beta = gsl_vector_alloc(i_k);
    gsl_vector_view vv_tmp;



    for(i=0;i<i_J;i++){

        i_ni = a_firm_data[i].i_ni;
        v_ret = a_firm_data[i].v_ret;
        m_factors = a_firm_data[i].m_factors;

        // v_beta = (alpha_i, beta_i)'
        v_beta->data[0] = a_firm_data[i].d_alpha;
        vv_tmp = gsl_vector_subvector(v_beta,1,i_k-1);
 //       gsl_vector_memcpy(&vv_tmp.vector, a_firm_data[i].v_beta);
        d_s2 = a_firm_data[i].d_s2;


        // form covariate matrix X = (ones, m_factors)
        m_Xi = gsl_matrix_alloc(i_ni,i_k);

        v_ones = gsl_vector_alloc(i_ni);
        gsl_vector_set_all(v_ones,1.0);
        gsl_matrix_set_col(m_Xi,0,v_ones);

        mv_Fi = gsl_matrix_submatrix(m_factors,0,0,i_ni,i_k-1);
        mv_Xi = gsl_matrix_submatrix(m_Xi,0,1,i_ni,i_k-1);
        gsl_matrix_memcpy(&mv_Xi.matrix,&mv_Fi.matrix);

        //  compute X'X and X'y
        olsg(v_ret, m_Xi, v_XiTyi, m_XiTXi);

        // draw sigma^2_i | r_i - alpha, beta_i
        d_s2 = linear_gibbs_sigma2(v_ret, m_Xi, v_beta, d_nu_y0, d_s2_y0);

        // draw beta_i | r_i - alpha, sigma^2_i
        linear_gibbs_beta(m_XiTXi, v_XiTyi, v_beta0, d_s2, m_Lambda0, v_beta);

        a_firm_data[i].d_s2 = d_s2;
        a_firm_data[i].d_alpha = v_beta->data[0];
        vv_tmp = gsl_vector_subvector(v_beta,1,i_k-1);
        gsl_vector_memcpy(a_firm_data[i].v_beta,&vv_tmp.vector);

        gsl_vector_free(v_ones);
        gsl_matrix_free(m_Xi);
    }


    gsl_vector_free(v_beta);
    gsl_vector_free(v_XiTyi);
    gsl_matrix_free(m_XiTXi);
    return 0;
}


// draw from pi(mu_a, s2_a | v_alpha) where v_alpha = (alpha_1,...,alpha_J) and
// pi( alpha ) ~ N( mu_a, s2_a ) and
// mu_a ~ N(m0, s2_a/kappa0) and s2_a ~ IG(v0/2, s0 v0/2)
int firm_mu_var_alpha(struct str_firm_data *a_firm_data, int i_J, double *d_mu_a, double *d_s2_a,
                      double d_m0, double d_kappa0, double d_v0, double d_s0)
{
    int i;
    gsl_vector *v_alpha = gsl_vector_alloc(i_J);
    gsl_vector *v_zero = gsl_vector_alloc(i_J);
    gsl_vector *v_ones = gsl_vector_alloc(i_J);

    gsl_vector_set_zero(v_zero);
    gsl_vector_set_all(v_ones,1.0);

    for (i=0;i<i_J;i++){
         v_alpha->data[i] = a_firm_data[i].d_alpha;
    }

    *d_mu_a = sample_gibbs_mean(v_alpha, v_zero,
                    d_m0, d_kappa0, *d_s2_a);


    *d_s2_a = sample_gibbs_var(v_alpha, v_ones,
                    *d_mu_a, d_v0, d_s0 * d_v0);


    return 0;
}


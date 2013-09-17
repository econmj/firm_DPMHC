#include "DPMN_polya_smplr.h"


int DPMN_alloc(struct str_DPMN *ptr_DPMN_data, int i_T, int i_maxK, double d_m0, double d_kappa0,
    double d_nu0, double d_s20, double d_a, double d_b){

    ptr_DPMN_data->d_m0 = d_m0;
    ptr_DPMN_data->d_kappa0 = d_kappa0;
    ptr_DPMN_data->d_nu0 = d_nu0;
    ptr_DPMN_data->d_s20 = d_s20;
    ptr_DPMN_data->d_a = d_a;
    ptr_DPMN_data->d_b = d_b;

    ptr_DPMN_data->v_y = gsl_vector_alloc(i_T);
    ptr_DPMN_data->m_DPtheta = gsl_matrix_alloc(i_maxK,2);
    gsl_matrix_set_all(ptr_DPMN_data->m_DPtheta,0.0);
    ptr_DPMN_data->vi_S = gsl_vector_int_alloc(i_T);
    ptr_DPMN_data->vi_n = gsl_vector_int_alloc(i_maxK);

  return 0;
}

int DPMN_initialize(struct str_DPMN *ptr_DPMN_data, int i_draws){
    int i,j;
    int i_T = (ptr_DPMN_data->vi_S)->size;
    if (i_T == 0){
        fprintf(stderr,"Error in DPMN_initialize(): DPMN_alloc() has not been called.\n");
        exit(1);
    }
    //int i_maxK = (ptr_DPMN_data->vi_n)->size;
    gsl_vector *v_probs = gsl_vector_alloc(i_T);
    int i_mndraw;

    ptr_DPMN_data->i_ctr = 0;
    ptr_DPMN_data->m_DPmcmc = gsl_matrix_alloc(i_draws,2); // for draws of i_k and d_DP_alpha

    int i_K = 0;
    double d_m0 = ptr_DPMN_data->d_m0;
    double d_kappa0 = ptr_DPMN_data->d_kappa0;
    double d_nu0 = ptr_DPMN_data->d_nu0;
    double d_s20 = ptr_DPMN_data->d_s20;
    gsl_vector_int *vi_S = ptr_DPMN_data->vi_S;
    gsl_vector_int *vi_n = ptr_DPMN_data->vi_n;
    gsl_vector *v_phi_i = gsl_vector_alloc(2); // phi_i = (eta_i, lambda^2_i)
    gsl_vector_int_set_all(vi_n,0);
    gsl_vector_int_set_all(vi_S,0);

    // draw DP precision parameter d_DPalpha ~ Gamma(a,b)
    double d_DPalpha;
    d_DPalpha = ran_gamma(rng, ptr_DPMN_data->d_a, ptr_DPMN_data->d_b);
    ptr_DPMN_data->d_DPalpha = d_DPalpha;

    // draw phi_i for i=0
    i_K++;
    vset_int(vi_S,0,0);
    vset_int(vi_n,0,1);
    normal_invgamma_univ(v_phi_i, d_nu0, d_nu0*d_s20, d_m0, d_kappa0);
    gsl_matrix_set_row(ptr_DPMN_data->m_DPtheta,0,v_phi_i);


    for(i=1;i<i_T;i++){

        // cluster probs
        vset(v_probs,0, d_DPalpha/(d_DPalpha + i) );
        for(j=1;j<=i;j++){
            vset(v_probs,i,1./(d_DPalpha + i) );
        }

        // if i_mndraw < 0 => phi_i ~ G_0, otherwise phi_i = m_DPtheta[i_mndraw,.]
        i_mndraw = (int)ran_multinomial(v_probs,i)-1;
        if(i_mndraw < 0){
            i_K++;
            vset_int(vi_S,i,i_K-1);
            vset_int(vi_n,i_K-1,1);
            normal_invgamma_univ(v_phi_i, d_nu0, d_nu0*d_s20, d_m0, d_kappa0);
            gsl_matrix_set_row(ptr_DPMN_data->m_DPtheta,i_K-1,v_phi_i);
        } else {
            vset_int(vi_S,i,i_mndraw);
            (vi_n->data[i_mndraw])++;
        }
   }
    ptr_DPMN_data->i_K = i_K;

    gsl_vector_free(v_probs);
    gsl_vector_free(v_phi_i);

  return 0;
}

int DPMN_polya_smplr(struct str_DPMN *ptr_DPMN_data){

    DPMN_vS_K_polya_smplr(ptr_DPMN_data);

    DPMN_theta_polya_smplr(ptr_DPMN_data);

    DPMN_DPalpha_smplr(ptr_DPMN_data);

    return 0;
}



int DPMN_vS_K_polya_smplr(struct str_DPMN *ptr_DPMN_data){

    int i,ii,j;
    int i_n = (ptr_DPMN_data->v_y)->size;
    int i_K =  ptr_DPMN_data->i_K;
    int i_si;
    int i_nj;
    int i_draw;
    double d_yi,d_etaj,d_lambda2j,d_scale,d_sum,d_tmp;
    double d_s21,d_mu1;
    double d_m0 = ptr_DPMN_data->d_m0;
    double d_kappa0 = ptr_DPMN_data->d_kappa0;
    double d_nu0 = ptr_DPMN_data->d_nu0;
    double d_s20 = ptr_DPMN_data->d_s20;
    double d_alpha = ptr_DPMN_data->d_DPalpha;


    gsl_vector *v_y = ptr_DPMN_data->v_y;
    gsl_matrix *m_theta = ptr_DPMN_data->m_DPtheta;
    gsl_vector_int *vi_S = ptr_DPMN_data->vi_S;
    gsl_vector_int *vi_n = ptr_DPMN_data->vi_n;
    gsl_vector *v_probs = gsl_vector_alloc(i_n);
    gsl_vector *v_phi_i = gsl_vector_alloc(2);

    for(i=0;i<i_n;i++)
    {
        //drop (eta_{s_i}, lambda^2_{s_i})
        i_si = vget_int(vi_S,i);
        vi_n->data[ i_si ]--;
        i_nj = vi_n->data[i_si];

        if (i_nj == 0) { // if n_{s_i} = 0 drop theta_{s_i} from m_theta
            if (i_si != i_K-1) { // swap and relabel theta_{s_i} with theta_K, n_{s_i} = n_K and s_i' = K for all i':s_{i'} = K

               gsl_matrix_swap_rows(m_theta, i_si, i_K-1 );
               gsl_vector_int_swap_elements( vi_n, i_si, i_K-1 );

               for(ii=0;ii<i_n;ii++){
                   if (vget_int( vi_S, ii) == i_K-1){
                       vset_int( vi_S, ii, i_si );
                   }
                }

            }
            i_K--;
        }

        // compute probability of each cluster
        d_yi = vget(v_y,i);
        gsl_vector_set_all(v_probs,0.);
        d_sum = 0.;

        for(j=0;j<i_K;j++){
            d_etaj = mget(m_theta,j,0);
            d_lambda2j = mget(m_theta,j,1);
            i_nj = vget_int(vi_n,j);
            d_tmp = i_nj*den_nor( d_yi, d_etaj, d_lambda2j);
            vset(v_probs, j, d_tmp);
            d_sum += d_tmp;
        }

        // prob of a new cluster
        j = i_K;
        d_scale = d_s20 * (1. + d_kappa0)/ d_kappa0;
        d_tmp = d_alpha * den_st( d_yi, d_m0, d_scale, d_nu0);
        vset(v_probs,j,d_tmp);
        d_sum += d_tmp;

        gsl_vector_scale(v_probs,1./d_sum);



        // draw cluster
        i_draw = (int) ran_multinomial( v_probs, i_K);

        if (i_draw == i_K){
            d_mu1 = d_kappa0 * d_m0 + d_yi/(1. + d_kappa0);
            d_s21 = (d_nu0 * d_s20 + d_kappa0/(1. + d_kappa0) * (d_m0 - d_yi) * (d_m0 - d_yi) )/
                    (1. + d_nu0);
            normal_invgamma_univ(v_phi_i, d_nu0 + 1., (d_nu0 + 1.) * d_s21, d_mu1, 1. + d_kappa0);
            gsl_matrix_set_row(m_theta,i_draw,v_phi_i);
            vset_int(vi_S, i, i_draw);
            vset_int(vi_n, i_draw, 1);
            i_K++;

        } else {
            vset_int(vi_S, i, i_draw);
            vi_n->data[i_draw]++;
        }

    }  // end of i loop

    ptr_DPMN_data->i_K = i_K;

    gsl_vector_free(v_probs);
    gsl_vector_free(v_phi_i);
    return 0;
}



int DPMN_theta_polya_smplr(struct str_DPMN *ptr_DPMN_data)
{
   int i,j;
   int i_n = ptr_DPMN_data->v_y->size;
   int i_K = ptr_DPMN_data->i_K;
   int i_nj;
   double d_muj,d_s2j,d_nuj,d_kappaj,d_yhatj,d_ssj;

   double d_m0 = ptr_DPMN_data->d_m0;
   double d_kappa0 = ptr_DPMN_data->d_kappa0;
   double d_nu0 = ptr_DPMN_data->d_nu0;
   double d_s20 = ptr_DPMN_data->d_s20;

   gsl_vector *v_y = ptr_DPMN_data->v_y;
   gsl_vector_int *vi_S = ptr_DPMN_data->vi_S;
   gsl_vector_int *vi_n = ptr_DPMN_data->vi_n;
   gsl_matrix *m_theta = ptr_DPMN_data->m_DPtheta;

   gsl_vector_view vv_theta_j;

 //  printf("\ni_K = %d\n",i_K);

   for(j=0;j<i_K;j++){

        d_yhatj = 0.;
        i_nj = 0;
        for(i=0;i<i_n;i++){
            if( vget_int(vi_S,i) == j ){
                d_yhatj += vget(v_y,i);
                i_nj++;
            }
        }
        if (vget_int(vi_n,j) != i_nj){
            fprintf(stderr,"Error in  DPMN_theta_polya_smplr(): vi_n[%d] does not equal i_nj\n", j);
            exit(1);
        }
        d_yhatj /= i_nj;
        d_nuj = d_nu0 + i_nj;
        d_kappaj = d_kappa0 + i_nj;

        d_ssj = 0;
        for(i=0;i<i_n;i++){
            if (vget_int(vi_S,i) == j ){
                d_ssj += (vget(v_y,i) - d_yhatj) * (vget(v_y,i) - d_yhatj);
            }
        }

        d_muj = (d_kappa0 * d_m0 + i_nj * d_yhatj)/d_kappaj;
        d_s2j = (d_nu0 * d_s20
                 + d_ssj
                 + i_nj * d_kappa0 * (d_m0 - d_yhatj) * (d_m0 -d_yhatj)/d_kappaj)
                 /d_nuj;

        vv_theta_j = gsl_matrix_row(m_theta,j);
        normal_invgamma_univ( &vv_theta_j.vector, d_nuj, d_s2j * d_nuj, d_muj, d_kappaj );
       // printf("%d: eta = %g lambda^2 = %g\n",j, mget(m_theta,j,0), mget(m_theta,j,1) );
     }

    return 0;
}

int DPMN_DPalpha_smplr(struct str_DPMN *ptr_DPMN_data){

    int i_n = (ptr_DPMN_data->v_y)->size;
    int i_K = ptr_DPMN_data->i_K;
    double d_a = ptr_DPMN_data->d_a;
    double d_b = ptr_DPMN_data->d_b;
    double d_alpha = ptr_DPMN_data->d_DPalpha;

    if (d_a > 0){
        double d_eta = gsl_ran_beta( rng, d_alpha + 1., i_n);
        double d_lneta = log(d_eta);
        double d_prob = (d_a + i_K - 1.)/(i_n * (d_b - d_lneta) + d_a + i_K - 1.);

        if (gsl_rng_uniform(rng) <= d_prob){
            ptr_DPMN_data->d_DPalpha = ran_gamma(rng, d_a + i_K, d_b - d_lneta);
        } else {
            ptr_DPMN_data->d_DPalpha = ran_gamma(rng, d_a + i_K - 1., d_b - d_lneta);
        }
    } else if ((d_a == 0) && (d_b > 0)) {
        double d_sig = 1.0;
        double d_alpha_p = d_alpha + gsl_ran_gaussian_ziggurat(rng,d_sig);
        double d_mh = log(DPMN_z_unif_prior(d_alpha_p, d_b))
                        - log(DPMN_z_unif_prior(d_alpha,d_b));
        d_mh += DPMN_DPalpha_ln_p_k(i_K, d_alpha_p, i_n) - DPMN_DPalpha_ln_p_k(i_K, d_alpha, i_n);

        if (d_mh >= 0.0 || gsl_rng_uniform(rng) <= exp(d_mh)){
            ptr_DPMN_data->d_DPalpha = d_alpha_p;
        }
    }

    return 0;
}

/* Prior for alpha given z = b/(b+alpha) ~ U(0,1) */
double DPMN_z_unif_prior(double d_alpha, double d_b)
{
    double d_wrk = d_alpha + d_b;
    return d_b/(d_wrk * d_wrk);
}


/********************************************************************************
log of the Escobar and West (1994) sampling distribution for the number of mixture
clusters k in a DP prior DP(alpha) having i_N observations
*********************************************************************************/
double DPMN_DPalpha_ln_p_k(int i_K, double d_DPalpha, int i_n)
{
  double d_pk = (double) i_K * log(d_DPalpha);
  d_pk += gsl_sf_lngamma(d_DPalpha) - gsl_sf_lngamma(d_DPalpha + i_n);

  return d_pk;
}




int DPMN_initialize_pred_dens(struct str_DPMN *ptr_DPMN_data)
{
    int i;
    int i_nsteps = 1000;

    gsl_vector *v_y = ptr_DPMN_data->v_y;

    ptr_DPMN_data->m_pred_den = gsl_matrix_alloc(i_nsteps,2);
    gsl_matrix *m_pred_den = ptr_DPMN_data->m_pred_den;
    gsl_matrix_set_zero(m_pred_den);

    double d_std = sqrt(var(v_y));
    double d_min = gsl_vector_min(v_y) - d_std;
    double d_max = gsl_vector_max(v_y) + d_std;

    double d_step = (d_max - d_min)/(double) i_nsteps;
    for(i=0;i<i_nsteps;i++){
        mset(m_pred_den,i,0, d_min + (double)i * d_step);
    }

   return 0;
}


int DPMN_compute_pred_dens(struct str_DPMN *ptr_DPMN_data)
{
    int i,j;
    double d_tmp,d_scale,d_ynp1,d_etaj,d_lambda2j;

    double d_m0 = ptr_DPMN_data->d_m0;
    double d_kappa0 = ptr_DPMN_data->d_kappa0;
    double d_nu0 = ptr_DPMN_data->d_nu0;
    double d_s20 = ptr_DPMN_data->d_s20;
    double d_alpha = ptr_DPMN_data->d_DPalpha;

    gsl_matrix *m_DPtheta = ptr_DPMN_data->m_DPtheta;
    gsl_matrix *m_pred_den = ptr_DPMN_data->m_pred_den;
    int i_obs = m_pred_den->size1;
    int i_n = ptr_DPMN_data->v_y->size;
    int i_K = ptr_DPMN_data->i_K;
    gsl_vector_int *vi_n = ptr_DPMN_data->vi_n;
    int i_nj;

    d_scale = d_s20 * (1. + d_kappa0)/ d_kappa0;
    for(i=0;i<i_obs;i++){

        d_ynp1 = mget(m_pred_den,i,0);

        d_tmp = d_alpha/(d_alpha + i_n) * den_st( d_ynp1, d_m0, d_scale, d_nu0);

        for(j=0;j<i_K;j++){

            i_nj = vget_int(vi_n,j);
            d_etaj = mget(m_DPtheta,j,0);
            d_lambda2j = mget(m_DPtheta,j,1);
            d_tmp += i_nj/(d_alpha + i_n) * gsl_ran_gaussian_pdf( d_ynp1 - d_etaj, sqrt(d_lambda2j) );

        }

        mset(m_pred_den,i,1, mget(m_pred_den,i,1) + d_tmp);
    }

    return 0;
}

int DPMN_store_draw(struct str_DPMN *ptr_DPMN_data){

    int i_K = ptr_DPMN_data->i_K;
    double d_DPalpha = ptr_DPMN_data->d_DPalpha;

    gsl_matrix *m_DPmcmc = ptr_DPMN_data->m_DPmcmc;
    int i_ctr = ptr_DPMN_data->i_ctr;

    if ((i_ctr >= 0) &&(i_ctr < m_DPmcmc->size1)){
        mset(m_DPmcmc,i_ctr,0,i_K);
        mset(m_DPmcmc,i_ctr,1,d_DPalpha);
    } else {
    return 1;
    }

    (ptr_DPMN_data->i_ctr)++;
    return 0;
}

int DPMN_write_DPMN_draw(FILE *f_outfile, struct str_DPMN DPMN_data){

    int j,i_K,i_nj;
    double d_etaj,d_lambda2j,d_DPalpha;
    gsl_vector_int *vi_n = DPMN_data.vi_n;
    gsl_matrix *m_DPtheta = DPMN_data.m_DPtheta;

    i_K = DPMN_data.i_K;
    d_DPalpha = DPMN_data.d_DPalpha;

    fprintf(f_outfile,"%4.0d %14.8lg\n",i_K,d_DPalpha);

    for(j=0;j<i_K;j++){

        i_nj = vget_int(vi_n,j);
        d_etaj = mget(m_DPtheta,j,0);
        d_lambda2j = mget(m_DPtheta,j,1);

        fprintf(f_outfile, "%8.0d %14.8lg %14.8lg\n",i_nj, d_etaj, d_lambda2j);
    }



    return 0;
}

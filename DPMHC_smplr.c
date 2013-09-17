#include "DPMHC_smplr.h"

int DPMHC_alloc(struct str_DPMHC *ptr_DPMHC_data, size_t i_T, size_t i_maxK, double d_m0, double d_s2m,
    double d_A, double d_a, double d_b)
{
    ptr_DPMHC_data->d_m0 = d_m0;
    ptr_DPMHC_data->d_s2m = d_s2m;
    ptr_DPMHC_data->d_A = d_A;
    ptr_DPMHC_data->d_a = d_a;
    ptr_DPMHC_data->d_b = d_b;

    ptr_DPMHC_data->v_y = gsl_vector_alloc(i_T);
    ptr_DPMHC_data->m_DPtheta = gsl_matrix_alloc(i_maxK,3);
    gsl_matrix_set_all(ptr_DPMHC_data->m_DPtheta,0.0);
    ptr_DPMHC_data->vi_S = gsl_vector_int_alloc(i_T);
    ptr_DPMHC_data->v_u = gsl_vector_alloc(i_T);
    ptr_DPMHC_data->v_e = gsl_vector_alloc(i_T);
    ptr_DPMHC_data->vi_n = gsl_vector_int_alloc(i_maxK);
    ptr_DPMHC_data->v_v = gsl_vector_alloc(i_maxK);
    ptr_DPMHC_data->v_w = gsl_vector_alloc(i_maxK);


    return 0;
}


int DPMHC_init(struct str_DPMHC *ptr_DPMHC_data, int i_draws){

    int j;
    int i_T = (ptr_DPMHC_data->vi_S)->size;
    if (i_T == 0){
        fprintf(stderr,"Error in DPMHC_init(): DPMHC_alloc() has not been called.\n");
        exit(1);
    }

    int i_K;
    gsl_matrix *m_DPtheta = ptr_DPMHC_data->m_DPtheta;

    ptr_DPMHC_data->m_DPmcmc = gsl_matrix_alloc(i_draws, 2); // for draw of i_K and d_DPalpha

    // initialize slice truction to K = 4 and one alive cluster, i_m = 1
    i_K = ptr_DPMHC_data->i_K = 4;
    ptr_DPMHC_data->i_m = 1;
    gsl_vector_int_set_all(ptr_DPMHC_data->vi_S,0);
    vset_int(ptr_DPMHC_data->vi_n,0,i_T);

    // draw DP precision parameter d_DPalpha ~ Gamma(a,b)
    double d_DPalpha;
    d_DPalpha = ran_gamma(rng, ptr_DPMHC_data->d_a, ptr_DPMHC_data->d_b);
    ptr_DPMHC_data->d_DPalpha = d_DPalpha;

    // Draw initial mixture locations for K clusters
    for(j = 0; j < i_K; j++){

        mset(m_DPtheta, j, 0,
                ptr_DPMHC_data->d_m0 + gsl_ran_gaussian_ziggurat(rng, sqrt(ptr_DPMHC_data->d_s2m)));

        mset(m_DPtheta, j, 1,
                gsl_ran_gaussian_ziggurat(rng, ptr_DPMHC_data->d_A));

        mset(m_DPtheta, j, 2,
                gsl_ran_gamma(rng, 0.5, 0.5) );

    }

    return 0;
}


int DPMHC_smplr(struct str_DPMHC *ptr_DPMHC_data, int i_J, struct str_firm_data *a_firm_data)
{
    DPMHC_v_smplr(ptr_DPMHC_data);
    DPMHC_u_smplr(ptr_DPMHC_data);
    DPMHC_K(ptr_DPMHC_data);

    DPMHC_S_smplr(ptr_DPMHC_data);
    DPMHC_mu_smplr(ptr_DPMHC_data);
    DPMHC_tau_smplr(ptr_DPMHC_data);

    DPMHC_xi_smplr(ptr_DPMHC_data, i_J, a_firm_data);

    return 0;
}


int DPMHC_v_smplr(struct str_DPMHC *ptr_DPMHC_data)
{
    gsl_vector *v_v = ptr_DPMHC_data->v_v;
    gsl_vector *v_w = ptr_DPMHC_data->v_w;
    gsl_vector_int *vi_S = ptr_DPMHC_data->vi_S;
    int i_K = ptr_DPMHC_data->i_K;
    double d_DPalpha = ptr_DPMHC_data->d_DPalpha;

    size_t i_T = vi_S->size;
    int i,j,i_si;
    double d_vj,d_prod;
    //   printf("inside sample_v K=%d\n",K);
    gsl_vector *v_a = gsl_vector_alloc ( (size_t) i_K );
    gsl_vector *v_b = gsl_vector_alloc ( (size_t) i_K );

    gsl_vector_set_all ( v_a, 1.0 );
    gsl_vector_set_all ( v_b, d_DPalpha);



    for(i=0;i<i_T;i++){
        i_si = vget_int(vi_S,i);
        (v_a->data[ i_si ]) += 1.0;
        for(j = 0; j < i_si;j++){
            (v_b->data[ j ]) += 1.0;
        }
    }

   /* pvec(a); */
   /* pvec(b); */


    /* take draws and form w */
    for(j=0;j<i_K;j++){
        d_vj = gsl_ran_beta ( rng , vget(v_a,j) , vget(v_b,j) );
        vset( v_v, j, d_vj );
        /* w_j */
        if( j == 0 ){
            vset(v_w,j, d_vj );
            d_prod = (1.0-d_vj);
        }else{
            vset(v_w,j, d_vj*d_prod );
            d_prod *= (1.0-d_vj);
        }
    }

    gsl_vector_free (v_a);
    gsl_vector_free (v_b);

    return 0;
}


int DPMHC_u_smplr(struct str_DPMHC *ptr_DPMHC_data)
{
    gsl_vector_int *vi_S = ptr_DPMHC_data->vi_S;
    gsl_vector *v_w = ptr_DPMHC_data->v_w;
    gsl_vector *v_u = ptr_DPMHC_data->v_u;

    size_t i_T = vi_S->size;
    int i,i_si;
    double d_ui,d_w_si;

    for(i=0;i<i_T;i++){
        i_si = vget_int( vi_S, i);
        d_w_si = vget( v_w, i_si);
        dw:
            d_ui = d_w_si * gsl_rng_uniform (rng);

            if( d_ui == 0.0 ){
            printf("\n\n ******** u_i = 0.0 \n");
            goto dw;
            }

        vset( v_u, i, d_ui);
    }

    return 0;
}

int DPMHC_K(struct str_DPMHC *ptr_DPMHC_data)
{
    int i_K = ptr_DPMHC_data->i_K;
    gsl_vector *v_u = ptr_DPMHC_data->v_u;
    gsl_vector *v_v = ptr_DPMHC_data->v_v;
    gsl_vector *v_w  = ptr_DPMHC_data->v_w;
    gsl_matrix *m_DPtheta = ptr_DPMHC_data->m_DPtheta;
    double d_DPalpha = ptr_DPMHC_data->d_DPalpha;

    int K_tmp, K_new,j;
    double a,v_j,w_j,csum,min_u;
    //gsl_vector_view theta_j;

  //int k_asset_number = P -> size1; /* number of assets in model */

    K_tmp = i_K;
    min_u = gsl_vector_min ( v_u );
    a = 1.0 - min_u;

    if( a == 1.0 )
        printf("**********min_u = %g *************\n",min_u);

    csum = 0.0;
    j=0;

    while ( csum <= a ){

        /* check if new v_j,w_j and theta_j should be generated */
        if( j >= K_tmp ){

            v_j = gsl_ran_beta ( rng , 1.0, d_DPalpha );
            vset( v_v, j, v_j);

            w_j = v_j * (vget( v_w, j-1 )/vget(v_v,j-1))*(1.0-vget(v_v,j-1));
            vset( v_w, j, w_j);

        /* generate new mu, xi, tau from prior G_0 */
            mset(m_DPtheta, j, 0,
                ptr_DPMHC_data->d_m0 + gsl_ran_gaussian_ziggurat(rng, sqrt(ptr_DPMHC_data->d_s2m)));

            mset(m_DPtheta, j, 1,
                gsl_ran_gaussian_ziggurat(rng, ptr_DPMHC_data->d_A));

            mset(m_DPtheta, j, 2,
                gsl_ran_gamma(rng, 0.5, 0.5) );
        }

        csum += vget(v_w,j);
        K_new = j + 1;
        j++;
    }

    ptr_DPMHC_data->i_K = K_new;

    return 0;
}

int DPMHC_S_smplr(struct str_DPMHC *ptr_DPMHC_data)
{
    int i_K = ptr_DPMHC_data->i_K;
    gsl_vector *v_y = ptr_DPMHC_data->v_y;
    gsl_vector *v_w = ptr_DPMHC_data->v_w;
    gsl_vector *v_u = ptr_DPMHC_data->v_u;
    gsl_vector_int *vi_S = ptr_DPMHC_data->vi_S;
    gsl_vector_int *vi_n = ptr_DPMHC_data->vi_n;
    gsl_matrix *m_DPtheta = ptr_DPMHC_data->m_DPtheta;




    size_t i_n=v_y->size;
    gsl_vector *v_p = gsl_vector_calloc (i_K);
    gsl_vector_int *vi_cnt = gsl_vector_int_calloc (i_K);
    int i,j,st,ism;
    double sm,dn;

    double d_yi,d_muj,d_xij,d_tauj;

    gsl_vector_int_set_zero ( vi_cnt );

    gsl_vector_int_set_zero ( vi_n );

    for(i=0;i<i_n;i++){

    //printf("i=%d\n",i);

        d_yi = v_y->data[i];


        sm=0.0;
        for(j=0;j<i_K;j++){

            if( vget(v_w,j) <= vget(v_u,i) ){
                vset(v_p,j,0.0);
        //      printf("w_j,u_i %g %g",vget(w,j),vget(u,i));
            }else{
        //  printf("K=%d\n",i_K);
        //        dn = den_norm_prec ( vget(y,i), mget(theta,j,0), mget(theta,j,1) );

                d_muj = mget(m_DPtheta,j,0);
                d_xij = mget(m_DPtheta,j,1);
                d_tauj = mget(m_DPtheta,j,2);

                dn = exp( log_nor(d_yi, d_muj, fabs(d_xij) * fabs(d_xij) / d_tauj  ));
                vset(v_p,j,dn);
                sm+=dn;
        //printf("p>0\n");
            }
        }

        /* standardize */
        gsl_vector_scale( v_p, 1.0/sm );

      //      pvec(p);

        /* take draw */
        st = (int)ran_multinomial( v_p, i_K-1 );
        vset_int(vi_S,i,st);
        vset_int(vi_cnt,st,1);
        (vi_n->data[st])++;
    }



    /* find number of clusters with at least one observation assigned to it */
    ism=0;
    for(j=0;j<i_K;j++){
        if( vget_int(vi_cnt,j) == 1)ism++;
    }

    ptr_DPMHC_data->i_m = ism;

    gsl_vector_free (v_p);
    gsl_vector_int_free (vi_cnt);


    return 0;
}



int DPMHC_mu_smplr(struct str_DPMHC *ptr_DPMHC_data)
{
    int i,j;
    int i_n = ptr_DPMHC_data->v_y->size;
    int i_K = ptr_DPMHC_data->i_K;
    int i_nj;
    double d_muj,d_s2j,d_yhatj;

    double d_xij,d_tauj;

    double d_m0 = ptr_DPMHC_data->d_m0;
    double d_s2m = ptr_DPMHC_data->d_s2m;

   gsl_vector *v_y = ptr_DPMHC_data->v_y;
   gsl_vector_int *vi_S = ptr_DPMHC_data->vi_S;
   gsl_vector_int *vi_n = ptr_DPMHC_data->vi_n;
   gsl_matrix *m_theta = ptr_DPMHC_data->m_DPtheta;


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
        d_xij = mget(m_theta,j,1);
        d_tauj = mget(m_theta,j,2);
        d_yhatj /= (d_xij * d_xij / d_tauj);

        d_s2j = 1./( 1./d_s2m + i_nj/ (d_xij * d_xij / d_tauj) );
        d_muj = d_s2j * ( d_m0/d_s2m + d_yhatj );


        d_muj = d_muj + gsl_ran_gaussian_ziggurat(rng, sqrt(d_s2j));

        mset(m_theta,j,0, d_muj);


       // printf("%d: eta = %g lambda^2 = %g\n",j, mget(m_theta,j,0), mget(m_theta,j,1) );
     }





    return 0;
}


int DPMHC_tau_smplr(struct str_DPMHC *ptr_DPMHC_data)
{
    int i,j;
    int i_n = ptr_DPMHC_data->v_y->size;
    int i_K = ptr_DPMHC_data->i_K;
    int i_nj;
    double d_muj,d_yi;
    double d_yhat;

    double d_xij,d_tauj;



    gsl_vector *v_y = ptr_DPMHC_data->v_y;
    gsl_vector_int *vi_S = ptr_DPMHC_data->vi_S;
    gsl_vector_int *vi_n = ptr_DPMHC_data->vi_n;
    gsl_matrix *m_theta = ptr_DPMHC_data->m_DPtheta;


 //  printf("\ni_K = %d\n",i_K);
    for(j=0;j<i_K;j++){

        d_muj = mget(m_theta,j, 0);
        d_xij = mget(m_theta,j, 1);

        d_yhat = 0.;
        i_nj = 0;
        for(i=0;i<i_n;i++){
            if( vget_int(vi_S,i) == j ){
                d_yi = vget(v_y,i);
                d_yhat += (d_yi/fabs(d_xij) - d_muj/fabs(d_xij)) * (d_yi/fabs(d_xij) - d_muj/fabs(d_xij));
                i_nj++;
            }
        }
        if (vget_int(vi_n,j) != i_nj){
            fprintf(stderr,"Error in  DPMN_tau_smplr(): vi_n[%d] does not equal i_nj\n", j);
            exit(1);
        }

        d_tauj = gsl_ran_gamma(rng, 0.5 + (double)i_nj/2.0, 0.5 + d_yhat/2.0);

        mset(m_theta,j, 2, d_tauj);


       // printf("%d: eta = %g lambda^2 = %g\n",j, mget(m_theta,j,0), mget(m_theta,j,1) );
     }





    return 0;
}

int DPMHC_xi_smplr(struct str_DPMHC *ptr_DPMHC_data, int i_J, struct str_firm_data *a_firm_data)
{
    int j,i;
    int i_K = ptr_DPMHC_data->i_K;
    int i_n = ptr_DPMHC_data->v_y->size; // is this the same as i_J???
    if (i_n != i_J){
        fprintf(stderr,"Error in  DPMN_xi_smplr(): DPMHC.v_y length does not equal i_J\n");
        exit(1);
    }
    double d_sumT_si;

    gsl_vector_int *vi_S = ptr_DPMHC_data->vi_S;
    gsl_matrix *m_theta = ptr_DPMHC_data->m_DPtheta;
    double d_A = ptr_DPMHC_data->d_A;
    double d_mu_si, d_tau_si, d_ei;
    double d_mean_j, d_var_j;
    double d_xi_j;

    int i_Ti;
    int i_factors = (a_firm_data[0].v_beta)->size;
    gsl_vector *v_ret;
    gsl_matrix *m_factors;
    gsl_vector *v_betai;
    gsl_vector *v_rstar_i;
    gsl_matrix *m_Xi;
    gsl_matrix_view mv_Xi;
    double d_rstar_j;
    double d_s2i;


    for(j=0;j<i_K;j++){

        d_mu_si = mget(m_theta,j,0);
        d_tau_si = mget(m_theta,j,2);


        d_rstar_j = 0.;
        d_sumT_si = 0;
        for(i=0;i<i_J;i++){
            if( vget_int(vi_S,i) == j ){

                d_ei = vget(ptr_DPMHC_data->v_e,i);

                m_factors = a_firm_data[i].m_factors;
                i_Ti = a_firm_data[i].i_ni;
                d_s2i = a_firm_data[i].d_s2;
                m_Xi = gsl_matrix_alloc(i_Ti,i_factors);
                mv_Xi = gsl_matrix_submatrix(m_factors,0,0,i_Ti,i_factors);
                gsl_matrix_memcpy(m_Xi,&mv_Xi.matrix);

                v_betai = a_firm_data[i].v_beta;
                v_ret = a_firm_data[i].v_ret;
                v_rstar_i = gsl_vector_alloc(i_Ti);
                gsl_vector_memcpy(v_rstar_i,v_ret);
                gsl_blas_dgemv(CblasNoTrans, -1.0, m_Xi, v_betai, 1.0, v_rstar_i);

                gsl_vector_add_constant(v_rstar_i, -d_mu_si);
                gsl_vector_scale(v_rstar_i, 1./(sqrt(d_tau_si) * d_ei) );

                d_rstar_j += sum(v_rstar_i);

                d_sumT_si += i_Ti/(d_s2i/(d_tau_si * d_ei * d_ei) );

                gsl_matrix_free(m_Xi);
                gsl_vector_free(v_rstar_i);
            }
        }

        d_var_j = 1./( 1./(d_A * d_A) + d_sumT_si);
        d_mean_j = d_rstar_j * d_var_j;

        d_xi_j = d_mean_j + gsl_ran_gaussian_ziggurat(rng, sqrt(d_var_j) );

        mset(m_theta, j, 1, d_xi_j);

       // printf("%d: eta = %g lambda^2 = %g\n",j, mget(m_theta,j,0), mget(m_theta,j,1) );
     }


    return 0;
}

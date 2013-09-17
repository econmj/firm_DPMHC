#include "read_mf_data.h"

struct str_firm_data* read_mf_data(char *c_filename, int *p_i_n_firms){

    FILE *ptr_file = fopen(c_filename,"r");
    if (!ptr_file){
        printf("ERROR: File, %s not there.\n",c_filename);
        exit(1);
    }

    char buf[1024];
    char *tok;
    int i=0,j;
    int i_ncol=10;
    int i_n_firms=0;
    int i_id,i_firmid=0,i_firm0;

    gsl_vector *v_y;
    gsl_matrix *m_X;
    gsl_vector *v_XTy;
    gsl_matrix *m_XTX;

    // first line of ptr_file header
    fgets(buf,sizeof(buf),ptr_file);

    while (fgets(buf,sizeof(buf),ptr_file)){
        tok = strtok(buf,",");
        i_id = atoi(tok);

        if(i == 0){
            i_firm0 = i_id;  // first firm id
        }

        if (i_id != i_firmid){
            i_firmid = i_id;
            i_n_firms++;
        }

        i++;
    }

    double d_tmp;
    gsl_matrix *m_tmp = gsl_matrix_alloc(i,i_ncol - 2);
    gsl_vector_view vv_tmp;
    gsl_matrix_view mv_tmp;

    *p_i_n_firms = i_n_firms;
    struct str_firm_data *p_str_firm_data = (struct str_firm_data*) malloc(i_n_firms * sizeof(struct str_firm_data));

    rewind(ptr_file);
    /*
    fclose(ptr_file);
    ptr_file = fopen(c_filename,"r");
    */
    fgets(buf,sizeof(buf),ptr_file); // read header

    i=0;
    i_n_firms=0;
    i_firmid=i_firm0;

    while (fgets(buf,sizeof(buf),ptr_file)){

        tok = strtok(buf,",");
        i_id = atoi(tok);

        if ( (i_id != i_firmid) && (i != 0) ){

            p_str_firm_data[i_n_firms].i_id = i_firmid;
            p_str_firm_data[i_n_firms].i_ni = i;
            p_str_firm_data[i_n_firms].i_si = 0;
            p_str_firm_data[i_n_firms].v_ret = gsl_vector_alloc(i);
            p_str_firm_data[i_n_firms].m_factors = gsl_matrix_alloc(i,i_ncol - 3);
            p_str_firm_data[i_n_firms].v_XTy = gsl_vector_alloc(i_ncol-3);
            p_str_firm_data[i_n_firms].m_XTX = gsl_matrix_alloc(i_ncol-3,i_ncol-3);
            p_str_firm_data[i_n_firms].v_alpha = gsl_vector_alloc(i);

            p_str_firm_data[i_n_firms].d_rho = 0.0;
            p_str_firm_data[i_n_firms].d_mu = 0.0;
            p_str_firm_data[i_n_firms].d_w2 = 1.0;

            vv_tmp = gsl_matrix_subcolumn(m_tmp,0,0,i);
            gsl_vector_scale( &vv_tmp.vector, 1200);  // convert monthly returns into yearly percentage returns
            gsl_vector_memcpy(p_str_firm_data[i_n_firms].v_ret,&vv_tmp.vector);

            mv_tmp = gsl_matrix_submatrix(m_tmp,0,1,i,i_ncol - 3);
            gsl_matrix_scale(&mv_tmp.matrix, 1200); // factors scaled to annual returns
            gsl_matrix_memcpy(p_str_firm_data[i_n_firms].m_factors,&mv_tmp.matrix);

            // compute XTy and XTX
            v_y = p_str_firm_data[i_n_firms].v_ret;
            m_X = p_str_firm_data[i_n_firms].m_factors;
          //  v_XTy = p_str_firm_data[i_n_firms].v_XTy;
          //  m_XTX = p_str_firm_data[i_n_firms].m_XTX;
          //  olsg(v_y, m_X, v_XTy, m_XTX);


            i_firmid = i_id;
            i=0;
            i_n_firms++;
        }

        tok = strtok(NULL,","); // date

        for(j=0;j<i_ncol - 2;j++){
            tok = strtok(NULL,",");
            d_tmp = atof(tok);
            gsl_matrix_set(m_tmp,i,j,d_tmp);
        }

        i++;
    }
    fclose(ptr_file);

    // fill in struct of last firm
    p_str_firm_data[i_n_firms].i_id = i_firmid;
    p_str_firm_data[i_n_firms].i_ni = i;
    p_str_firm_data[i_n_firms].i_si = 0;
    p_str_firm_data[i_n_firms].v_ret = gsl_vector_alloc(i);
    p_str_firm_data[i_n_firms].m_factors = gsl_matrix_alloc(i,i_ncol - 3);
    p_str_firm_data[i_n_firms].v_XTy = gsl_vector_alloc(i_ncol-3);
    p_str_firm_data[i_n_firms].m_XTX = gsl_matrix_alloc(i_ncol-3,i_ncol-3);
    p_str_firm_data[i_n_firms].v_alpha = gsl_vector_alloc(i);

    p_str_firm_data[i_n_firms].d_rho = 0.0;
    p_str_firm_data[i_n_firms].d_mu = 0.0;
    p_str_firm_data[i_n_firms].d_w2 = 1.0;


    vv_tmp = gsl_matrix_subcolumn(m_tmp,0,0,i);
    gsl_vector_memcpy(p_str_firm_data[i_n_firms].v_ret,&vv_tmp.vector);

    mv_tmp = gsl_matrix_submatrix(m_tmp,0,1,i,i_ncol - 3);
    gsl_matrix_memcpy(p_str_firm_data[i_n_firms].m_factors,&mv_tmp.matrix);

    // compute XTy and XTX
    v_y = p_str_firm_data[i_n_firms].v_ret;
    m_X = p_str_firm_data[i_n_firms].m_factors;
    //v_XTy = p_str_firm_data[i_n_firms].v_XTy;
    //m_XTX = p_str_firm_data[i_n_firms].m_XTX;
    //olsg(v_y, m_X, v_XTy, m_XTX);

    gsl_matrix_free(m_tmp);
    return p_str_firm_data;
}

void alloc_mcmc_matrix(struct str_firm_data* a_firm_data, int i_n_firms, int i_draws, int i_factors)
{
    int i;

    for(i=0; i<i_n_firms; i++){
        a_firm_data[i].d_alpha = 0;
        a_firm_data[i].d_s2 = 1.0;
        a_firm_data[i].v_beta = gsl_vector_alloc(i_factors);
        gsl_vector_set_all(a_firm_data[i].v_beta,0.0);
        a_firm_data[i].m_mcmc_draws = gsl_matrix_alloc(i_draws,i_factors+2); // row store of draws [alpha, beta', s2]
     }

}

int write_mf_draws(struct str_firm_data *a_firm_data, int i_n_firms)
{
    int i;
    int i_id;
    int i_status;

    FILE *outputfile = fopen("test","w");

    char s_dirname[8] = "results";
    char *s_dirfirm;
    //asprintf(&s_dirname,"%s","results/");

    i_status = mkdir (s_dirname, S_IRWXU | S_IRWXG );

    for(i=0;i<i_n_firms;i++){
        i_id = a_firm_data[i].i_id;

        asprintf(&s_dirfirm, "%s/%d",s_dirname,i_id);

        i_status = mkdir( s_dirfirm, S_IRWXU | S_IRWXG);


         //write_posterior_summary(a_firm_data[i].m_mcmc_draws, s_dirname, outputfile,5);

        write_posterior_summary( a_firm_data[i].m_mcmc_draws, s_dirfirm, -1);

    }

    return 0;
}

int write_mfs_alpha_draws(struct str_firm_data *a_firm_data, int i_n_firms)
{
    int i,j;
    size_t i_draws = a_firm_data[0].m_mcmc_draws->size1;

    FILE *outputfile = fopen("mfs_alpha_draws.txt","w");

    for(i=0;i<i_draws;i++){
        for (j=0;j<i_n_firms-1;j++){

            fprintf(outputfile,"%12.8g,", mget(a_firm_data[j].m_mcmc_draws,i,0) );

        }

        fprintf(outputfile," %12.8g\n", mget(a_firm_data[i_n_firms-1].m_mcmc_draws,i,0) );
    }


  return 0;
}

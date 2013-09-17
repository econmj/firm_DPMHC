#include "distr.h"



/* Normal pdf */
double den_nor ( double x, double m, double s2 )
{
  double inv_std,e,p;

  inv_std=1.0/sqrt(s2);
  e=(x-m)*inv_std;
  p=gsl_ran_ugaussian_pdf ( e )*inv_std;
  return p;
}




/* Normal pdf based on precision */
double den_norm_prec ( double x, double m, double prec )
{
  double inv_std,e,p;

  inv_std=sqrt(prec);
  e=(x-m)*inv_std;
  p=gsl_ran_ugaussian_pdf ( e )*inv_std;
  return p;
}



/* return the log value of the Normal density at x, with mean mu, and var s2 */
double log_nor(double x, double mu, double s2)
{
  double u = x-mu;
  double p = -0.5*log(2*M_PI*s2) - 0.5*u*u / s2;
  return p;
}


/* return t-density evaluated at x, with mean m, inverse scale inv_s2, and degree of
   freedom df */
double den_st_prec ( double x, double m, double inv_s2, double df)
{
  double inv_s,e,p;
  
  inv_s=sqrt(inv_s2);
  e=( x - m )*inv_s;
  p=gsl_ran_tdist_pdf (e , df)*inv_s;
  return p;
}


/* return t-density with mean mu, and var s2*df/(df-2) if nu>2 */
double den_st(double x, double mu, double s2, double df)
{
  double inv_stdev = 1.0/sqrt(s2);
  double y = (x-mu)*inv_stdev;
  return gsl_ran_tdist_pdf(y,df) * inv_stdev;
}


/* computes log density of multivariate normal */
/*  */
double log_mnv_den(const gsl_vector *x, const gsl_vector *M, const gsl_matrix *V)
{
  size_t i,n=x->size;
  double qf,log_det_L;
  gsl_matrix *chol = gsl_matrix_alloc(n,n);
  gsl_vector *y = gsl_vector_alloc(n);

  gsl_matrix_memcpy(chol,V);
  gsl_linalg_cholesky_decomp(chol);

  /* compute determinant of L */
  log_det_L=0.0;
  for(i=0;i<n;i++)
    log_det_L += log(mget(chol,i,i));

  gsl_vector_memcpy(y,x);
  gsl_vector_sub (y, M);	/* y = (x-M) */
  cholesky_partial_solve(chol,y); /* y = chol^{-1}*y */
  gsl_blas_ddot (y,y,&qf);

  gsl_matrix_free(chol);
  gsl_vector_free(y);

  return -0.5*( n*log(2.0*M_PI) + qf) - log_det_L;
}


/* computes log density of multivariate t density St(x|M,V,df)*/
/* mean = M, Var(x) = V*df/(df-2) */
double log_mvt_den(const gsl_vector *x, const gsl_vector *M, double df, const gsl_matrix *V)
{
  size_t i,j,n=x->size;
  double qf,log_det_L,c_term,b_term;
  gsl_matrix *chol = gsl_matrix_alloc(n,n);
  gsl_vector *y = gsl_vector_alloc(n);
 
  gsl_matrix_memcpy(chol,V);
  gsl_linalg_cholesky_decomp(chol);
 
  /* compute determinant of L */
  log_det_L=0.0;
  for(i=0;i<n;i++)
    log_det_L += log(mget(chol,i,i));
 
  gsl_vector_memcpy(y,x);
  gsl_vector_sub (y, M);        /* y = (x-M) */
  cholesky_partial_solve(chol,y); /* y = chol^{-1}*y */
  gsl_blas_ddot (y,y,&qf);
  
  c_term = gsl_sf_lngamma( 0.5*(df+n) ) - gsl_sf_lngamma( 0.5*df ) - n*0.5*log( M_PI*df );
  b_term = -0.5*(df+n)*log(1.0 + qf/df);
 
 
  gsl_matrix_free(chol);
  gsl_vector_free(y);

  return c_term + b_term - log_det_L;
}


/* multivariate normal density function    */
/*
*       n       dimension of the random vetor
*       mean    vector of means of size n
*       var     variance matrix of dimension n x n
*/
double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mn, const gsl_matrix *vr)
{
int s;
double ax,ay;
gsl_vector *ym, *xm;
gsl_matrix *work = gsl_matrix_alloc(n,n),
           *winv = gsl_matrix_alloc(n,n);
gsl_permutation *p = gsl_permutation_alloc(n);

gsl_matrix_memcpy( work, vr );
gsl_linalg_LU_decomp( work, p, &s );
gsl_linalg_LU_invert( work, p, winv );
ax = gsl_linalg_LU_det( work, s );
gsl_matrix_free( work );
gsl_permutation_free( p );

xm = gsl_vector_alloc(n);
gsl_vector_memcpy( xm, x);
gsl_vector_sub( xm, mn );
ym = gsl_vector_alloc(n);
gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
gsl_matrix_free( winv );
gsl_blas_ddot( xm, ym, &ay);
gsl_vector_free(xm);
gsl_vector_free(ym);
ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );

return ay;
}

/****************************************************************************
Generalized Student t density found in Zellner (1998) Appendix B.5, p.396 
with P = 1 (scalar, so p = 1) and T is a random 1 x q vector
****************************************************************************/
double gen_Student_t(gsl_vector *v_t, gsl_matrix *m_Q, int i_n)
{
  size_t i,q=v_t->size;
  double d_kinv,d_num,d_den;
  double d_det_Q,d_f;

  if (i_n <= q){
    printf("Error in gen_Student_t( ): i_n must be greater than length of v_t\n");
    exit(1);
  }
  
  if ((q != m_Q->size1 ) || (q != m_Q->size2)){
    printf("Error in gen_Student_t(): the row/colums of m_Q must equal length of v_t\n");
    exit(1);
  }

  gsl_matrix *m_Chol = gsl_matrix_alloc(q,q);
  gsl_matrix *m_tt = gsl_matrix_alloc(q,q);


  d_kinv = pow(M_PI,(double)q/2.0);

  d_num = d_den = 1.0;
  for(i=1;i<=q;i++){
    d_num *= gsl_sf_gamma( 0.5*(i_n - 1 - i + 1) );
    d_den *= gsl_sf_gamma( 0.5*(i_n - i + 1) );
  }
  d_kinv *= d_num/d_den;

      
  /* determinant of m_S0 */
  gsl_matrix_memcpy(m_Chol,m_Q);
  gsl_linalg_cholesky_decomp(m_Chol);
  d_det_Q = 1.0;
  for(i=0;i<q;i++){
    d_det_Q *= mget(m_Chol,i,i);
  }
  d_det_Q *= d_det_Q;

  d_f = pow(d_det_Q,0.5*(i_n - 1));

  gsl_matrix_memcpy(m_tt,m_Q);
  gsl_blas_dger(1.0,v_t,v_t,m_tt); // m_tt = Q + t*t'
  
  gsl_matrix_memcpy(m_Chol,m_tt);
  gsl_linalg_cholesky_decomp(m_Chol);
  d_det_Q = 1.0;
  for(i=0;i<q;i++){
    d_det_Q *= mget(m_Chol,i,i);  // nb d_det_Q = |Q + T*T'|
  }
  d_det_Q *= d_det_Q;

  d_f *= pow( d_det_Q, -0.5*i_n );

  d_f = d_f/d_kinv;

  return d_f;
}


/* Inverse-Gamma(v,s) pdf */
/* x^{-nu-1} exp(-s/x)    */
double dens_inv_gamma(const double d_x, const double d_nu, const double d_s)
{
    double d_f;
    
    d_f = pow(d_s,d_nu)/gsl_sf_gamma(d_nu);
    d_f *= pow(d_x,-d_nu-1.0);
    d_f *= exp(-d_s/d_x);

    return d_f;
}

/* Log Transform of the Inverse-Gamma(v,s) pdf */
/* x^{-nu-1} exp(-s/x)    */
double ln_dens_inv_gamma(const double d_x, const double d_nu, const double d_s)
{
    double d_lnf;
    
    d_lnf = d_nu*log(d_s) - gsl_sf_lngamma(d_nu);
    d_lnf += (-d_nu-1.)*log(d_x);
    d_lnf -= d_s/d_x;

    return d_lnf;
}

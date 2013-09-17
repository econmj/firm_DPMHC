#include "linear_model.h"


extern gsl_rng *rng;


/* Compute Gibbs sampling pieces X'*y, and X'*X */
void olsg(gsl_vector *y, gsl_matrix *X, gsl_vector *XTy, gsl_matrix *XTX)
{
  gsl_vector_set_all(XTy,0.0);
  gsl_blas_dgemv (CblasTrans,1.0,X,y, 0.0 , XTy);

  gsl_matrix_set_all(XTX,0.0);
  /* XTX stored in lower triangle only */
  gsl_blas_dsyrk (CblasLower, CblasTrans,1.0,X,0.0,XTX);
}



/* beta | s2 ~ N(M,V^{-1})  */
/* V = tau*XTX + T0, M = V^{-1}*(tau*XTy + T0*b0) */
void linear_gibbs_beta(const gsl_matrix *XTX, const gsl_vector *XTy,
		       const gsl_vector *b0, const double s2,
		       const gsl_matrix *T0, gsl_vector *draw)
{
  size_t k=XTy->size;
  double tau=1./s2;
  gsl_matrix *V=gsl_matrix_alloc(k,k);
  gsl_matrix *VI=gsl_matrix_alloc(k,k);
  gsl_vector *tmp=gsl_vector_alloc(k);
  gsl_vector *M=gsl_vector_alloc(k);

  /* compute V = tau*XTX + T0 */
  gsl_matrix_memcpy(V,XTX);
  gsl_matrix_scale(V,tau);
  gsl_matrix_add(V,T0);

  /* compute V inverse = VI */
  gsl_linalg_cholesky_decomp (V);
  gsl_matrix_set_identity(VI);
  gsl_blas_dtrsm (CblasLeft, CblasLower,CblasNoTrans,CblasNonUnit,1.0,V,VI);
  gsl_blas_dtrsm (CblasLeft, CblasLower,CblasTrans,CblasNonUnit,1.0,V,VI);


  /* form T0*b0 + tau*XTy */
  gsl_vector_memcpy(tmp,XTy);
  gsl_blas_dsymv (CblasLower,1.0,T0,b0,tau,tmp);
  /* form V^{-1}*( T0*b0 + tau*XTy ) */
  gsl_vector_set_all(M,0.0);
  gsl_blas_dsymv (CblasLower,1.0,VI,tmp,0.0,M);

  ran_mvn(M,VI,draw);

  gsl_matrix_free(V);
  gsl_matrix_free(VI);
  gsl_vector_free(tmp);
  gsl_vector_free(M);
}


/* compute draw from sigma^2 | beta ~ IG( (T+v)/2 , (s0 +
   (y-X*beta)'(y-X*beta))/2 ) **/
double linear_gibbs_sigma2(const gsl_vector *y, const gsl_matrix *X,
		       const gsl_vector *beta, const double v0,
		       const double s0)
{
  size_t nobs=y->size;
  gsl_vector *u=gsl_vector_alloc(nobs);
  double ss;

  /* form u=y-X*beta */
  gsl_vector_memcpy(u,y);
  gsl_blas_dgemv (CblasNoTrans,-1.0,X,beta,1.0,u);
  gsl_blas_ddot (u,u,&ss);
  gsl_vector_free(u);
  return ran_invgamma( rng, .5*(nobs+v0), .5*(s0+ss) );
}




/* Gibbs sampling of beta | sigma2 for linear model with X  n X 1 */
/* prior is N(m,s2/v2) */
/* s2 is regression model innovation variance */
double sample_gibbs_mean(gsl_vector *y, gsl_vector *x,
			 double m, double v2, double s2)
{
  size_t i,n=y->size;
  double xt,yt,x2,xy,var_beta,mean_beta;

/*   x2=0.0; */
/*   xy=0.0; */
/*   for(i=0;i<n;i++){ */
/*     xt=vget(x,i); */
/*     yt=vget(y,i); */
/*     x2 += xt*xt; */
/*     xy += xt*yt; */
/*   } */
  gsl_blas_ddot(y,x,&xy);
  gsl_blas_ddot(x,x,&x2);


  var_beta = s2*v2 / (v2*x2 + s2);
  mean_beta = var_beta * ( xy/s2 + m/v2 );

  return mean_beta + gsl_ran_gaussian_ziggurat(rng,sqrt(var_beta));
}


/* Gibbs sampling of sigma2 | beta for linear model with X  n X 1 */
/* prior is IG(v0/2,s0/2) */
/* Conditional Posterior is IG((v0+nobs)/2,(ss+s0)/2) */
/* ss= sum (y-x*beta)^2 */
double sample_gibbs_var(gsl_vector *y, gsl_vector *x,
			 double beta, double v0, double s0)
{
  size_t i,nobs=y->size;
  double ss,u;

  ss = 0.0;
  for(i=0;i<nobs;i++){
    u = ( vget(y,i) - vget(x,i)*beta );
    ss += u*u;
  }

  return ran_invgamma( rng, .5*(nobs+v0), .5*(s0+ss) );
}

gsl_vector* linear_ols_beta(gsl_vector *v_y, gsl_matrix *m_X){
    size_t i_k = m_X->size2;

    gsl_vector *v_XTy = gsl_vector_alloc(i_k);
    gsl_vector *v_betahat = gsl_vector_alloc(i_k);
    gsl_matrix *m_XTX = gsl_matrix_alloc(i_k,i_k);
    gsl_matrix *m_invXTX = gsl_matrix_alloc(i_k,i_k);

    gsl_vector_set_all(v_XTy,0);
    gsl_vector_set_all(v_betahat,0);
    gsl_matrix_set_all(m_XTX,0);

    olsg(v_y,m_X,v_XTy,m_XTX);

    gsl_linalg_cholesky_decomp (m_XTX);
    gsl_matrix_set_identity(m_invXTX);
    gsl_blas_dtrsm (CblasLeft, CblasLower,CblasNoTrans,CblasNonUnit,1.0,m_XTX,m_invXTX);
    gsl_blas_dtrsm (CblasLeft, CblasLower,CblasTrans,CblasNonUnit,1.0,m_XTX,m_invXTX);

    gsl_vector_set_all(v_betahat,0.0);
    gsl_blas_dsymv (CblasLower,1.0,m_invXTX,v_XTy,0.0,v_betahat);

    gsl_vector_free(v_XTy);
    gsl_matrix_free(m_XTX);
    gsl_matrix_free(m_invXTX);

    return v_betahat;
}



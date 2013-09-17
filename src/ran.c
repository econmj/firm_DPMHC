/* Random number generators in addition to GSL
   Seed is assumed to be set in main and external */

#include "ran.h"

extern gsl_rng *rng;



/* GSL gamma uses scale, BS use inverse  scale */
/* This is a scale version of Gamma */
/* G(a,b) density with E(x) = a/b */
double ran_gamma(gsl_rng *r, double v, double s)
{
return gsl_ran_gamma_mt (r,v,1./s);
}

// NOTE IG(v,s) has E(x) = s/(v-1)
double ran_invgamma(gsl_rng *r,double v,double s)
{
  return 1./gsl_ran_gamma_mt (r,v,1./s);
}


/* Return a single draw from normal-gamma density
   It is h ~ G(v/2,s/2), b|h ~ N ( b_0, (h*tau)^{-1} )
   The results is placed in the gsl vector x[0]=h, x[1]=b */

void normal_gamma_univ(gsl_vector *x,double v, double s,double b_0, double tau)
{
  double h,beta,sigma;

  h=ran_gamma(rng, .5*v, .5*s);
  sigma=1./sqrt(h*tau);
  beta=b_0+gsl_ran_gaussian_ziggurat(rng,sigma);
  vset(x,0,beta);
  vset(x,1,h);
}

/* Return a single draw from normal-gamma density
   It is sigma^2 == 1/h ~ IG(v/2,s/2), b|sigma^2~ N ( b_0,  sigma^2/tau )
   The results is placed in the gsl vector x[0]=h, x[1]=b */

void normal_invgamma_univ(gsl_vector *x,double v, double s,double b_0, double tau)
{
  double h,beta,sigma;

  h=ran_gamma(rng, .5*v, .5*s);
  sigma=1./sqrt(h*tau);
  beta=b_0+gsl_ran_gaussian_ziggurat(rng,sigma);
  vset(x,0,beta);
  vset(x,1,1./h);
}

// A single draw from discrete density P(X=i)=p[i], i=0,...,n
// Returns i given gsl vector of probs. (NOTE SUPPORT IS n+1 DISCRETE POINTS
// NOT n)
size_t ran_multinomial(const gsl_vector *p,const size_t n)
{
  size_t i,vec_size;
  double sm,u;

  vec_size=p->size;
  if(vec_size < n){
    printf("\n **Error in ran_multinomial VECTOR DIM ERROR ** \n\n");
    exit(EXIT_FAILURE);
  }
  if(n < 1){
    printf("\n **Nothing to do in ran_multinomial** \n\n");
    return 0;
  }

  u=gsl_rng_uniform(rng);
  sm=0.;
  for(i=0;i<=n;i++){
    sm+=vget(p,i);
    if(u < sm){
      return i;
     }
  }
  printf("\n **Error in ran_multinomial** sm=%lf \n\n",sm);
  //  dgmpnt_gsl_vec(sm,stdout);
  exit(EXIT_FAILURE);
}


/* int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mn,  */
/* 	    const gsl_matrix *v, gsl_vector *result) */
/* { */
/* /\* multivariate normal distribution random number generator *\/ */
/* /\* */
/* *	n	dimension of the random vetor */
/* *	mean	vector of means of size n */
/* *	v	variance matrix of dimension n x n */
/* *	result	output variable with a sigle random vector normal distribution generation */
/* *\/ */
/*   int k; */
/*   gsl_matrix *work = gsl_matrix_alloc(n,n); */

/*   gsl_matrix_memcpy(work,v); */
/*   printf("copy v to work, before decomp. PRINT work \n");  */
/*   pmat(work); */


/*   gsl_linalg_cholesky_decomp(work); */
/*   //  printf("Cholesky decomposition \n");  */
/*   //  pmat(work); */
/*   //  printf("%f\n",mget(work,1,1)); */

/*   //  exit(0); */

/* for(k=0; k<n; k++) */
/* 	gsl_vector_set( result, k, gsl_ran_ugaussian(r) ); */

/* gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result); */
/* gsl_vector_add(result,mn); */

/* gsl_matrix_free(work); */

/* return 0; */
/* } */


/* Draw from the Multivariate Normal Distr */
/* mu=mean, V=Covariance matrix, draw=random draw */
int ran_mvn(const gsl_vector *mu, const gsl_matrix *V, gsl_vector *draw)
{
  size_t i,j,n=mu->size;
  double sm;
  gsl_matrix *Lower = gsl_matrix_alloc(n,n);

  gsl_matrix_memcpy(Lower,V);
  gsl_linalg_cholesky_decomp(Lower);

  for(i=0;i<n;i++)
    vset(draw,i,gsl_ran_ugaussian(rng));

  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, Lower, draw);
  gsl_vector_add(draw,mu);
  gsl_matrix_free(Lower);
  return 0;
}

/* Alternative (SLOWER)  Draw from the Multivariate Normal Distr */
int ran_mvn_alt(const gsl_vector *mu, const gsl_matrix *V, gsl_vector *draw)
{
  size_t n=mu->size;
  int i,j;
  double sm;
  gsl_matrix *Lower = gsl_matrix_alloc(n,n);

  gsl_matrix_memcpy(Lower,V);
  gsl_linalg_cholesky_decomp(Lower);

  for(i=0;i<n;i++)
    vset(draw,i,gsl_ran_ugaussian(rng));

  for(i=n-1;i >= 0;i--){
    sm=vget(mu,i);
    for(j=0;j<=i;j++){
      sm+=mget(Lower,i,j)*vget(draw,j);
    }
    vset(draw,i,sm);
  }

  gsl_matrix_free(Lower);

  return 0;
}


/* Draw from the Multivariate Student-t Distr */
/* mu=mean, V=Covariance matrix, df=degree of freedom, draw=random draw */
/* Generated as t = mu + CHOL(V)*(1.0/v)*z* */
/* z ~ MVN(mu,V), and v ~ G(df/2,df/2) */
int ran_mvt(const gsl_vector *mu, const gsl_matrix *V, double df, gsl_vector *draw)
{
  size_t i,j,n=mu->size;
  double sm,nu,v,inv_v;
  gsl_matrix *Lower = gsl_matrix_alloc(n,n);

  nu = df*0.5;
  //v = gsl_ran_gamma (rng ,nu , 1.0/nu );
  v = ran_gamma(rng, nu , nu);
  inv_v=sqrt( 1.0 / v );
  gsl_matrix_memcpy(Lower,V);
  gsl_linalg_cholesky_decomp(Lower);

  for(i=0;i<n;i++)
    vset(draw,i,gsl_ran_ugaussian(rng)*inv_v );

  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, Lower, draw);
  gsl_vector_add(draw,mu);
  gsl_matrix_free(Lower);
  return 0;
}

int ran_wish(const gsl_matrix *m_S0, const int i_v, gsl_matrix *m_Draw)
{
  size_t i;
  size_t p = m_S0->size1;

  if (p != m_S0->size2){
    printf("Error in ran_wish - m_S0 not symmetrical\n");
    exit(1);
  }

  if ((p != m_Draw->size1) || (p != m_Draw->size2)){
    printf("Error in ran_wish - m_Draw must have the same dim as m_S0");
    exit(1);
  }

  if (i_v < p){
    printf("Error in ran_wish - i_v must be greater than or equal to dim(m_S0)\n");
    exit(1);
  }

  gsl_matrix *m_y = gsl_matrix_alloc(p,(size_t) i_v);
  gsl_vector *v_zero = gsl_vector_alloc(p);
  gsl_vector *v_e = gsl_vector_alloc(p);
  gsl_vector_set_zero(v_zero);

  for(i=0;i<i_v;i++){
    ran_mvn(v_zero,m_S0,v_e);                         //v_e ~ N(0,m_S0)
    gsl_matrix_set_col(m_y,i,v_e);
  }

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, m_y, m_y, 0.0, m_Draw);

  gsl_vector_free(v_zero);
  gsl_vector_free(v_e);
  gsl_matrix_free(m_y);
  return 0;
}


/********************************************************************************
ran_inv_wish - generates a draw from the Inverse-Wishart(m_S0,i_v) distribution
where m_S0 is a p X p positive definite symmetrical matrix and i_v >= p.
********************************************************************************/
int ran_inv_wish(const gsl_matrix *m_S0, const int i_v, gsl_matrix *m_Draw)
{
  size_t p = m_S0->size1;

  if (p != m_S0->size2){
    printf("Error in ran_wish - m_S0 not symmetrical\n");
    exit(1);
  }

  if ((p != m_Draw->size1) || (p != m_Draw->size2)){
    printf("Error in ran_wish - m_Draw must have the same dim as m_S0");
    exit(1);
  }

  if (i_v < p){
    printf("Error in ran_wish - i_v must be greater than or equal to dim(m_S0)\n");
    exit(1);
  }

  gsl_matrix *m_S0_inv = gsl_matrix_alloc(p,p);
  gsl_matrix_memcpy(m_S0_inv, m_S0);
  inv_pd_sym( m_S0_inv );

  ran_wish(m_S0_inv,i_v,m_Draw);
  inv_pd_sym( m_Draw );

  return 0;
}


double ran_trunc_normal(const double d_mu, const double d_s2,
			const double d_a, double d_b)
{
  if (d_a >= d_b){
    printf("Error in ran_trunc_normal(): d_a and d_b must not be equal and d_a must be less than d_b\n");
    exit(1);
  }


  double d_ax = (d_a - d_mu)/sqrt(d_s2);
  double d_bx = (d_b - d_mu)/sqrt(d_s2);

  double d_draw,d_xstar,d_u;

  d_xstar = GSL_MIN(fabs(d_ax),fabs(d_bx));


  // Case (a,infty)
  if ( (gsl_isinf(d_bx) == 1) && (d_ax >= 0.0) ){

    if (d_ax <= 0.45){  // normal rejection

      d_draw = gsl_ran_gaussian_ziggurat(rng,1.0);
      while (d_draw < d_ax){
	d_draw = gsl_ran_gaussian_ziggurat(rng,1.0);
      }

    }
    else {  // exponential rejection

	d_draw = d_ax + gsl_ran_exponential(rng, 1./d_ax);
	d_u = gsl_ran_flat(rng,0,1);
	while ( d_u > exp( d_ax*d_draw - 0.5*(d_draw*d_draw + d_ax*d_ax) ) )  {
	  /* nb the proability of acceptance of Geweke (1991) hsould have lambda = -a
	     this has been correct for in the above while ( ) condition */
	  d_draw = d_ax + gsl_ran_exponential(rng, 1./d_ax);
	  d_u = gsl_ran_flat(rng,0,1);
	}

    }
    return d_mu + sqrt(d_s2)*d_draw;
  }


  // Case (-infty,b) b <= 0
  if ( (gsl_isinf(d_ax) == -1) && (d_bx <= 0.0) ){

    if (d_bx >= -0.45){  // normal rejection
      d_draw = gsl_ran_gaussian_ziggurat(rng,1.0);
      while (d_draw > d_bx){
	d_draw = gsl_ran_gaussian_ziggurat(rng,1.0);
      }

    }
    else {  // exponential rejection

	d_draw = d_bx - gsl_ran_exponential(rng, -1./d_bx);
	d_u = gsl_ran_flat(rng,0,1);
	while ( d_u > exp( d_bx*d_draw - 0.5*(d_draw*d_draw + d_bx*d_bx) ) ) {
	  /* nb the proability of acceptance of Geweke (1991) hsould have lambda = -a
	     this has been correct for in the above while ( ) condition */
	  d_draw = d_bx - gsl_ran_exponential(rng, -1./d_bx);
	  d_u = gsl_ran_flat(rng,0,1);
	}

    }
    return d_mu + sqrt(d_s2)*d_draw;
  }

  // Case (a,b) where 0 elem (a,b)
  if ( (d_ax < 0) && (d_bx > 0) ){

    if ( (d_ax < -0.375) ||
	 ( d_bx  > 0.375) ){

      d_draw = gsl_ran_gaussian_ziggurat(rng,1.0);
      while ( (d_draw < d_ax) || (d_draw > d_bx) ){
	d_draw = gsl_ran_gaussian_ziggurat(rng,1.0);
      }

    }
    else {

      d_draw = gsl_ran_flat(rng,d_ax,d_bx);
      d_u = gsl_ran_flat(rng,0,1);
      while (d_u > gsl_ran_gaussian_pdf(d_draw,1.0)/gsl_ran_gaussian_pdf(d_xstar,1.0)){
	d_draw = gsl_ran_flat(rng,d_ax,d_bx);
	d_u = gsl_ran_flat(rng,0,1);
      }

    }
    return d_mu + sqrt(d_s2)*d_draw;
  }

  // Case (a,b) where 0 <= a < b < infty
  if ( (gsl_isinf(d_bx) != 1) && (d_bx > 0) && (d_ax >= 0) ){

    if (gsl_ran_gaussian_pdf(d_ax,1.0)/gsl_ran_gaussian_pdf(d_bx,1.0) <= 2.18){

      d_draw = gsl_ran_flat(rng,d_ax,d_bx);
      d_u = gsl_ran_flat(rng,0,1);
      while (d_u > gsl_ran_gaussian_pdf(d_draw,1.0)/gsl_ran_gaussian_pdf(d_xstar,1.0)){
	d_draw = gsl_ran_flat(rng,d_ax,d_bx);
	d_u = gsl_ran_flat(rng,0,1);
      }

    }
    else {
      if (d_ax <= 0.725){

	d_draw = gsl_ran_gaussian_ziggurat(rng,1.0);
	while ( (d_draw < d_ax) || (d_draw > d_bx) ){
	  d_draw = gsl_ran_gaussian_ziggurat(rng,1.0);
	}
	d_draw = fabs(d_draw);

      }
      else {

	d_draw = d_ax + gsl_ran_exponential(rng, 1./d_ax);
	d_u = gsl_ran_flat(rng,0,1);
	while ( (d_draw > d_bx) ||
		( d_u > exp( d_ax*d_draw - 0.5*(d_draw*d_draw + d_ax*d_ax) ) ) ) {
	  /* nb the proability of acceptance of Geweke (1991) hsould have lambda = -a
	     this has been correct for in the above while ( ) condition */
	  d_draw = d_ax + gsl_ran_exponential(rng, 1./d_ax);
	  d_u = gsl_ran_flat(rng,0,1);
	}

      }
    }
    return d_mu + sqrt(d_s2)*d_draw;
  }

  // Case (a,b) where -infty < a < b <= 0
  if ( (gsl_isinf(d_ax) != -1) && (d_ax < 0) && (d_bx <= 0) ){

    if (gsl_ran_gaussian_pdf(d_bx,1.0)/gsl_ran_gaussian_pdf(d_ax,1.0) <= 2.18){

      d_draw = gsl_ran_flat(rng,d_ax,d_bx);
      d_u = gsl_ran_flat(rng,0,1);
      while (d_u > gsl_ran_gaussian_pdf(d_draw,1.0)/gsl_ran_gaussian_pdf(d_xstar,1.0)){
	d_draw = gsl_ran_flat(rng,d_ax,d_bx);
	d_u = gsl_ran_flat(rng,0,1);
      }

    }
    else {
      if (d_bx >= -0.725){

	d_draw = gsl_ran_gaussian_ziggurat(rng,1.0);
	while ( (fabs(d_draw) < -d_bx) || (fabs(d_draw) > -d_ax) ){
	  d_draw = gsl_ran_gaussian_ziggurat(rng,1.0);
	}

      }
      else {

	d_draw = d_bx - gsl_ran_exponential(rng, -1./d_bx);
	d_u = gsl_ran_flat(rng,0,1);
	while ( (d_draw < d_ax) ||
		( d_u > exp( d_bx*d_draw - 0.5*(d_draw*d_draw + d_bx*d_bx) ) ) ) {
	  /* nb the proability of acceptance of Geweke (1991) hsould have lambda = -a
	     this has been correct for in the above while ( ) condition */
	  d_draw = d_bx - gsl_ran_exponential(rng, -1./d_bx);
	  d_u = gsl_ran_flat(rng,0,1);
	}
      }
    }

    return d_mu + sqrt(d_s2)*d_draw;
  }

}


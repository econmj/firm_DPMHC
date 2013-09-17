#include "matrix.h"

#include "jensen_util.h"
// NOTE ran_gamma_mt (r,alpha,theta) has E(x)=alpha*theta


/* copy x[s:s+n-1] to y[t:t+n-1] */
void vec_copy(gsl_vector *x , int s, int n, gsl_vector *y, int t)
{
  int i,j,p;
  for(j=0;j<n;j++){
    i=s+j;
    p=t+j;
    vset(y,p, vget(x,i) );
  }
}


double sum(const gsl_vector *y)
{
  size_t i,n;
  double mn;

  n=y->size;
  mn=0.;
  for(i=0;i<n;i++){
	mn+=vget(y,i);
  }
  return mn;
}

double sumr(const gsl_vector *y,const size_t st,const size_t ed)
{
  size_t i;
  double mn;

  mn=0.;
  for(i=st;i<ed+1;i++){
	mn+=vget(y,i);
  }
  return mn;
}

double mean(const gsl_vector *y)
{
  size_t i,n;
  double mn;

  n=y->size;
  mn=0.;
  for(i=0;i<n;i++){
    mn+=vget(y,i);
  }
  return mn/((double)n);
}

/* sample mean of each column of data in matrix m */
/* returns a vector */
void mean_matrix_col(const gsl_matrix *A , gsl_vector *sample_mean)
{
  size_t i;
  size_t n=A->size1;

  gsl_vector_set_zero (sample_mean);
  for(i=0;i<n;i++){
    gsl_vector_const_view row = gsl_matrix_const_row (A, i);
    gsl_vector_add (sample_mean, &row.vector);
  }
  gsl_vector_scale (sample_mean,1./((double)n));
}

double var(const gsl_vector *y)
{
  size_t i,n;
  double mn,mn2;

  n=y->size;
  mn=0.;
  mn2=0.;
  for(i=0;i<n;i++){
    mn+=vget(y,i)/((double)n);
    mn2+=pow(vget(y,i),2.)/((double)n);
  }
  return mn2-mn*mn;
}

void acf(gsl_vector *y, long n, gsl_vector *rho)
{
  size_t nobs=y->size;
  size_t i,k;

  if(n > rho->size){
    printf("\n **acf()**  rho too small \n");
    exit(EXIT_FAILURE);
  }

  double m=mean(y);
  double inv_v=1.0/(var(y)*nobs);

  for(k=1;k<=n;k++){
    double tmp=0.0;
    for(i=k;i<nobs;i++)
      tmp += (vget(y,i)-m)*(vget(y,i-k)-m)*inv_v;
    vset(rho,k-1,tmp);
  }
}
  
  



/* computes sample Covariance matrix for column data */
void cov_matrix_col(const gsl_matrix *A, gsl_matrix *cov)
{
  size_t i,j,s;
  size_t n=A->size1;
  size_t k=A->size2;
  double tmp,mean_i,mean_j;
  gsl_vector *sample_mean=gsl_vector_alloc(k);

  mean_matrix_col(A ,sample_mean);

  for(i=0;i<k;i++){
    for(j=i;j<k;j++){
      tmp=0.;
      mean_i=vget(sample_mean,i);
      mean_j=vget(sample_mean,j);
      for(s=0;s<n;s++){
	tmp+=(mget(A,s,i)-mean_i)*(mget(A,s,j)-mean_j) /((double)n);
      }
      mset(cov,i,j,tmp);
      mset(cov,j,i,tmp);
    }
  }

  gsl_vector_free(sample_mean);
  
}


//Sort and return a and b (a<b) quantiles called q_a, and q_b for data vector y
void den_interval(const gsl_vector *y,const double a,double * q_a, 
		  const double b,double * q_b)
{
  gsl_vector *ysorted;
  size_t n,i,lower,upper;

  n=y->size;
  ysorted=gsl_vector_alloc(n);
  gsl_vector_memcpy(ysorted,y);
  gsl_sort_vector (ysorted);
  lower=(size_t)n*a-1;
  upper=(size_t)n*b-1;
  if(lower >n-1 || upper>n-1)
    {
      printf("INVALID ARGUEMENTS a, and b in den_interval\n");
      abort();
    }
  *q_a=vget(ysorted,lower);
  *q_b=vget(ysorted,upper);
  gsl_vector_free (ysorted);

}


//Computes X'*X and places results in XTX (upper left hand corner if XTX is a
//larger matrix
void mat_xtx(const gsl_matrix *x, gsl_matrix *xtx)
{
  size_t n_x,k_x,n_xtx,k_xtx,i,j,k;
  double sm;

  n_x=x->size1;
  k_x=x->size2;

  n_xtx=xtx->size1;
  k_xtx=xtx->size2;
  if( (n_xtx < k_x) || (k_xtx < k_x))
    {
      printf("\n\n **DIMENSIONS INCORRECT IN mat_xtx()** \n\n");
      abort();
    }


  for(j=0;j<k_x;j++)
    for(i=j;i<k_x;i++)
      {
	sm=0.;
	for(k=0;k<n_x;k++)
	  {
	    sm+=mget(x,k,i)*mget(x,k,j);
	  }
	mset(xtx,i,j,sm);
	mset(xtx,j,i,sm);
      }


}

/* print GSL matrix  */
void pmat(const gsl_matrix *x)
{
  size_t n,m,i,j;
  n=x->size1;
  m=x->size2;  
  printf("\n");
  for(i=0;i<n;i++)
    {
      for(j=0;j<m;j++)
	{
	  printf("%14.7e ",mget(x,i,j));
	}
      printf("\n");
    }
  printf("\n");
}

/* print GSL matrix with range X[a1:a2][b1:b2] */
void pmatr(const gsl_matrix *x, size_t a1,size_t a2,size_t b1, size_t b2)
{
  size_t n,m,i,j;
  printf("\n");
  for(i=a1;i<=a2;i++)
    {
      for(j=b1;j<=b2;j++)
	{
	  printf("%14.7e ",mget(x,i,j));
	}
      printf("\n");
    }
  printf("\n");
}


/* print GSL vector  */
void pvec(const gsl_vector *x)
{
  size_t n,i;
  n=x->size;
  printf("\n");
  for(i=0;i<n;i++)
    {
      printf("%14.7e\n",vget(x,i));
    }
  printf("\n");
}

/* print GSL vector_int  */
void pvec_int(const gsl_vector_int *x)
{
  size_t n,i;
  n=x->size;
  printf("\n");
  for(i=0;i<n;i++)
    {
      printf("%d\n",vget_int(x,i));
    }
  printf("\n");
}



/* print GSL vector with range st:ed */
void pvecr(const gsl_vector *x,size_t st,size_t ed)
{
  size_t i;

  printf("\n");
  for(i=st;i<ed+1;i++){
    printf("%14.7e\n",vget(x,i));
  }
  printf("\n");
}


void mypause()
{
     char dummy[10];
     //     printf("Press any key to continue\n");
     fgets(dummy , 10 ,stdin);
}

/* return nonparametric univariate density estimate at x given data vector*/
double dens_est(double x, gsl_vector *data, double h)
{
  double den;
  size_t n,i;

  n=data->size;

  den=0.;
  for(i=0;i<n;i++){
    den += epan( (x - vget(data,i))/h )/(n*h);
  }

  return den;
}


/* Epanechnikov Kernel */
double epan( double u)
{
  if( fabs(u) < 1.){
    return 0.75*(1.-pow(u,2));
  }else{
    return 0.;
  }
}

/* inverse of a positive definite symmtric matrix */
/* returned in A */
void inv_pd_sym( gsl_matrix *A )
{
  size_t n=A->size1;
  gsl_matrix *Lower=gsl_matrix_alloc(n,n);
  gsl_matrix_memcpy(Lower,A);
  gsl_linalg_cholesky_decomp(Lower);
  gsl_matrix_set_identity(A);
  gsl_blas_dtrsm (CblasLeft, CblasLower,CblasNoTrans,CblasNonUnit,1.0,Lower,A);
  gsl_blas_dtrsm (CblasLeft, CblasLower,CblasTrans,CblasNonUnit,1.0,Lower,A);
  gsl_matrix_free(Lower);
}

/* Solve for x in A*x = y, where A is the cholesky decomposition and a lower
   triangular matrix */
/* answer returned in y */
/* To COMPUTE QUADRATIC FORM y'*inv(V)*y */
/* 1. gsl_linalg_cholesky_decomp(chol); */
/* 2. cholesky_partial_solve(gsl_matrix *chol, gsl_matrix *y)   */
/* 3. Dot product <y,y> */
void  cholesky_partial_solve(const gsl_matrix *chol, gsl_vector  *y)  
{
  size_t n=chol->size1;
  int i,j;
  double sm;

  for(i=0;i<n;i++){
      sm=vget(y,i);
      for(j=0;j<i;j++){
	sm += -mget(chol,i,j)*vget(y,j);
      }
      sm = sm/mget(chol,i,i);
      vset(y,i,sm);
  }
}


gsl_matrix* matrix_drop_row(gsl_matrix *m_A, size_t j)
{
  size_t i,ii=0;
  size_t i_n = m_A->size1;
  if (j >= i_n){
    printf("Error in matrix_drop_row(): row argument, j, must be less than m_A->size1\n");
    exit(1);
  }

  gsl_vector *v_tmpi = gsl_vector_alloc(m_A->size2);
  gsl_matrix *m_Aj = gsl_matrix_alloc(i_n-1,m_A->size2);


  for(i=0;i<i_n;i++){

    if (i != j){
      gsl_matrix_get_row(v_tmpi,m_A,i);
      gsl_matrix_set_row(m_Aj,ii,v_tmpi);
      ii++;
    }
  }

  gsl_vector_free(v_tmpi);
  return m_Aj;
}


gsl_matrix* matrix_drop_col(gsl_matrix *m_A, size_t j)
{
  size_t i,ii=0;
  size_t i_n = m_A->size2;
  if (j >= i_n){
    printf("Error in matrix_drop_col(): col argument, j, must be less than m_A->size2\n");
    exit(1);
  }

  gsl_vector *v_tmpi = gsl_vector_alloc(m_A->size1);
  gsl_matrix *m_Aj = gsl_matrix_alloc(m_A->size1,i_n-1);


  for(i=0;i<i_n;i++){

    if (i != j){
      gsl_matrix_get_col(v_tmpi,m_A,i);
      gsl_matrix_set_col(m_Aj,ii,v_tmpi);
      ii++;
    }
  }

  gsl_vector_free(v_tmpi);
  return m_Aj;
}


gsl_matrix* matrix_drop_row_col(gsl_matrix *m_A, size_t j1, size_t j2)
{
  
  gsl_matrix *m_Ajj = gsl_matrix_alloc(m_A->size1 - 1, m_A->size2 - 1);
  gsl_matrix *m_tmpAj;

  m_tmpAj = matrix_drop_col(m_A,j2);
  m_Ajj = matrix_drop_row(m_tmpAj,j1);

  gsl_matrix_free(m_tmpAj);
  return m_Ajj;
}

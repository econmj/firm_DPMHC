#include "jensen_util.h"


void gsl_gnuplot_vector(gsl_vector *v_X)
{
  FILE *ftmp = fopen("ftmp.dat","w");
  gsl_vector_fprintf(ftmp,v_X,"%g");
  fclose(ftmp);

  FILE *gplot = popen("gnuplot -persist", "w");
  fprintf(gplot,"plot \"%s\" u 1 w l\n","ftmp.dat"); //,"ftmp.dat");
  fflush(gplot);
  pclose(gplot);

  FILE *tmp = popen("/bin/rm ftmp.dat","w");
  fflush(tmp);
  pclose(gplot);
}


/********************************************************/
/* TITLE: DGMPNT_GSL_VEC                                */
/*                                                      */
/* PURPOSE:                                             */
/* PRINT A GSL VECTOR                                   */
/*                                                      */
/* USAGE: dgmpnt_gsl_vec(v_A,outputfile);               */
/*                                                      */
/* ARGUMENTS                                            */
/* v_A  : AN GSL VECTOR (TYPE gsl_vector*)              */
/* outputfile : FILE TO WHICH THE MATRIX IS PRINTED     */
/*              (TYPE *FILE)                            */
/********************************************************/

void dgmpnt_gsl_vec(gsl_vector *v_A,FILE *outputfile)
{
  double f;
  int m=1,n,i,j,start,stop;
  int lnsize=80;                     /* lnsize can be changed up to 133 */
  int maxcol,ipad,ipad10,ipad1,lnsiz;
  char *fmt;

  n = v_A->size;
  lnsiz = lnsize;
  if((lnsize < 72) || (lnsize > 133)) lnsiz=133;
  maxcol = (lnsiz-8)/12;
  if(n < maxcol)  maxcol=n;
  ipad = (lnsiz-8-12*maxcol)/2+1;
  ipad10 = ipad/10;
  ipad1 = ipad-10*ipad10;

  start=1;
 L1: stop = start-1+maxcol;
  if(stop > n) stop=n;
  

  fprintf(outputfile,"\n\n");
  fprintf(outputfile,"      ");
  for(j=start; j<=stop; j++)
    {
      fprintf(outputfile,"      COL%3d",j);
    }
  fprintf(outputfile,"\n");

  for(i=1; i<=m; i++)
    {
      fprintf(outputfile,"ROW%3d",i);
      for(j=start; j<=stop; j++)
	{
	    f = fabs( gsl_vector_get(v_A,j-1) );
	    if( f < 1.e+8 ) fmt = "%12.0f";
	    if( f < 1.e+5 ) fmt = "%12.1f";
	    if( f < 1.e+4 ) fmt = "%12.2f";
	    if( f < 1.e+3 ) fmt = "%12.3f";
	    if( f < 1.e+2 ) fmt = "%12.4f";
	    if( f < 1.e+1 ) fmt = "%12.5f";
	    if( f < 1.e+0 ) fmt = "%12.6f";
	    if( f < 1.e-1 ) fmt = "%12.7f";
	    if( f < 1.e-2 ) fmt = "%12.8f";
	    if( f < 1.e-4 ) fmt = "%12.4e";
	    if( f < 1.e-55) fmt = "%12.1f";
	    fprintf(outputfile,fmt, gsl_vector_get(v_A,j-1));
	}
      fprintf(outputfile,"\n");
    }

  if(stop == n) return;
  start = stop+1;
  goto L1;
}


/********************************************************/
/* TITLE: DGMPNT_GSL_MATRIX                             */
/*                                                      */
/* PURPOSE:                                             */
/* PRINT A gls_matrix                                   */
/*                                                      */
/* USAGE: dgmpnt_gsl_matrix(m_A,outputfile);            */
/*                                                      */
/* ARGUMENTS                                            */
/* m_A  : AN m BY n MATRIX STORED COLUMNWISE            */
/*      (TYPE gsl_matrix*)                              */
/* outputfile : FILE TO WHICH THE MATRIX IS PRINTED     */
/*              (TYPE *FILE)                            */
/********************************************************/
                                                                             
void dgmpnt_gsl_matrix(gsl_matrix *m_A,FILE *outputfile)
{   
    double f;   
    int m,n,i,j,start,stop;   
    int lnsize=80;                     /* lnsize can be changed up to 133 */
    int maxcol,ipad,ipad10,ipad1,lnsiz;
    char *fmt;
                                                                             
    m = m_A->size1;
    n = m_A->size2;
    lnsiz = lnsize;
    if((lnsize < 72) || (lnsize > 133)) lnsiz=133;
    maxcol = (lnsiz-8)/12;
    if(n < maxcol)  maxcol=n;
    ipad = (lnsiz-8-12*maxcol)/2+1;
    ipad10 = ipad/10;
    ipad1 = ipad-10*ipad10;
                                                                             
    start=1;
    L1: stop = start-1+maxcol;
    if(stop > n) stop=n;
                                                                             
                                                                             
    fprintf(outputfile,"\n\n");
    fprintf(outputfile,"      ");
    for(j=start; j<=stop; j++)
      {
	  fprintf(outputfile,"      COL%3d",j);
      }
    fprintf(outputfile,"\n");
                                                                             
    for(i=1; i<=m; i++)
      {
	  fprintf(outputfile,"ROW%3d",i);
	  for(j=start; j<=stop; j++)
            {
		f = fabs( gsl_matrix_get(m_A,i-1,j-1) );
		if( f < 1.e+8 ) fmt = "%12.0f";
		if( f < 1.e+5 ) fmt = "%12.1f";
		if( f < 1.e+4 ) fmt = "%12.2f";
		if( f < 1.e+3 ) fmt = "%12.3f";
		if( f < 1.e+2 ) fmt = "%12.4f";
		if( f < 1.e+1 ) fmt = "%12.5f";
		if( f < 1.e+0 ) fmt = "%12.6f";
		if( f < 1.e-1 ) fmt = "%12.7f";
		if( f < 1.e-2 ) fmt = "%12.8f";
		if( f < 1.e-4 ) fmt = "%12.4e";
		if( f < 1.e-55) fmt = "%12.1f";
		fprintf(outputfile,fmt,gsl_matrix_get(m_A,i-1,j-1));
	    }
	  fprintf(outputfile,"\n");
      }
 
    if(stop == n) return;
    start = stop+1;
    goto L1;
}





/***************************************************************************
This code was generously provided to me by John Lamb from GSL
  discussion forum Estimate the Hessian of a function f. The vector x
  and matrix H should satisfy * x->size == H->size1 == H->size2 and h
  should be positive. The function checks these values and returns
  GSL_EINVAL or GSL_EDOM if they are not satisfied.
 
  The Hessian is estimated with central difference quotient estimates and uses
  \f$4n^2-2n+4\f$ function calls.
  
  @param x A vector.
  @param f A function that takes as its first argument a vector of the same 
  size as x. The void* is used to pass any additional function parameters.
  @param params Additional parameters to be given to each function call.
  @param h A small positive constant used in the calculations. See standard 
  texts on numerical analysis for typical values.
  @param H Stores the resulting Hessian.
  ************************************************************************* */

int hessian( const gsl_vector* x, double (*f)( const gsl_vector*, void* ),
	     void* params, double h, gsl_matrix* H )
{
  /* variables */
  size_t n;              /* size of x (and H) */
  size_t row, col;       /* for indexing H    */
  double pp, pm, mp, mm; /* temporaries for calculations */
  double cc, twh2, fh2, th;
  double xrow, xcol;     /* local storage */
  gsl_vector_view view;  /* last row of H */
  
  /* carry out some sanity checks */
  if( x == 0 || f == 0 || H == 0 )
    return GSL_EINVAL;
  n = x->size;
  if( H->size1 != n || H->size2 != n )
    return GSL_EINVAL;
  if( h <= 0 )
    return GSL_EDOM;
  
  /* set up */
  view = gsl_matrix_row( H, n - 1 );
  gsl_vector_memcpy( &view.vector, x );  /* use a copy of x */
  th = 2 * h;
  fh2 = th * th;
  twh2 = 3 * fh2;
  cc = f( &view.vector, params );
  
  /* calculate upper triangle */
  for( row = 0; row < n; ++row ){
    for( col = row; col < n; ++col ){
      /* copy element to be modified */
      xrow = gsl_vector_get( &view.vector, row );
      if( row == col ){
	gsl_vector_set( &view.vector, row, xrow + th );
	pp = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow + h );
	pm = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow - th );
	mm = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow - h );
	mp = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow );
	/* set matrix entry */
	gsl_matrix_set( H, row, col, (16 * (pm + mp) - (pp + 30 * cc + mm)) / twh2 );
	/* very last call will overwrite view and view should not then be changed */
      } else {
	/* copy element to be modified */
	xcol = gsl_vector_get( &view.vector, col );
	gsl_vector_set( &view.vector, row, xrow + h );
	gsl_vector_set( &view.vector, col, xcol + h );
	pp = f( &view.vector, params );
	gsl_vector_set( &view.vector, col, xcol - h );
	pm = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow - h );
	mm = f( &view.vector, params );
	gsl_vector_set( &view.vector, col, xcol + h );
	mp = f( &view.vector, params );
	gsl_vector_set( &view.vector, row, xrow );
	gsl_vector_set( &view.vector, col, xcol );
	/* set matrix entry */
	gsl_matrix_set( H, row, col, (pp + mm - (pm + mp)) / fh2 );
      }
    }
  }

  /* calculate lower triangle */
  for( row = 1; row < n; ++row ){ /* will work even for 1 x 1 matrix */
    for( col = 0; col < row; ++col ){
      gsl_matrix_set( H, row, col, gsl_matrix_get( H, col, row ) );
    }
  }

  return 0;
}





void Quasi_Newton_Hessian(gsl_matrix *m_QuasiHess, 
			  gsl_vector *v_grad_new, gsl_vector *v_x_new,
			  gsl_vector *v_grad_old, gsl_vector *v_x_old)
{
    /***********************************************************************/
    /* See p. 55 of Bonnas, J.F., Gilbert, J.C., Lemarechal, C. and        */
    /*   Sagastizabal, C.A., "Numerical Optimization: Theoretical and      */
    /*  Practical Aspects" BFGS formula for Quasi-Newton method:           */
    /*  Wnew = W - (sy'W+Wys')/(y*s) + [1+(y*(Wy))/(y*s)]ss'/(y*s)         */
    /*	where s = v_x_new - v_x_old; y = v_grad_new - v_grad_old           */
    /***********************************************************************/

  double d_ys, d_yWy, d_coef1;
  size_t dim = v_x_old->size;
  
  gsl_vector *v_y = gsl_vector_alloc(dim);
  gsl_vector *v_s = gsl_vector_alloc(dim);
  gsl_matrix *m_sy = gsl_matrix_alloc(dim, dim);
  gsl_matrix *m_syW = gsl_matrix_alloc(dim, dim);
  gsl_vector *v_Wy = gsl_vector_alloc(dim);
  gsl_matrix *m_ss = gsl_matrix_alloc(dim, dim);
  gsl_matrix_set_zero (m_sy);	gsl_matrix_set_zero (m_ss);
	

  gsl_vector_memcpy(v_y,v_grad_new); gsl_vector_memcpy(v_s,v_x_new);
  gsl_vector_sub (v_y, v_grad_old); /*y = v_grad_new - v_grad_old */
  gsl_vector_sub (v_s, v_x_old); /*s = v_x_new - v_x_old*/
  
  gsl_blas_ddot (v_y, v_s, &d_ys);/*(y*s) */
  if (d_ys>1.0e-16) /*This check is necessary to ensure numerical stability*/
    {
      gsl_blas_dger (1.0, v_s, v_y, m_sy);/*sy'*/
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, m_sy, m_QuasiHess, 0.0, m_syW);/*sy'W*/
      gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, m_QuasiHess, m_sy, 1.0, m_syW);/*sy'W+Wys'*/
      gsl_matrix_scale (m_syW, -1.0/d_ys);/*-(sy'W+Wys')/(y*s)*/
      
      gsl_blas_dgemv (CblasNoTrans, 1.0, m_QuasiHess, v_y, 0.0, v_Wy);/*Wy*/
      gsl_blas_ddot (v_y, v_Wy, &d_yWy); /*y*(Wy)*/
      d_coef1 = (1.0+d_yWy/d_ys)/d_ys;/*[1+(y*(Wy))/(y*s)]/(y*s)*/
      gsl_blas_dger (1.0, v_s, v_s, m_ss);/*ss'*/
      gsl_matrix_scale (m_ss, d_coef1);/*[1+(y*(Wy))/(y*s)]ss'/(y*s)*/
      
      gsl_matrix_add (m_syW, m_ss);/*-(sy'W+Wys')/(y*s)+[1+(y*(Wy))/(y*s)]ss'/(y*s)*/
      gsl_matrix_add (m_QuasiHess, m_syW);/*Wnew*/
    }
  else
    {
      /*              printf("\n (y*s) = %.8g \n",d_ys); */
    }

  gsl_vector_free (v_y);
  gsl_vector_free (v_s);
  gsl_vector_free (v_Wy);
  gsl_matrix_free (m_sy);
  gsl_matrix_free (m_syW);
  gsl_matrix_free (m_ss);
}


/***************************************************************************/
/* Purpose: Makes a random draw from a multivariate Student-t distribution */
/* t_{nu}(v_mu,m_Sigma) where nu is the degrees of freedom.  The random    */                
/* vector x ~ t_{nu}(v_mu,m_Sigma) has the properties, E(x) = v_mu,        */
/* and Cov(x) = (nu/(nu-2)) m_Sigma (nu > 2)                               */
/* See Gelman, Carlin, Stern, and Rubin (1995) Bayesian Data Analysis,     */
/* Chapman and Hall, p. 476, 481.                                          */
/***************************************************************************/
gsl_vector *jensen_ran_multv_tdist(const gsl_rng *gsl_r, gsl_vector *v_mu,
				   gsl_matrix *m_Sigma, double d_nu)
{
    int i,j,i_k;
    double d_x,d_mu;
    gsl_vector *v_z,*v_x;
    gsl_matrix *m_Chol;

    i_k = v_mu->size;
    if ((i_k != m_Sigma->size1) || (i_k != m_Sigma->size2))
       {
	   fprintf(stderr,"\n\nError in jensen_ran_multv_tdist(). Mean vector and Covariance matrix input must be of the same dimension.\n");
	   //error(1);
	
       }

    v_z = gsl_vector_alloc(i_k);
    v_x = gsl_vector_alloc(i_k);
    m_Chol = gsl_matrix_alloc(i_k,i_k);
    gsl_matrix_memcpy(m_Chol,m_Sigma);

    /* d_x ~ X^2(d_nu) */
    d_x = gsl_ran_chisq( gsl_r, d_nu);

    /* choleski decomposition of m_Sigma */

    if ( gsl_linalg_cholesky_decomp(m_Chol) )
      {
	  fprintf(stderr,"\n\nError in jensen_ran_multv_tdist(). Covariance must be postive definite\n");

	  for (i=0;i<i_k;i++)
	    v_x->data[i] = 0.0/0.0;   // should return nan
	  return v_x;
	  //	  error(1);
	
      }
    
    /* v_z ~ N(0,I) */
    for (i=0;i<i_k;i++)
        {
	    v_z->data[i] =  gsl_ran_gaussian( gsl_r,1.0);
	    v_x->data[i] = v_mu->data[i];
	    for (j=0;j<=i;j++)
	        {
		    v_x->data[i] += gsl_matrix_get(m_Chol,i,j)*v_z->data[j]*
			            sqrt(d_nu/d_x);
		}
	}

    gsl_vector_free(v_z);
    gsl_matrix_free(m_Chol);

    return v_x;
}



/******************************************************************/
/* Date: November 23, 2005                                        */
/* Function: Student_t_kernel( )                                  */
/* Purpose: Calculates the kernel (pdf w/ no normalizing constant)*/
/* of multivariate Student-t distribution, with mean v_mu and     */
/* covariance matrix: 
/*                m_Sigma * d_nu/(d_nu-2)                         */
/* See Gelman et al, 1995 Bayesian Data Analysis
/*                                                                */
/*                                                                */
/* f(v_x|nu,v_mu,m_Sigma)                                         */
/* \propto [1 + (v_x-v_mu)'(Sigma*nu)^{-1}(v_x-v_mu)]^{-(N+nu)/2} */
/*                                                                */
/*                                                                */
/* OUTPUT: (double) pdf value                                     */
/*                                                                */
/* INPUTS:                                                        */
/*        v_x (gsl_vector*)    : random vector                    */
/*        v_mu (gsl_vector*)   : distribution mean                */
/*        m_Sigma (gsl_matrix*): matrix                           */ 
/*        m_nu (double)        : degrees of freedom               */
/******************************************************************/

double Student_t_kernel(gsl_vector *v_x, gsl_vector *v_mu, 
			gsl_matrix *m_Sigma, double d_nu)
{
    int i_N;
    double d_kern,d_f;
    gsl_vector *v_demeaned,*v_tmp;
    
    i_N = v_x->size;
    v_demeaned = gsl_vector_alloc(i_N);
    v_tmp = gsl_vector_alloc(i_N);

    gsl_vector_memcpy(v_demeaned,v_x);

    /* v_demeaned = v_x - v_mu */
    gsl_vector_sub (v_demeaned, v_mu);


    /* inv(m_Sigma) */
    gsl_matrix *m_Tmp = gsl_matrix_alloc(i_N,i_N);
    gsl_matrix *m_Siginv = gsl_matrix_alloc(i_N,i_N);
    gsl_permutation *p = gsl_permutation_alloc(i_N);
    int s;
    
    gsl_matrix_memcpy(m_Tmp,m_Sigma);
    gsl_linalg_LU_decomp(m_Tmp,p,&s);
    gsl_linalg_LU_invert(m_Tmp,p,m_Siginv);

    /*  m_Siginv * (v_x - v_mu) */
    gsl_blas_dgemv(CblasTrans, 1.0, m_Siginv, v_demeaned, 0.0, v_tmp);
    
    /* (v_x - v_mu)' * v_temp1 */
    gsl_blas_ddot(v_demeaned, v_tmp, &d_kern);
    
    d_kern /= d_nu;


    gsl_vector_free(v_demeaned);
    gsl_vector_free(v_tmp);
    gsl_matrix_free(m_Tmp);
    gsl_matrix_free(m_Siginv);
    gsl_permutation_free(p);

    d_f =  pow(1.0 + d_kern, -0.5*(d_nu + (double) i_N));

    return d_f;
}



double find_machine_precision(void)
{
  double d_numb,d_val;
  
  d_val = 1.0;
  d_numb = 1.0;

  while ((d_val + d_numb * d_numb) != d_val)
    d_numb /= 10.0;

  return d_numb;
}


void pnt_gsl_matrix(gsl_matrix *m_A,FILE *outputfile)
{   
    double f;   
    int m,n,i,j;
    char *fmt;
                                                                             
    m = m_A->size1;
    n = m_A->size2;
                                                                             
                                                                             
    for(i=1; i<=m; i++)
      {
	  for(j=1; j<=n; j++)
            {
		f = fabs( gsl_matrix_get(m_A,i-1,j-1) );
		if( f < 1.e+8 ) fmt = "%12.0f";
		if( f < 1.e+5 ) fmt = "%12.1f";
		if( f < 1.e+4 ) fmt = "%12.2f";
		if( f < 1.e+3 ) fmt = "%12.3f";
		if( f < 1.e+2 ) fmt = "%12.4f";
		if( f < 1.e+1 ) fmt = "%12.5f";
		if( f < 1.e+0 ) fmt = "%12.6f";
		if( f < 1.e-1 ) fmt = "%12.7f";
		if( f < 1.e-2 ) fmt = "%12.8f";
		if( f < 1.e-4 ) fmt = "%12.4e";
		if( f < 1.e-55) fmt = "%12.1f";
		fprintf(outputfile,fmt,gsl_matrix_get(m_A,i-1,j-1));
	    }
	  fprintf(outputfile,"\n");
      }
 
    return;
}

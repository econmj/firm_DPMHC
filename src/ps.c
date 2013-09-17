#include "ps.h"


extern gsl_rng *rng;


void posterior_summary(const gsl_matrix *theta, FILE *ofile, long M)
{

  size_t T=theta->size1;
  size_t npar=theta->size2;
  gsl_vector *tmp=gsl_vector_alloc(T);
  int i,j;
  double median,lower,upper;

  printf("\n Writing MCMC draws to out\n\n");
  FILE *file = fopen("out","w");
  for(i=0;i<T;i++){
    for(j=0;j<npar;j++)
      fprintf(file,"%14.6e ",mget(theta,i,j));
    fprintf(file,"\n");
  }

  fprintf(ofile,"\n\n Posterior Summary \n");
  fprintf(ofile,"\n T=%lu\n\n",T);
  fprintf(ofile,"\n      Mean          Median         Stdev           0.95 DI\n\n");
  for(i=0;i<npar;i++){
    gsl_matrix_get_col( tmp, theta, i);
    gsl_sort_vector(tmp);
    median=gsl_stats_median_from_sorted_data(tmp->data,tmp->stride,tmp->size);
    lower=gsl_stats_quantile_from_sorted_data(tmp->data,tmp->stride,tmp->size,0.025);
    upper=gsl_stats_quantile_from_sorted_data(tmp->data,tmp->stride,tmp->size,0.975);

    fprintf(ofile,"%2d %14.7e %14.7e %14.7e (%14.7e,%14.7e)\n"
	   ,i,mean(tmp),median,sqrt(var(tmp)),lower,upper);
  }

  long tau;
  if( M < 0 )
    tau=1000;
  else
    tau=M;

  gsl_vector *rho=gsl_vector_alloc(tau);

  fprintf(ofile,"\n                                ACF");
  fprintf(ofile,"\n      NSE          Ineff        1            50           100          500\n");
  for(i=0;i<npar;i++){
    gsl_matrix_get_col( tmp, theta, i);
    acf(tmp,tau,rho);

    /* write out ACF for each parameter */
    char file_name[20] = "acf.dat";
    sprintf( &file_name[7],"%d",i);
    FILE *fp_acf = fopen( file_name, "w");

    for(j=0;j<tau;j++)
      fprintf(fp_acf,"%g\n",vget(rho,j));

    fclose(fp_acf);


    /* get inefficiency factor using Newey-West estimate of Long-run var*/
    double ineff=1.0;
    for(j=0;j<tau-1;j++){
      ineff += 2.0*(tau-j-1)/tau*vget(rho,j);
      }
    /* numerical standard error for posterior mean */
    double nse=sqrt(var(tmp)*ineff/T);

    fprintf(ofile,"%2d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n"
	    ,i,nse,ineff,vget(rho,0),vget(rho,49),vget(rho,99),vget(rho,499));

    /* produce kernel density plot for each parameter */
    char file_name2[20] = "den.dat";

    sprintf( &file_name2[7],"%d",i);
    FILE *fp_den = fopen( file_name2, "w");

    double stdev = sqrt(var(tmp));
    lower = gsl_vector_min(tmp) - stdev;
    upper = gsl_vector_max(tmp) + stdev;


    den_est_file(tmp, lower , upper ,100, fp_den, -1.0);


  }

  fprintf(ofile,"\n\n");
  gsl_vector_free(rho);
  gsl_vector_free(tmp);
}

/* kernel density estimate over range [Lower,Upper) for Points equals spaced.
   h=bandwidth, h<=0 gives default value . Written to file */
void den_est_file (gsl_vector *y , double Lower , double Upper, int Points , FILE *fp, double h)
{
  int i,j;
  size_t n = y->size;
  double d,f,x,x_j,mn,s2,stdev;

  d = (Upper - Lower)/(double)Points;

  mn = mean(y);
  s2 = var(y);
  stdev = sqrt(s2);


  if( h <= 0.0 )
    h = 0.9 *stdev*pow( (double)n, -0.2);

  for(i=0;i<Points;i++){
    x = Lower + (double)i * d;
    f = dens_est( x, y, h);
    fprintf(fp,"%g %g %g\n",x,f,den_nor(x,mn,s2));
  }
}


void write_posterior_summary(const gsl_matrix *theta, char *s_dirname, long M)
{

  size_t T=theta->size1;
  size_t npar=theta->size2;
  gsl_vector *tmp=gsl_vector_alloc(T);
  int i,j;
  double median,lower,upper;

  char *s_filename;


  //printf("\n Writing MCMC draws to out\n\n");
  asprintf(&s_filename,"%s/out",s_dirname);

  FILE *file = fopen(s_filename,"w");
  for(i=0;i<T;i++){
    for(j=0;j<npar;j++)
      fprintf(file,"%14.6e ",mget(theta,i,j));
    fprintf(file,"\n");
  }
  fclose(file);

  asprintf(&s_filename,"%s/mcmc_summary.txt",s_dirname);
  FILE *ofile = fopen(s_filename,"w");

  fprintf(ofile,"\n\n Posterior Summary \n");
  fprintf(ofile,"\n T=%lu\n\n",T);
  fprintf(ofile,"\n      Mean          Median         Stdev           0.95 DI\n\n");
  for(i=0;i<npar;i++){
    gsl_matrix_get_col( tmp, theta, i);
    gsl_sort_vector(tmp);
    median=gsl_stats_median_from_sorted_data(tmp->data,tmp->stride,tmp->size);
    lower=gsl_stats_quantile_from_sorted_data(tmp->data,tmp->stride,tmp->size,0.025);
    upper=gsl_stats_quantile_from_sorted_data(tmp->data,tmp->stride,tmp->size,0.975);

    fprintf(ofile,"%2d %14.7e %14.7e %14.7e (%14.7e,%14.7e)\n"
	   ,i,mean(tmp),median,sqrt(var(tmp)),lower,upper);
  }

  long tau;
  if( M < 0 )
    tau=1000;
  else
    tau=M;

  gsl_vector *rho=gsl_vector_alloc(tau);

  fprintf(ofile,"\n                                ACF");
  fprintf(ofile,"\n      NSE          Ineff        1            50           100          500\n");
  for(i=0;i<npar;i++){
    gsl_matrix_get_col( tmp, theta, i);
    acf(tmp,tau,rho);

    /* write out ACF for each parameter */
    asprintf(&s_filename,"%s/acf.dat%d",s_dirname,i);
    FILE *fp_acf = fopen( s_filename, "w");

    for(j=0;j<tau;j++)
      fprintf(fp_acf,"%g\n",vget(rho,j));

    fclose(fp_acf);


    /* get inefficiency factor using Newey-West estimate of Long-run var*/
    double ineff=1.0;
    for(j=0;j<tau-1;j++){
      ineff += 2.0*(tau-j-1)/tau*vget(rho,j);
      }
    /* numerical standard error for posterior mean */
    double nse=sqrt(var(tmp)*ineff/T);

    fprintf(ofile,"%2d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n"
	    ,i,nse,ineff,vget(rho,0),vget(rho,49),vget(rho,99),vget(rho,499));

    /* produce kernel density plot for each parameter */
    asprintf(&s_filename,"%s/den.dat%d",s_dirname,i);
    FILE *fp_den = fopen( s_filename, "w");

    double stdev = sqrt(var(tmp));
    lower = gsl_vector_min(tmp) - stdev;
    upper = gsl_vector_max(tmp) + stdev;


    den_est_file(tmp, lower , upper ,100, fp_den, -1.0);
    fclose(fp_den);

  }

  fprintf(ofile,"\n\n");
  fclose(ofile);

  gsl_vector_free(rho);
  gsl_vector_free(tmp);
}

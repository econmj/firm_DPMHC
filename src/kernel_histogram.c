#include "kernel_histogram.h"


int smooth_density_vector(gsl_vector *v_x){

    char file_name[20] = "smooth_den.dat";
    FILE *fp_den = fopen( file_name, "w");

    double d_stdev = sqrt(var(v_x));
    double d_low = gsl_vector_min(v_x) - d_stdev;
    double d_high = gsl_vector_max(v_x) + d_stdev;

    den_est_file( v_x, d_low, d_high, 100, fp_den, -1);

    fclose(fp_den);

    return 0;
}

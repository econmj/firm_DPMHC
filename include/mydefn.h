// GSL shortforms

#define vget gsl_vector_get
#define vset gsl_vector_set 
#define vfree gsl_vector_free
#define vget_int gsl_vector_int_get
#define vset_int gsl_vector_int_set 
#define vfree_int gsl_vector_int_free

#define mset gsl_matrix_set
#define mget gsl_matrix_get
#define mfree gsl_matrix_free
#define mset_int gsl_matrix_int_set
#define mget_int gsl_matrix_int_get
#define mfree_int gsl_matrix_int_free

/* subvector view */
#define vv(a) &(a).vector
/* submatrix view */
#define mv(b) &(b).matrix

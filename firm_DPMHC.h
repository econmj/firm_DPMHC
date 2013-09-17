#ifndef FIRM_DPMHC_H_INCLUDED
#define FIRM_DPMHC_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include "DPMHC_smplr.h"
#include "read_mf_data.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


int firm_DPMHC(struct str_firm_data *a_firm_data, int i_J,
                    struct str_DPMHC *ptr_DPMHC_alpha);



#endif // FIRM_DPMHC_H_INCLUDED

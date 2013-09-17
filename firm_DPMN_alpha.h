#ifndef FIRM_DPMN_ALPHA_H_INCLUDED
#define FIRM_DPMN_ALPHA_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include "DPMN_polya_smplr.h"
#include "read_mf_data.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


int firm_DPMN_alpha(struct str_firm_data *a_firm_data, int i_J,
                    struct str_DPMN *ptr_DPMN_alpha);


#endif // FIRM_DPMN_ALPHA_H_INCLUDED

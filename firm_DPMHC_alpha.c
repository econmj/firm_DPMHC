#include "firm_DPMHC_alpha.h"

int firm_DPMHC_alpha(struct str_firm_data *a_firm_data, int i_J, struct str_DPMHC *ptr_DPMHC_alpha){

    int i;
    for(i=0;i<i_J;i++){
        vset(ptr_DPMHC_alpha->v_y, i, a_firm_data[i].d_alpha);
    }

    DPMHC_smplr(ptr_DPMHC_alpha, a_firm_data);

    return 0;
}





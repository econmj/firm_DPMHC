#include "firm_DPMHC.h"

int firm_DPMHC(struct str_firm_data *a_firm_data, int i_J, struct str_DPMHC *ptr_DPMHC_alpha){

    int i;
    for(i=0;i<i_J;i++){
        vset(ptr_DPMHC_alpha->v_y, i, a_firm_data[i].d_alpha);
    }

    DPMHC_smplr(ptr_DPMHC_alpha, i_J, a_firm_data);

    return 0;
}





#include "firm_DPMN_alpha.h"

int firm_DPMN_alpha(struct str_firm_data *a_firm_data, int i_J,
                    struct str_DPMN *ptr_DPMN_alpha)
{
    int i;
    for(i=0;i<i_J;i++){
        vset(ptr_DPMN_alpha->v_y, i, a_firm_data[i].d_alpha);
    }

/*
    FILE *testfile = fopen("alpha.out","w");
    for(i=0;i<i_J;i++){
        fprintf(testfile, "%8.12f\n",vget(ptr_DPMN_alpha->v_y,i));
    }
    fclose(testfile);
    FILE *gplot = popen("gnuplot -persist", "w");
    fprintf(gplot, "plot \"%s\" u 1 w l\n", "alpha.out");
    fflush(gplot);
    close(gplot);
    //exit(1);
*/

    DPMN_polya_smplr(ptr_DPMN_alpha);



    return 0;
}

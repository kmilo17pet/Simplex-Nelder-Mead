#include <stdio.h>
#include <stdlib.h>
#include "nmsimplex.h"

double rosen(double v[]){
    double x=v[1],y=v[0];
    double a=1.0, b=100.0;
    return pow(a-x , 2.0) + b*pow(y-pow(x,2.0),2.0);
}

void constraints(double *x, size_t n){
    
}

nmsimplex_data_t MinDat = {.ConstraintFcn=NULL ,.epsilon=DBL_EPSILON, .Iterations=10000};
int main(void) {
    double params[] = {2, 4};
    double min;
    min=nmsimplex(&MinDat, rosen, params, 2);
    printf("\r\nMinimal:%g\r\nIterations:%d\r\nF.Eval:%d ", min, MinDat.Iterations, MinDat.FcnEvaluations);
    printf("\r\nResult:%g %g\r\n",params[0],params[1]);
    return (EXIT_SUCCESS);
}
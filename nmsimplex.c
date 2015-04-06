/* 
 * File:   nmsimplex.h
 * Author: Eng. Juan Camilo Gómez Cadavid MSc.
 *
 * Created on 31 de diciembre de 2011, 08:52 AM
 */
/*
 * algorithm presented in Margaret H. Wright's paper on Direct Search Methods.
 */
#include "nmsimplex.h"

static double __nmsimplex_fcn_evaluate(nmsimplex_data_t *obj, double* data);

/*============================================================================*/
static double __nmsimplex_fcn_evaluate(nmsimplex_data_t *obj, double* data){
    if (obj->ConstraintFcn) obj->ConstraintFcn(data, obj->n);
    obj->FcnEvaluations++;
    return obj->objfunc(data);
}
/*============================================================================*/
double nmsimplex(nmsimplex_data_t *obj, double (*objfunc)(double*), double* seed, size_t n){	
    int vs,vh,vg;   // vertexes: vs: smallest, vh: next smallest, vg: largest	
    int i,j,row;
    double **v;     /* holds vertices of simplex */
    double pn,qn;   /* values used to create initial simplex */
    double *f;      /* value of function at each vertex */
    double fr,fe,fc;      // function evaluation at: fr: reflection, fe: expansion, fc: contraction
    double *vr,*ve,*vc,*vm;  // coordinates: vr-reflection, ve-expansion, vc-contraction, vm-centroid
    double favg,s;
    obj->n  = n;
    obj->objfunc = objfunc;
    if(obj->alpha<=0.0) obj->alpha = 1.0;
    if(obj->gamma<1.0) obj->gamma = 2.0;
    if(obj->beta<=0.0 || obj->beta>1.0) obj->beta = 0.5;
    if(obj->delta<=0.0 || obj->delta>1.0) obj->delta = 0.5;
    if(obj->epsilon<=0.0 || obj->epsilon> 1E-8 )  obj->epsilon = 1E-12;
    if(obj->InitSimplexscale<=0.0)  obj->InitSimplexscale = 1.0;
    if(obj->MaxIterations<=0.0)  obj->MaxIterations = 1000;
    v =  (double**)calloc((n+1) , sizeof(double*));
    for (i=0; i<=n; i++) v[i] = (double*)calloc(n , sizeof(double)); // allocate the columns of the arrays
    f =  (double*)calloc((n+1) , sizeof(double));
    vr = (double*)calloc(n , sizeof(double));
    ve = (double*)calloc(n , sizeof(double));  
    vc = (double*)calloc(n , sizeof(double));  
    vm = (double*)calloc(n , sizeof(double));  
    obj->Iterations = obj->FcnEvaluations = 0;  
    #define __SQRT2   (1.4142135623730950488016887242097)
    pn = obj->InitSimplexscale*(sqrt(n+1)-1+n)/(n*__SQRT2); // create the initial simplex : assume one of the vertices is 0,0
    qn = obj->InitSimplexscale*(sqrt(n+1)-1)/(n*__SQRT2);
    for (i=1; i<=n; i++) {
        v[0][i-1] = seed[i-1];
        for (j=0;j<n;j++)   v[i][j] = ((i-1 == j)? pn : qn) + seed[j];
    }
    
    if (obj->ConstraintFcn)  obj->ConstraintFcn(v[j],n);
    for (j=0; j<=n; j++) f[j] = objfunc(v[j]); // find the initial function values 
    obj->FcnEvaluations = j;
   	
    for (obj->Iterations=1;obj->Iterations<=obj->MaxIterations;obj->Iterations++) { // begin the main loop of the minimization
        for (vh=1, vg=0, vs=0, j=0; j<=n; j++) {
            if(f[j] > f[vg]){ //find the index of the largest and 2nd largest values
                vh = vg;
                vg = j;
            }            
            else if(f[j]>f[vh] && f[i]!=f[vg]) vh = j; 
            if (f[j] < f[vs]) vs = j; // find the index of the smallest value      
        }
        
        for ( j=0; j<=n-1; j++) { // compute the centroid vm[j]
            vm[j]=0.0;
            for (i=0; i<=n; i++) vm[j] += v[i][j];    
            vm[j] = (vm[j]-v[vg][j])/n;  //largest excluded from centroid
        }
		
        for (j=0; j<=n-1; j++) vr[j] = vm[j] + obj->alpha*(vm[j]-v[vg][j]); // reflect vg to new vertex vr
        fr = __nmsimplex_fcn_evaluate(obj, vr);
        if (fr < f[vh] && fr >= f[vs]){
            for (j=0; j<=n-1; j++) v[vg][j] = vr[j];
            f[vg] = fr;
        }
        
        if ( fr < f[vs]){ // expansion
            for (j=0; j<=n-1; j++) ve[j] = vm[j] + obj->gamma*(vr[j]-vm[j]);
            fe = __nmsimplex_fcn_evaluate(obj, ve);
            if (fe < f[vs]/*fr*/){
                for (j=0; j<=n-1; j++) v[vg][j] = ve[j];
                f[vg] = fe;
            }
            else{
                for (j=0; j<=n-1; j++) v[vg][j] = vr[j];
                f[vg] = fr;
            }
        }
       
        if (fr >= f[vh]){ // check if contraction is necessary
            if (fr < f[vg] && fr >= f[vh]) { 
                for (j=0; j<=n-1; j++) vc[j] = vm[j]+obj->beta*(vr[j]-vm[j]); // outside contraction
            }
            else{ 
                for (j=0; j<=n-1; j++) vc[j] = vm[j]-obj->beta*(vm[j]-v[vg][j]); //inside contraction
            }
            fc = __nmsimplex_fcn_evaluate(obj, vc);
            if (fc <= f[vg]) {
                for (j=0;j<=n-1;j++) v[vg][j] = vc[j];
                f[vg] = fc;
            }
            else{ // contraction failed and shrinkage will be the next attemp
                for (row=0; row<=n; row++) {
                    if (row != vs) {
                        for (j=0; j<=n-1; j++) v[row][j] = v[vs][j]+ obj->delta*(v[row][j]-v[vs][j]); // shrinkage
                    }
                }
                f[vg] = __nmsimplex_fcn_evaluate(obj, v[vg]);
                f[vh] = __nmsimplex_fcn_evaluate(obj, v[vh]);
            }
        }
		
        for (favg=0.0, j=0; j<=n; j++) favg += f[j];
        favg /= (n+1);

        for (s=0.0, j=0; j<=n; j++) s += pow((f[j]-favg),2.0)/n;
        s = sqrt(s);
        if (s < obj->epsilon) break; 
    }
    	
    for (vs=0, j=0; j<=n; j++) {
        if (f[j] < f[vs]) vs = j; // find the index of the smallest value
    }
    for (j=0; j<n; j++) seed[j] = v[vs][j]; // put optimization result into "seed"
      
    obj->FcnEvaluations++;
    free(f);free(vr);free(ve);free(vc);free(vm);
    for (i=0; i<=n; i++) free(v[i]);
    free(v);
    return objfunc(v[vs]);
}
/*============================================================================*/
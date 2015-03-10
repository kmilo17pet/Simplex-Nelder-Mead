/* 
 * File:   nmsimplex.h
 * Author: juanc_000
 *
 * Created on 31 de diciembre de 2014, 08:52 AM
 */

#ifndef NMSIMPLEX_H
#define	NMSIMPLEX_H

#ifdef	__cplusplus
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

typedef struct{
    void (*ConstraintFcn)(double*, size_t); //ConstraintFunction
    unsigned long MaxIterations;
    double alpha; // Reflection 
    double beta; // Contraction 
    double gamma; // Expansion
    double delta; // Shrinkage
    double epsilon; 
    double InitSimplexscale;
    int FcnEvaluations;
    int Iterations;
    int n;
    double (*objfunc)(double*);
}nmsimplex_data_t;    

double nmsimplex(nmsimplex_data_t *obj, double (*objfunc)(double*), double *seed, size_t n);

#ifdef	__cplusplus
}
#endif

#endif	/* NMSIMPLEX_H */


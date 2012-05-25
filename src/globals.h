#ifndef _GLOBALS_H_
#define _GLOBALS_H_


#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef float single;

/* Define the simParam structure */
typedef struct {
    single dt;
    int Nt, Ny, Nx;
    single dy, dx;
    int Ns, Np, Npxy;
    single *P;
    single *Pxy;
    single *bcMtop, *bcMbot, *bcMrig, *bcMlef;
    int useRK4, useGPU, useGPUnum;
} simParam;


// prototype for GPU function
void differentiateStateVectorField_d( single *Sprime, single *S, simParam sp );

#endif // #ifndef  _GLOBALS_H_

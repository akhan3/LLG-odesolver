/* ==========================================================================
 * odeStepComp.c
 *
 * This is the top level for compiling odeStepComp mex-Funtion.
 *
 * Project:     Massive ODE Solver
 * Author:      Aamir Ahmed Khan (akhan3@nd.edu)
 * Copyright:   2012, University of Notre Dame
 * =========================================================================*/


#include "mex.h"
#include "string.h"
#include "assert.h"
#include "math.h"
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
    int useRK4, useGPU;
} simParam;


/* IMPORTANT: Don't mess with this function!!! Just use it */
/*! Picks the correct neighbors based on the location
 *      Picks the neighbor if exists otherwise pick the corresponding boundary condition
 *  \param S_top Return pointer for top neighbor vector
 *  \param S_bot Return pointer for bottom neighbor vector
 *  \param S_rig Return pointer for right neighbor vector
 *  \param S_lef Return pointer for left neighbor vector
 *  \param S State vector
 *  \param sp Simulation parameters
 *  \param ix index number on x-axis
 *  \param iy index number on y-axis */
void pickNeighbors( single **S_top, single **S_bot, single **S_rig, single **S_lef,
                    single *S, simParam sp, int ix, int iy )
{
    /* Here, top means maximum iy coordinate, not the first matrix row
     *       This code follows xy-coordinate axes convention */
    // if at top-most row, use top boundary condition
    *S_top = (iy == sp.Ny-1)  ?  &sp.bcMtop[ix*sp.Ns]  :  &S[(iy+1)*sp.Ns+ ix   *sp.Ny*sp.Ns]; // +y (wrt to xy-cood axes)
    // if at bottom-most row, use bottom boundary condition
    *S_bot = (iy == 0)        ?  &sp.bcMbot[ix*sp.Ns]  :  &S[(iy-1)*sp.Ns+ ix   *sp.Ny*sp.Ns]; // -y (wrt to xy-cood axes)
    // if at right-most column, use right boundary condition
    *S_rig = (ix == sp.Nx-1)  ?  &sp.bcMrig[iy*sp.Ns]  :  &S[ iy   *sp.Ns+(ix+1)*sp.Ny*sp.Ns]; // +x
    // if at left-most column, use left boundary condition
    *S_lef = (ix == 0)        ?  &sp.bcMlef[iy*sp.Ns]  :  &S[ iy   *sp.Ns+(ix-1)*sp.Ny*sp.Ns]; // -x
}


/* Including the user defined source file containing the model */
#include "modelingEquations.c"


/*! Takes one ODE step by Euler's method.
 *  \param Snext State at the next time instant
 *  \param sp Simulation parameters
 *  \param S Current state */
void eulerStep( single *Snext, single *S, simParam sp )
{
    int Nxy = sp.Ny * sp.Nx;    // size of array
    /* Allocate memory for Sprime vector field */
    single *Sprime = (single*)calloc( sp.Ns*Nxy, sizeof(single) );
    /* Calculate the slopes by executing differential equations */
    differentiateStateVectorField( Sprime, S, sp );
    /* Advance to the next time instant */
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int ixy = 0; ixy < Nxy; ++ixy ) {
        int i = ixy*sp.Ns;
        Snext[i+0] = S[i+0] + Sprime[i+0] * sp.dt;
        Snext[i+1] = S[i+1] + Sprime[i+1] * sp.dt;
        Snext[i+2] = S[i+2] + Sprime[i+2] * sp.dt;
        /* compensate the Snext vector after update */
        compensateStateVector(&Snext[i], &S[i], sp, ixy);
    }
    /* clean-up */
    // must deallocate all the reserved memory, otherwise huge performance penalty!!!
    free(Sprime);
}

/*! Takes one ODE step by Runge-Kutta's 4th order method.
 *  \param Snext State at the next time instant
 *  \param sp Simulation parameters
 *  \param S Current state */
void rk4Step( single *Snext, single *S, simParam sp )
{
    int Nxy = sp.Ny * sp.Nx;    // size of array
    /* Allocate memory for Sprime and cumulative-slope vector field */
    single *Sprime = (single*)calloc( sp.Ns*Nxy, sizeof(single) );
    single *slope = (single*)calloc( sp.Ns*Nxy, sizeof(single) );

    /* k1: Calculate the slopes by executing differential equations */
    {
        single *k1 = Sprime;    // k1 <= Sprime (just a reference. Memory is still allocated only once)
        differentiateStateVectorField( k1, S, sp );
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for( int ixy = 0; ixy < Nxy; ++ixy ) {
            int i = ixy*sp.Ns;
            /* Advance half the step by k1 */
            Snext[i+0] = S[i+0] + k1[i+0] * .5*sp.dt;
            Snext[i+1] = S[i+1] + k1[i+1] * .5*sp.dt;
            Snext[i+2] = S[i+2] + k1[i+2] * .5*sp.dt;
            /* Compensate the Snext vector after update */
            compensateStateVector(&Snext[i], &S[i], sp, ixy);
            /* Initialize the cumulative slope by k1/6 */
            slope[i+0] = k1[i+0] / 6.0;
            slope[i+1] = k1[i+1] / 6.0;
            slope[i+2] = k1[i+2] / 6.0;
        }
    }
    /* k2: Calculate the slopes by executing differential equations */
    {
        single *k2 = Sprime;    // k2 <= Sprime (just a reference. Memory is still allocated only once)
        differentiateStateVectorField( k2, Snext, sp );
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for( int ixy = 0; ixy < Nxy; ++ixy ) {
            int i = ixy*sp.Ns;
            /* Advance half the step by k2 */
            Snext[i+0] = S[i+0] + k2[i+0] * .5*sp.dt;
            Snext[i+1] = S[i+1] + k2[i+1] * .5*sp.dt;
            Snext[i+2] = S[i+2] + k2[i+2] * .5*sp.dt;
            /* Compensate the Snext vector after update */
            compensateStateVector(&Snext[i], &S[i], sp, ixy);
            /* Update the cumulative slope by adding k2/3 */
            slope[i+0] += k2[i+0] / 3.0;
            slope[i+1] += k2[i+1] / 3.0;
            slope[i+2] += k2[i+2] / 3.0;
        }
    }
    /* k3: Calculate the slopes by executing differential equations */
    {
        single *k3 = Sprime;    // k3 <= Sprime (just a reference. Memory is still allocated only once)
        differentiateStateVectorField( k3, Snext, sp );
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for( int ixy = 0; ixy < Nxy; ++ixy ) {
            int i = ixy*sp.Ns;
            /* Advance the full step by k3 */
            Snext[i+0] = S[i+0] + k3[i+0] * sp.dt;
            Snext[i+1] = S[i+1] + k3[i+1] * sp.dt;
            Snext[i+2] = S[i+2] + k3[i+2] * sp.dt;
            /* Compensate the Snext vector after update */
            compensateStateVector(&Snext[i], &S[i], sp, ixy);
            /* Update the cumulative slope by adding k3/3 */
            slope[i+0] += k3[i+0] / 3.0;
            slope[i+1] += k3[i+1] / 3.0;
            slope[i+2] += k3[i+2] / 3.0;
        }
    }
    /* k4: Calculate the slopes by executing differential equations */
    {
        single *k4 = Sprime;    // k4 <= Sprime (just a reference. Memory is still allocated only once)
        differentiateStateVectorField( k4, Snext, sp );
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for( int ixy = 0; ixy < Nxy; ++ixy ) {
            int i = ixy*sp.Ns;
            /* Complete the cumulative slope by finally adding k4/6 */
            slope[i+0] += k4[i+0] / 6.0;
            slope[i+1] += k4[i+1] / 6.0;
            slope[i+2] += k4[i+2] / 6.0;
            /* Advance the full step by cumulative skope for the next time instant */
            Snext[i+0] = S[i+0] + slope[i+0] * sp.dt;
            Snext[i+1] = S[i+1] + slope[i+1] * sp.dt;
            Snext[i+2] = S[i+2] + slope[i+2] * sp.dt;
            /* Compensate the Snext vector after final update */
            compensateStateVector(&Snext[i], &S[i], sp, ixy);
        }
    }

    /* clean-up */
    // must deallocate all the reserved memory, otherwise huge performance penalty!!!
    free(Sprime);
    free(slope);
}


/*! Converts Matlab's simParam structure to an equivalent C structure.
 *  \param sp_in Matlab's version of structure
 *  \return Equivalent C structure */
simParam parseSimParam( const mxArray *sp_in )
{
    /* get the pointer to the nested boundary-condition structure */
    mxArray *bc = mxGetField(sp_in, 0, "boundCond");
    single *bcMtop = (single*)mxGetData( mxGetField(bc, 0, "S_top") );

    /* fill up the sp structure */
    simParam sp = {
        .dt = mxGetScalar(mxGetField(sp_in, 0, "dt")),
        .Nt = mxGetScalar(mxGetField(sp_in, 0, "Nt")),
        .Ny = mxGetScalar(mxGetField(sp_in, 0, "Ny")),
        .Nx = mxGetScalar(mxGetField(sp_in, 0, "Nx")),
        .dy = mxGetScalar(mxGetField(sp_in, 0, "dy")),
        .dx = mxGetScalar(mxGetField(sp_in, 0, "dx")),
        .Ns = mxGetScalar(mxGetField(sp_in, 0, "Ns")),
        .Np = mxGetScalar(mxGetField(sp_in, 0, "Np")),
        .Npxy = mxGetScalar(mxGetField(sp_in, 0, "Npxy")),

        .P   = (single*)mxGetData(mxGetField(sp_in, 0, "P")),
        .Pxy = (single*)mxGetData(mxGetField(sp_in, 0, "Pxy")),

        .bcMtop = (single*)mxGetData( mxGetField(bc, 0, "S_top") ),
        .bcMbot = (single*)mxGetData( mxGetField(bc, 0, "S_bot") ),
        .bcMrig = (single*)mxGetData( mxGetField(bc, 0, "S_rig") ),
        .bcMlef = (single*)mxGetData( mxGetField(bc, 0, "S_lef") ),

        .useRK4 = mxGetScalar(mxGetField(sp_in, 0, "useRK4")),
        .useGPU = mxGetScalar(mxGetField(sp_in, 0, "useGPU")),
    };

    // const mwSize *dimsP = mxGetDimensions(mxGetField(sp_in, 0, "P"));
    // const mwSize *dimsPxy = mxGetDimensions(mxGetField(sp_in, 0, "Pxy"));
    // sp.Np = dimsP[1];       // Np = length of P
    // sp.Npxy = dimsPxy[0];   // Npxy = length of first dimension of P

    return sp;
}


/*  The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* Validate the inputs and outputs */
    if(nrhs!=2)
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:invalidNumInputs",
                "Two inputs required.");
    else if(nlhs > 1)
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:maxlhs",
                "Too many output arguments.");

    /* Check for proper inputs */
    const mxArray *S_in = prhs[0];
    const mxArray *sp_in = prhs[1];
    if(!mxIsSingle(S_in) || mxGetNumberOfDimensions(S_in) != 3)
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:inputNotDouble",
                "First input (S) must be a 3-D array of class single.");
    else if(!mxIsStruct(sp_in))
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:inputNotStruct",
                "Second input (sp) must be a structure.");

    /* Assign the inputs */
    single *S = mxGetData( S_in );
    simParam sp = parseSimParam( sp_in );

    /* Create output array */
    const mwSize *dims = mxGetDimensions(S_in);
    mxArray *Snext_out = mxCreateNumericArray( mxGetNumberOfDimensions(S_in),
                                          mxGetDimensions(S_in),
                                          mxSINGLE_CLASS, mxREAL );
    single *Snext = mxGetData( Snext_out );

    /* Invoke the computation */
    if( sp.useRK4 )
        rk4Step( Snext, S, sp );
    else
        eulerStep( Snext, S, sp );

    /* Return the output array */
    plhs[0] = Snext_out;
    return;
}

/* ==========================================================================
 * odeStepComp.c
 *
 * This is the top level for compiling odeStepComp mex-Funtion.
 *
 * Project:     Massive ODE Solver
 * Author:      Aamir Ahmed Khan (akhan3@nd.edu)
 * Copyright:   2012, University of Notre Dame
 *==========================================================================*/


#include "mex.h"
#include "string.h"
#include "assert.h"
#include "math.h"

typedef float single;

/* Define the simParam structure */
typedef struct {
    single tf, dt;
    int Nt, Ny, Nx;
    single Ms, gamma, alpha;
    single *cCoupl, *cDemag;
    int useRK4, useGPU;
} simParam;


/*! Implements the cross product, C = AxB.
 *  All vectors must be 3 dimensional. */
void cross(single *C, const single *A, const single *B) {
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
}


/*! Evaluates the Landau-Lifshitz-Gilbert (LLG) vector ODE.
 *  \param Mprime Return pointer for computed dM/dt vector
 *  \param sp Simulation parameters
 *  \param M Magnetization vector
 *  \param H Field vector acting on M */
void LLG( single *Mprime, simParam sp, single *M, single *H )
{
    single MxH[3], MxMxH[3];
    cross(MxH, M, H);       // MxH
    cross(MxMxH, M, MxH);   // Mx(MxH)
    Mprime[0] = -sp.gamma * MxH[0] - (sp.alpha*sp.gamma/sp.Ms) * MxMxH[0];
    Mprime[1] = -sp.gamma * MxH[1] - (sp.alpha*sp.gamma/sp.Ms) * MxMxH[1];
    Mprime[2] = -sp.gamma * MxH[2] - (sp.alpha*sp.gamma/sp.Ms) * MxMxH[2];
}


/*! Takes one ODE step by Euler's method.
 *  \param H Return pointer for total effective field
 *  \param sp Simulation parameters
 *  \param M Current state
 *  \param Hext Current external field */
void Hfield( single *H, simParam sp, single *M, single *Hext )
{
    /* TODO: Some optimization room here */
    int Nxy = sp.Ny * sp.Nx;    // size of array
    /* Start with the external field */
    memcpy(H, Hext, 3*Nxy*sizeof(single));

    /* Iterate over all the dots */
    for( int i = 0; i < sp.Ny*sp.Nx; ++i ) {
        /* Interaction with itself (Demagnetizating field) */
        H[i*3+0] += -sp.cDemag[0] * M[i*3+0];
        H[i*3+1] += -sp.cDemag[1] * M[i*3+1];
        H[i*3+2] += -sp.cDemag[2] * M[i*3+2];
    }

    /* Iterate over all the dots except for the boundary,
     * in column-major-order, as the data is from MATLAB */
    for( int x = 1; x < sp.Nx-1; ++x ) {        // exclude boundary dots
        for( int y = 1; y < sp.Ny-1; ++y ) {    // exclude boundary dots
            /* TODO: resolve the ambiguity in matrix style axes and xy-coord axes */
            single *Mtop = &M[ (y+1)*3+ x   *sp.Ny*3 ]; // +y (a/c to xy-cood axes)
            single *Mbot = &M[ (y-1)*3+ x   *sp.Ny*3 ]; // -y
            single *Mrig = &M[  y   *3+(x+1)*sp.Ny*3 ]; // +x
            single *Mlef = &M[  y   *3+(x-1)*sp.Ny*3 ]; // -x
            /* Add the coupling field now */
            H[y*3+x*sp.Ny*3+0] += sp.cCoupl[0] * ( Mtop[0] + Mbot[0] + Mrig[0] + Mlef[0] );
            H[y*3+x*sp.Ny*3+1] += sp.cCoupl[1] * ( Mtop[1] + Mbot[1] + Mrig[1] + Mlef[1] );
            H[y*3+x*sp.Ny*3+2] += sp.cCoupl[2] * ( Mtop[2] + Mbot[2] + Mrig[2] + Mlef[2] );
        }
    }
    /* Iterate over the dots on the boundaries */
    for( int x = 0; x < sp.Nx; ++x ) {
        for( int y = 0; y < sp.Ny; ++y ) {
            // if

        }
    }
}

/*! Takes one ODE step by Euler's method.
 *  \param Mnext State at the next time instant
 *  \param sp Simulation parameters
 *  \param M Current state
 *  \param Hext Current external field */
void eulerStep( single *Mnext, simParam sp, single *M, single *Hext )
{
    // for( int i = 0; i < sp.Ny*sp.Nx; ++i ) {
        // Mnext[i*3+0] = M[i*3+0];
        // Mnext[i*3+1] = M[i*3+1];
        // Mnext[i*3+2] = M[i*3+2];
    // }
    // return;
    // mexPrintf("Now running Euler's step...\n", sp.tf);
    int Nxy = sp.Ny * sp.Nx;    // size of array
    /* Compute the field */
    single *H = (single*)calloc( 3*Nxy, sizeof(single) );
    Hfield( H, sp, M, Hext );
    /* Advance to the next time instant */
    for( int i = 0; i < Nxy; ++i ) {
        single Mprime[3];
        LLG( Mprime, sp, &M[i*3], &H[i*3] );
        Mnext[i*3+0] = M[i*3+0] + Mprime[0] * sp.dt;
        Mnext[i*3+1] = M[i*3+1] + Mprime[1] * sp.dt;
        Mnext[i*3+2] = M[i*3+2] + Mprime[2] * sp.dt;
        /* re-normalize the M vectors */
        // single magnitude = sqrt( Mnext[i*3+0]*Mnext[i*3+0] +
                                 // Mnext[i*3+1]*Mnext[i*3+1] +
                                 // Mnext[i*3+2]*Mnext[i*3+2] );
        // Mnext[i*3+0] *= sp.Ms / magnitude;
        // Mnext[i*3+1] *= sp.Ms / magnitude;
        // Mnext[i*3+2] *= sp.Ms / magnitude;
    }
    /* clean-up */
    free(H);
}


/*! Converts Matlab's simParam structure to an equivalent C structure.
 *  \param spIn Matlab's version of structure
 *  \return Equivalent C structure */
simParam parseSimParam( const mxArray *spIn )
{
    /* fill up the sp structure */
    simParam sp = {
        .tf = mxGetScalar(mxGetField(spIn, 0, "tf")),
        .dt = mxGetScalar(mxGetField(spIn, 0, "dt")),
        .Nt = mxGetScalar(mxGetField(spIn, 0, "Nt")),
        .Ny = mxGetScalar(mxGetField(spIn, 0, "Ny")),
        .Nx = mxGetScalar(mxGetField(spIn, 0, "Nx")),
        .Ms = mxGetScalar(mxGetField(spIn, 0, "Ms")),
        .gamma = mxGetScalar(mxGetField(spIn, 0, "gamma")),
        .alpha = mxGetScalar(mxGetField(spIn, 0, "alpha")),
        .cCoupl = (single*)mxGetData(mxGetField(spIn, 0, "cCoupl")),
        .cDemag = (single*)mxGetData(mxGetField(spIn, 0, "cDemag")),
        .useRK4 = mxGetScalar(mxGetField(spIn, 0, "useRK4")),
        .useGPU = mxGetScalar(mxGetField(spIn, 0, "useGPU"))
    };
    return sp;
    mxArray *cDemagArray = mxGetField(spIn, 0, "cDemag");
    mexPrintf("mxGetNumberOfDimensions(cDemagArray) = %d\n", mxGetNumberOfDimensions(cDemagArray));
    mexPrintf("mxGetClassName(cDemagArray) = %s\n", mxGetClassName(cDemagArray));
    mexPrintf("sp.cDemag = [%g %g %g]\n", sp.cDemag[0], sp.cDemag[1], sp.cDemag[2]);
    // mexPrintf("sp.Ms = %g\n", sp.Ms);
    // mexPrintf("sp.alpha = %g\n", sp.alpha);
    // mexPrintf("sp.gamma = %g\n", sp.gamma);
    // /* Create 3D arrays for M and Hext */
    // assert(mxGetNumberOfDimensions(MIn) == 3);   // must be a 3D array
    // const mwSize *dims = mxGetDimensions(MIn);
    // single *M = (single*)malloc(dims[0]*dims[1]*dims[2]*sizeof(single));
    // single *Hext = (single*)malloc(dims[0]*dims[1]*dims[2]*sizeof(single));
    // memcpy((void*)M, mxGetData(MIn), dims[0]*dims[1]*dims[2]*sizeof(single));
    // memcpy((void*)Hext, mxGetData(HextIn), dims[0]*dims[1]*dims[2]*sizeof(single));
    // memcpy(mxGetData(MOut), mxGetData(MIn), dims[0]*dims[1]*dims[2]*sizeof(single));
}



/*  The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* Validate the inputs and outputs */
    if(nrhs!=3)
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:invalidNumInputs",
                "Three inputs required.");
    else if(nlhs > 1)
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:maxlhs",
                "Too many output arguments.");

    /* Check for proper inputs */
    const mxArray *spIn = prhs[0];
    const mxArray *MIn = prhs[1];
    const mxArray *HextIn = prhs[2];
    if(!mxIsStruct(spIn))
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:inputNotStruct",
                "First input (sp) must be a structure.");
    else if(!mxIsSingle(MIn) || mxGetNumberOfDimensions(MIn) != 3)
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:inputNotDouble",
                "Second input (M) must be a 3-D array of class single.");
    else if(!mxIsSingle(HextIn) || mxGetNumberOfDimensions(HextIn) != 3)
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:inputNotDouble",
                "Third input (HextIn) must be a 3-D array of class single.");

    /* Assign the inputs */
    simParam sp = parseSimParam( spIn );
    single *M = mxGetData( MIn );
    single *Hext = mxGetData( HextIn );

    /* Create output array */
    const mwSize *dims = mxGetDimensions(MIn);
    mxArray *MnextOut = mxCreateNumericArray( mxGetNumberOfDimensions(MIn),
                                          mxGetDimensions(MIn),
                                          mxSINGLE_CLASS, mxREAL );
    single *Mnext = mxGetData( MnextOut );

    /* Invoke the computation */
    if( sp.useRK4 ) {
        // rk4Step( Mnext, sp, M, Hext );
        fprintf(stderr, "ERROR: RK4 method not implemented yet\n");
    }
    else {
        eulerStep( Mnext, sp, M, Hext );
    }

    /* Return the output array */
    plhs[0] = MnextOut;
    return;
}

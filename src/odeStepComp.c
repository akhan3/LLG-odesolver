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
    single Ms, gamma, alpha;
    single Aexch, Kanis, *anisVec, *couplVec, *demagVec;
    single *bcMtop, *bcMbot, *bcMrig, *bcMlef;
    int useRK4, useGPU, preserveNorm;
} simParam;


/*! Re-normalize the M vector */
void reNormalize(single *M, const simParam sp) {
    single magnitude = sqrt( M[0]*M[0] + M[1]*M[1] + M[2]*M[2] );
    M[0] *= sp.Ms / magnitude;
    M[1] *= sp.Ms / magnitude;
    M[2] *= sp.Ms / magnitude;
}


/*! Implements the dot product, C = A.B
 *  Both the vectors must be 3 dimensional. */
inline single dot(const single *A, const single *B) {
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}


/*! Implements the cross product, C = AxB
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
    // Mprime[0] = -sp.gamma * (M[1]*H[2]-M[2]*H[1]) - (sp.alpha*sp.gamma/sp.Ms) * ( M[0]*(M[1]*H[1]+M[2]*H[2]) - H[0]*(M[1]*M[1]+M[2]*M[2]) );
    // Mprime[1] = -sp.gamma * (M[2]*H[0]-M[0]*H[2]) - (sp.alpha*sp.gamma/sp.Ms) * ( M[1]*(M[2]*H[2]+M[0]*H[0]) - H[1]*(M[2]*M[2]+M[0]*M[0]) );
    // Mprime[2] = -sp.gamma * (M[0]*H[1]-M[1]*H[0]) - (sp.alpha*sp.gamma/sp.Ms) * ( M[2]*(M[0]*H[0]+M[1]*H[1]) - H[2]*(M[0]*M[0]+M[1]*M[1]) );
}

/* TODO: Implement spatially varying parameters
 * DONE: Implement RK4 Solver
 * DONE: Implement the use of boundary conditions
 * DONE: Implement Exchange field
 */


/*! Calculates the effictive H-field caused by several phenomena.
 *  \param H Return pointer for total effective field
 *  \param sp Simulation parameters
 *  \param M Current state
 *  \param Hext Current external field */
void Hfield( single *H, simParam sp, single *M, single *Hext )
{
    const single mu0 = 4*M_PI*1e-7;   // Vacuum permeability in SI units [N/A^2]
    /* TODO: Maybe some optimization room here */
    int Nxy = sp.Ny * sp.Nx;    // size of array

    /* Start with the external field */
    memcpy(H, Hext, 3*Nxy*sizeof(single));

    /* Add demagnetization field (interaction only with itself) */
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int i = 0; i < Nxy; ++i ) {
        H[i*3+0] += -sp.demagVec[0] * M[i*3+0];
        H[i*3+1] += -sp.demagVec[1] * M[i*3+1];
        H[i*3+2] += -sp.demagVec[2] * M[i*3+2];
    }

    /* Add anisotropy field (interaction only with itself) */
    const single anisConstant = 2.0 * sp.Kanis / (mu0 * sp.Ms *sp.Ms);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int i = 0; i < Nxy; ++i ) {
        single Mdotk = dot(&M[i*3], sp.anisVec);
        H[i*3+0] += anisConstant * Mdotk * sp.anisVec[0];
        H[i*3+1] += anisConstant * Mdotk * sp.anisVec[1];
        H[i*3+2] += anisConstant * Mdotk * sp.anisVec[2];
    }

    /* Add exchange and coupling fields from nearest neighbours
     * Iterate over all the dots in column-major-order
     * since the data is coming from MATLAB
     * This step is most time consuming (bottleneck) */
    const single exchangeConstant = 2.0 * sp.Aexch / (mu0 * sp.Ms *sp.Ms);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int x = 0; x < sp.Nx; ++x ) {
        for( int y = 0; y < sp.Ny; ++y ) {
            /* DONE: Very clearly document that top means maximum y coordinate amd not first matrix row
             *       This code follows xy-coordinate axes convention */
            single *Mtop, *Mbot, *Mrig, *Mlef;
            single *M0 = &M[y*3 + x*sp.Ny*3];   // itself
            /* Pick neighbours */
            // if at top-most row, use top boundary condition
            Mtop = (y == sp.Ny-1) ? &sp.bcMtop[x*3] :
                                    &M[ (y+1)*3+ x   *sp.Ny*3 ]; // +y (wrt to xy-cood axes)
            // if at bottom-most row, use bottom boundary condition
            Mbot = (y == 0)       ? &sp.bcMbot[x*3] :
                                    &M[ (y-1)*3+ x   *sp.Ny*3 ]; // -y (wrt to xy-cood axes)
            // if at right-most column, use right boundary condition
            Mrig = (x == sp.Nx-1) ? &sp.bcMrig[y*3] :
                                    &M[  y   *3+(x+1)*sp.Ny*3 ]; // +x
            // if at left-most column, use left boundary condition
            Mlef = (x == 0)       ? &sp.bcMlef[y*3] :
                                    &M[  y   *3+(x-1)*sp.Ny*3 ]; // -x
            /* Laplacian */
            single L_Mx = (Mrig[0] - 2*M0[0] + Mlef[0]) / (sp.dx*sp.dx) +
                          (Mtop[0] - 2*M0[0] + Mbot[0]) / (sp.dy*sp.dy);
            single L_My = (Mrig[1] - 2*M0[1] + Mlef[1]) / (sp.dx*sp.dx) +
                          (Mtop[1] - 2*M0[1] + Mbot[1]) / (sp.dy*sp.dy);
            single L_Mz = (Mrig[2] - 2*M0[2] + Mlef[2]) / (sp.dx*sp.dx) +
                          (Mtop[2] - 2*M0[2] + Mbot[2]) / (sp.dy*sp.dy);
            /* Now add these two field componets */
            H[y*3+x*sp.Ny*3+0] += exchangeConstant * L_Mx
                    + sp.couplVec[0] * (Mtop[0] + Mbot[0] + Mrig[0] + Mlef[0]);
            H[y*3+x*sp.Ny*3+1] += exchangeConstant * L_My
                    + sp.couplVec[1] * (Mtop[1] + Mbot[1] + Mrig[1] + Mlef[1]);
            H[y*3+x*sp.Ny*3+2] += exchangeConstant * L_Mz
                    + sp.couplVec[2] * (Mtop[2] + Mbot[2] + Mrig[2] + Mlef[2]);
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
    int Nxy = sp.Ny * sp.Nx;    // size of array
    /* Compute the field */
    /* TODO: Calloc not necessary */
    single *H = (single*)calloc( 3*Nxy, sizeof(single) );
    Hfield( H, sp, M, Hext );
    /* Advance to the next time instant */
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int i = 0; i < Nxy; ++i ) {
        single Mprime[3];
        LLG( Mprime, sp, &M[i*3], &H[i*3] );
        Mnext[i*3+0] = M[i*3+0] + Mprime[0] * sp.dt;
        Mnext[i*3+1] = M[i*3+1] + Mprime[1] * sp.dt;
        Mnext[i*3+2] = M[i*3+2] + Mprime[2] * sp.dt;
        /* re-normalize the M vectors */
        if(sp.preserveNorm)
            reNormalize(&Mnext[i*3], sp);
    }
    /* clean-up */
    free(H);
}


/*! Takes one ODE step by Runge-Kutta's 4th order method.
 *  \param Mnext State at the next time instant
 *  \param sp Simulation parameters
 *  \param M Current state
 *  \param Hext Current external field */
void rk4Step( single *Mnext, simParam sp, single *M, single *Hext )
{
    /* allocate memory for field and slope values */
    /* TODO: Calloc not necessary */
    int Nxy = sp.Ny * sp.Nx;    // size of array
    single *H = (single*)calloc( 3*Nxy, sizeof(single) );
    single *slope = (single*)calloc( 3*Nxy, sizeof(single) );
    /* For k1  */
    Hfield( H, sp, M, Hext );   // Compute the field
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int i = 0; i < Nxy; ++i ) {
        single Mprime[3];
        LLG( Mprime, sp, &M[i*3], &H[i*3] );
        Mnext[i*3+0] = M[i*3+0] + Mprime[0] * .5*sp.dt;
        Mnext[i*3+1] = M[i*3+1] + Mprime[1] * .5*sp.dt;
        Mnext[i*3+2] = M[i*3+2] + Mprime[2] * .5*sp.dt;
        if(sp.preserveNorm)
            reNormalize(&Mnext[i*3], sp);
        slope[i*3+0] = Mprime[0] / 6.0;
        slope[i*3+1] = Mprime[1] / 6.0;
        slope[i*3+2] = Mprime[2] / 6.0;
    }
    /* For k2 */
    /* TODO: small error because not interpolating Hext(t+dt/2) */
    Hfield( H, sp, Mnext, Hext );   // Compute the field
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int i = 0; i < Nxy; ++i ) {
        single Mprime[3];
        LLG( Mprime, sp, &Mnext[i*3], &H[i*3] );
        Mnext[i*3+0] = M[i*3+0] + Mprime[0] * .5*sp.dt;
        Mnext[i*3+1] = M[i*3+1] + Mprime[1] * .5*sp.dt;
        Mnext[i*3+2] = M[i*3+2] + Mprime[2] * .5*sp.dt;
        if(sp.preserveNorm)
            reNormalize(&Mnext[i*3], sp);
        slope[i*3+0] += Mprime[0] / 3.0;
        slope[i*3+1] += Mprime[1] / 3.0;
        slope[i*3+2] += Mprime[2] / 3.0;
    }
    /* For k3 */
    /* TODO: small error because not interpolating Hext(t+dt/2) */
    Hfield( H, sp, Mnext, Hext );   // Compute the field
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int i = 0; i < Nxy; ++i ) {
        single Mprime[3];
        LLG( Mprime, sp, &Mnext[i*3], &H[i*3] );
        Mnext[i*3+0] = M[i*3+0] + Mprime[0] * sp.dt;
        Mnext[i*3+1] = M[i*3+1] + Mprime[1] * sp.dt;
        Mnext[i*3+2] = M[i*3+2] + Mprime[2] * sp.dt;
        if(sp.preserveNorm)
            reNormalize(&Mnext[i*3], sp);
        slope[i*3+0] += Mprime[0] / 3.0;
        slope[i*3+1] += Mprime[1] / 3.0;
        slope[i*3+2] += Mprime[2] / 3.0;
    }
    /* For k4 */
    Hfield( H, sp, Mnext, Hext );   // Compute the field
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int i = 0; i < Nxy; ++i ) {
        single Mprime[3];
        LLG( Mprime, sp, &Mnext[i*3], &H[i*3] );
        slope[i*3+0] += Mprime[0] / 6.0;
        slope[i*3+1] += Mprime[1] / 6.0;
        slope[i*3+2] += Mprime[2] / 6.0;
        /* Fianlly advance to the next time instant */
        Mnext[i*3+0] = M[i*3+0] + slope[i*3+0] * sp.dt;
        Mnext[i*3+1] = M[i*3+1] + slope[i*3+1] * sp.dt;
        Mnext[i*3+2] = M[i*3+2] + slope[i*3+2] * sp.dt;
        /* re-normalize the M vectors */
        if(sp.preserveNorm)
            reNormalize(&Mnext[i*3], sp);
    }
    /* clean-up */
    free(H);
    free(slope);
}


/*! Converts Matlab's simParam structure to an equivalent C structure.
 *  \param spIn Matlab's version of structure
 *  \return Equivalent C structure */
simParam parseSimParam( const mxArray *spIn )
{
    /* get the pointer to the nested boundary-condition structure */
    mxArray *bc = mxGetField(spIn, 0, "boundCond");
    single *bcMtop = (single*)mxGetData( mxGetField(bc, 0, "Mtop") );

    /* fill up the sp structure */
    simParam sp = {
        .dt = mxGetScalar(mxGetField(spIn, 0, "dt")),
        .Nt = mxGetScalar(mxGetField(spIn, 0, "Nt")),
        .Ny = mxGetScalar(mxGetField(spIn, 0, "Ny")),
        .Nx = mxGetScalar(mxGetField(spIn, 0, "Nx")),
        .dy = mxGetScalar(mxGetField(spIn, 0, "dy")),
        .dx = mxGetScalar(mxGetField(spIn, 0, "dx")),

        .Ms = mxGetScalar(mxGetField(spIn, 0, "Ms")),
        .gamma = mxGetScalar(mxGetField(spIn, 0, "gamma")),
        .alpha = mxGetScalar(mxGetField(spIn, 0, "alpha")),
        .Aexch = mxGetScalar(mxGetField(spIn, 0, "Aexch")),
        .Kanis = mxGetScalar(mxGetField(spIn, 0, "Kanis")),
        .anisVec = (single*)mxGetData(mxGetField(spIn, 0, "anisVec")),
        .couplVec = (single*)mxGetData(mxGetField(spIn, 0, "couplVec")),
        .demagVec = (single*)mxGetData(mxGetField(spIn, 0, "demagVec")),

        .bcMtop = (single*)mxGetData( mxGetField(bc, 0, "Mtop") ),
        .bcMbot = (single*)mxGetData( mxGetField(bc, 0, "Mbot") ),
        .bcMrig = (single*)mxGetData( mxGetField(bc, 0, "Mrig") ),
        .bcMlef = (single*)mxGetData( mxGetField(bc, 0, "Mlef") ),

        .useRK4 = mxGetScalar(mxGetField(spIn, 0, "useRK4")),
        .useGPU = mxGetScalar(mxGetField(spIn, 0, "useGPU")),
        .preserveNorm = mxGetScalar(mxGetField(spIn, 0, "preserveNorm"))
    };
    return sp;
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
    if( sp.useRK4 )
        rk4Step( Mnext, sp, M, Hext );
    else
        eulerStep( Mnext, sp, M, Hext );

    /* Return the output array */
    plhs[0] = MnextOut;
    return;
}

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
    int numM;
    single *P;
    single *Pxy;
    single *bcMtop, *bcMbot, *bcMrig, *bcMlef;
    int useRK4, useGPU, preserveNorm;
} simParam;


/*! Re-normalize the M vector */
void reNormalize(single *M, const simParam sp) {
    single Ms = sp.P[0];
    single magnitude = sqrt( M[0]*M[0] + M[1]*M[1] + M[2]*M[2] );
    M[0] *= Ms / magnitude;
    M[1] *= Ms / magnitude;
    M[2] *= Ms / magnitude;
}


/*! Compensate the M vector */
void compensateM(single *M, single *Mprev, int ixy, const simParam sp) {
    if(sp.preserveNorm)
        reNormalize(M, sp);
    // check for frozen locations
    int frozen = (int)sp.Pxy[0*sp.Ny*sp.Nx + ixy];
    if(frozen) {
        M[0] = Mprev[0];
        M[1] = Mprev[1];
        M[2] = Mprev[2];
    }
}


/*! Evaluates the slopes from the vector differential equation.
 *  \param Mprime Return pointer for computed dM/dt vector
 *  \param sp Simulation parameters
 *  \param M Magnetization vector
 *  \param H Field vector acting on M */
void differentiateM( single *Mprime, int ixy, simParam sp, single *M, single *H )
{
    single Ms = sp.P[0];
    single gamma = sp.P[1];
    single alpha = sp.Pxy[1*sp.Ny*sp.Nx + ixy];
    Mprime[0] = -gamma * (M[1]*H[2]-M[2]*H[1]) - (alpha*gamma/Ms) * ( M[0]*(M[1]*H[1]+M[2]*H[2]) - H[0]*(M[1]*M[1]+M[2]*M[2]) );
    Mprime[1] = -gamma * (M[2]*H[0]-M[0]*H[2]) - (alpha*gamma/Ms) * ( M[1]*(M[2]*H[2]+M[0]*H[0]) - H[1]*(M[2]*M[2]+M[0]*M[0]) );
    Mprime[2] = -gamma * (M[0]*H[1]-M[1]*H[0]) - (alpha*gamma/Ms) * ( M[2]*(M[0]*H[0]+M[1]*H[1]) - H[2]*(M[0]*M[0]+M[1]*M[1]) );
}


/*! Calculates the effictive H-field caused by several phenomena.
 *  \param H Return pointer for total effective field
 *  \param sp Simulation parameters
 *  \param M Current state
 *  \param Hext Current external field */
void Hfield( single *H, simParam sp, single *M )
{
    /* TODO: Maybe some optimization room in this function */
    const single mu0 = 4*M_PI*1e-7;   // Vacuum permeability in SI units [N/A^2]
    int Nxy = sp.Ny * sp.Nx;    // size of array

    /* Start with the external field */
    // memcpy(H, Hext, 3*Nxy*sizeof(single));
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int ixy = 0; ixy < Nxy; ++ixy ) {
        single HextX = sp.Pxy[2*sp.Ny*sp.Nx + ixy];
        single HextY = sp.Pxy[3*sp.Ny*sp.Nx + ixy];
        single HextZ = sp.Pxy[4*sp.Ny*sp.Nx + ixy];
        H[ixy*3+0] = HextX;
        H[ixy*3+1] = HextY;
        H[ixy*3+2] = HextZ;
    }

    /* Add demagnetization field (interaction only with itself) */
    single demagX = sp.P[11];
    single demagY = sp.P[12];
    single demagZ = sp.P[13];
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int ixy = 0; ixy < Nxy; ++ixy ) {
        H[ixy*3+0] += -demagX * M[ixy*3+0];
        H[ixy*3+1] += -demagY * M[ixy*3+1];
        H[ixy*3+2] += -demagZ * M[ixy*3+2];
    }

    /* Add anisotropy field (interaction only with itself) */
    single Ms = sp.P[0];
    single Kanis = sp.P[4];
    single anisX = sp.P[5];
    single anisY = sp.P[6];
    single anisZ = sp.P[7];
    const single anisConstant = 2.0 * Kanis / (mu0 * Ms*Ms);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int ixy = 0; ixy < Nxy; ++ixy ) {
        // single Mdotk = dot(&M[ixy*3], sp.anisVec);
        single Mdotk = M[ixy*3+0]*anisX + M[ixy*3+1]*anisY + M[ixy*3+2]*anisZ;
        H[ixy*3+0] += anisConstant * Mdotk * anisX;
        H[ixy*3+1] += anisConstant * Mdotk * anisY;
        H[ixy*3+2] += anisConstant * Mdotk * anisZ;
    }

    /* Add exchange and coupling fields from nearest neighbours
     * Iterate over all the dots in column-major-order
     * since the data is coming from MATLAB
     * This step is most time consuming (bottleneck) */
    single Aexch = sp.P[3];
    single couplX = sp.P[8];
    single couplY = sp.P[9];
    single couplZ = sp.P[10];
    const single exchangeConstant = 2.0 * Aexch / (mu0 * Ms*Ms);
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
                    + couplX * (Mtop[0] + Mbot[0] + Mrig[0] + Mlef[0]);
            H[y*3+x*sp.Ny*3+1] += exchangeConstant * L_My
                    + couplY * (Mtop[1] + Mbot[1] + Mrig[1] + Mlef[1]);
            H[y*3+x*sp.Ny*3+2] += exchangeConstant * L_Mz
                    + couplZ * (Mtop[2] + Mbot[2] + Mrig[2] + Mlef[2]);
        }
    }
}


/*! Takes one ODE step by Euler's method.
 *  \param Mnext State at the next time instant
 *  \param sp Simulation parameters
 *  \param M Current state
 *  \param Hext Current external field */
void eulerStep( single *Mnext, simParam sp, single *M )
{
    int Nxy = sp.Ny * sp.Nx;    // size of array
    /* Compute the field */
    /* Allocate memory for as many auxiliary H-fields as necessary */
    single *H = (single*)calloc( sp.numM*Nxy, sizeof(single) ); // only 1 here
    Hfield( H, sp, M );
    /* Advance to the next time instant */
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int ixy = 0; ixy < Nxy; ++ixy ) {
        single Mprime[3];
        differentiateM( Mprime, ixy, sp, &M[ixy*3], &H[ixy*3] );
        Mnext[ixy*3+0] = M[ixy*3+0] + Mprime[0] * sp.dt;
        Mnext[ixy*3+1] = M[ixy*3+1] + Mprime[1] * sp.dt;
        Mnext[ixy*3+2] = M[ixy*3+2] + Mprime[2] * sp.dt;
        /* compensate the M vector after update */
        compensateM(&Mnext[ixy*3], &M[ixy*3], ixy, sp);
    }
    /* clean-up */
    free(H);
}

/* TODO: RK4 is messed up. Priority */
/*! Takes one ODE step by Runge-Kutta's 4th order method.
 *  \param Mnext State at the next time instant
 *  \param sp Simulation parameters
 *  \param M Current state
 *  \param Hext Current external field */
void rk4Step( single *Mnext, simParam sp, single *M )
{
    /* allocate memory for field and slope values */
    /* TODO: Calloc not necessary */
    int Nxy = sp.Ny * sp.Nx;    // size of array
    single *H = (single*)calloc( sp.numM*Nxy, sizeof(single) );
    single *slope = (single*)calloc( sp.numM*Nxy, sizeof(single) );
    /* For k1  */
    Hfield( H, sp, M );   // Compute the field
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int ixy = 0; ixy < Nxy; ++ixy ) {
        single Mprime[3];
        differentiateM( Mprime, ixy, sp, &M[ixy*3], &H[ixy*3] );
        Mnext[ixy*3+0] = M[ixy*3+0] + Mprime[0] * .5*sp.dt;
        Mnext[ixy*3+1] = M[ixy*3+1] + Mprime[1] * .5*sp.dt;
        Mnext[ixy*3+2] = M[ixy*3+2] + Mprime[2] * .5*sp.dt;
        /* compensate the M vector after update */
        compensateM(&Mnext[ixy*3], &M[ixy*3], ixy, sp);
        slope[ixy*3+0] = Mprime[0] / 6.0;
        slope[ixy*3+1] = Mprime[1] / 6.0;
        slope[ixy*3+2] = Mprime[2] / 6.0;
    }
    /* For k2 */
    /* TODO: small error because not interpolating Hext(t+dt/2) */
    Hfield( H, sp, Mnext );   // Compute the field
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int ixy = 0; ixy < Nxy; ++ixy ) {
        single Mprime[3];
        differentiateM( Mprime, ixy, sp, &Mnext[ixy*3], &H[ixy*3] );
        Mnext[ixy*3+0] = M[ixy*3+0] + Mprime[0] * .5*sp.dt;
        Mnext[ixy*3+1] = M[ixy*3+1] + Mprime[1] * .5*sp.dt;
        Mnext[ixy*3+2] = M[ixy*3+2] + Mprime[2] * .5*sp.dt;
        /* compensate the M vector after update */
        compensateM(&Mnext[ixy*3], &M[ixy*3], ixy, sp);
        slope[ixy*3+0] += Mprime[0] / 3.0;
        slope[ixy*3+1] += Mprime[1] / 3.0;
        slope[ixy*3+2] += Mprime[2] / 3.0;
    }
    /* For k3 */
    /* TODO: small error because not interpolating Hext(t+dt/2) */
    Hfield( H, sp, Mnext );   // Compute the field
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int ixy = 0; ixy < Nxy; ++ixy ) {
        single Mprime[3];
        differentiateM( Mprime, ixy, sp, &Mnext[ixy*3], &H[ixy*3] );
        Mnext[ixy*3+0] = M[ixy*3+0] + Mprime[0] * sp.dt;
        Mnext[ixy*3+1] = M[ixy*3+1] + Mprime[1] * sp.dt;
        Mnext[ixy*3+2] = M[ixy*3+2] + Mprime[2] * sp.dt;
        /* compensate the M vector after update */
        compensateM(&Mnext[ixy*3], &M[ixy*3], ixy, sp);
        slope[ixy*3+0] += Mprime[0] / 3.0;
        slope[ixy*3+1] += Mprime[1] / 3.0;
        slope[ixy*3+2] += Mprime[2] / 3.0;
    }
    /* For k4 */
    Hfield( H, sp, Mnext );   // Compute the field
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( int ixy = 0; ixy < Nxy; ++ixy ) {
        single Mprime[3];
        differentiateM( Mprime, ixy, sp, &Mnext[ixy*3], &H[ixy*3] );
        slope[ixy*3+0] += Mprime[0] / 6.0;
        slope[ixy*3+1] += Mprime[1] / 6.0;
        slope[ixy*3+2] += Mprime[2] / 6.0;
        /* Fianlly advance to the next time instant */
        Mnext[ixy*3+0] = M[ixy*3+0] + slope[ixy*3+0] * sp.dt;
        Mnext[ixy*3+1] = M[ixy*3+1] + slope[ixy*3+1] * sp.dt;
        Mnext[ixy*3+2] = M[ixy*3+2] + slope[ixy*3+2] * sp.dt;
        /* compensate the M vector after update */
        compensateM(&Mnext[ixy*3], &M[ixy*3], ixy, sp);
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
        .numM = mxGetScalar(mxGetField(spIn, 0, "numM")),

        .P   = (single*)mxGetData(mxGetField(spIn, 0, "P")),
        .Pxy = (single*)mxGetData(mxGetField(spIn, 0, "Pxy")),

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
    if(nrhs!=2)
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:invalidNumInputs",
                "Two inputs required.");
    else if(nlhs > 1)
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:maxlhs",
                "Too many output arguments.");

    /* Check for proper inputs */
    const mxArray *spIn = prhs[0];
    const mxArray *MIn = prhs[1];
    if(!mxIsStruct(spIn))
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:inputNotStruct",
                "First input (sp) must be a structure.");
    else if(!mxIsSingle(MIn) || mxGetNumberOfDimensions(MIn) != 3)
        mexErrMsgIdAndTxt( "MyToolbox:odeStepComp:inputNotDouble",
                "Second input (M) must be a 3-D array of class single.");

    /* Assign the inputs */
    simParam sp = parseSimParam( spIn );
    single *M = mxGetData( MIn );

    /* Create output array */
    const mwSize *dims = mxGetDimensions(MIn);
    mxArray *MnextOut = mxCreateNumericArray( mxGetNumberOfDimensions(MIn),
                                          mxGetDimensions(MIn),
                                          mxSINGLE_CLASS, mxREAL );
    single *Mnext = mxGetData( MnextOut );

    /* Invoke the computation */
    if( sp.useRK4 )
        rk4Step( Mnext, sp, M );
    else
        eulerStep( Mnext, sp, M );

    /* Return the output array */
    plhs[0] = MnextOut;
    return;
}

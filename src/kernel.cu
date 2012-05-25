/* ==========================================================================
 * kernel.cu
 *
 * This is the CUDA source file that is compiled separately by MATLAB routine
 * and then linked to by MEX.
 * User must edit by this file to model different problems.
 *
 * Project:     Massive ODE Solver
 * Author:      Aamir Ahmed Khan (akhan3@nd.edu)
 * Copyright:   2012, University of Notre Dame
 * =========================================================================*/

#include "globals.h"
#include <cutil_inline.h>


/* Device constant memory */
__device__ __constant__ single *bcMtop_c, *bcMbot_c, *bcMrig_c, *bcMlef_c;


/* IMPORTANT: Don't mess with this function!!! Just use it */
/*! Picks the correct neighbors based on the location
 *      Picks the neighbor if exists otherwise pick the corresponding boundary condition */
__device__
void pickNeighbors_d( single **S_top, single **S_bot, single **S_rig, single **S_lef,
                      single *bcMtop, single *bcMbot, single *bcMrig, single *bcMlef,
                      single *S, const int ix, const int iy,
                      const int Nx, const int Ny, const int Ns )
{
    /* Here, top means maximum iy coordinate, not the first matrix row
     *       This code follows xy-coordinate axes convention */
    // if at top-most row, use top boundary condition
    *S_top = (iy == Ny-1)  ?  &bcMtop[ix*Ns]  :  &S[(iy+1)*Ns+ ix   *Ny*Ns]; // +y (wrt to xy-cood axes)
    // if at bottom-most row, use bottom boundary condition
    *S_bot = (iy == 0)        ?  &bcMbot[ix*Ns]  :  &S[(iy-1)*Ns+ ix   *Ny*Ns]; // -y (wrt to xy-cood axes)
    // if at right-most column, use right boundary condition
    *S_rig = (ix == Nx-1)  ?  &bcMrig[iy*Ns]  :  &S[ iy   *Ns+(ix+1)*Ny*Ns]; // +x
    // if at left-most column, use left boundary condition
    *S_lef = (ix == 0)        ?  &bcMlef[iy*Ns]  :  &S[ iy   *Ns+(ix-1)*Ny*Ns]; // -x
}


/*! Kernel definition */
__global__ void
kernel( single *Sprime, single *S,
        single *bcMtop, single *bcMbot, single *bcMrig, single *bcMlef,
        const int Nx, const int Ny, const int Ns,
        const single dx, const single dy )
{
    const int n = blockIdx.x * blockDim.x + threadIdx.x;
    if(n >= Ny * Nx)
        return;

    // extract ix and iy indices from the single index n
    int iy = n % Ny;
    int ix = (n-iy) / Ny;

    /* Pick neighbours */
    single *S_top, *S_bot, *S_rig, *S_lef;  // pointers to neighbors
    // this function must be called to take care of boundary conditions and pick proper neighbors
    pickNeighbors_d(&S_top, &S_bot, &S_rig, &S_lef, bcMtop, bcMbot, bcMrig, bcMlef, S, ix, iy, Nx, Ny, Ns);


    /* Execute the differential equations */
    // Variable order : S[0] - Hx        S[1] - Hy       S[2] - Ez
    // Calculate Laplacian
    single *S_cen = &S[n * Ns];   // pointer to variable
    single DHx = -1.0f*(S_top[2] - 2.0f*S_cen[2] + S_bot[2]) / (dy*dy);
    single DHy = +1.0f*(S_lef[2] - 2.0f*S_cen[2] + S_rig[2]) / (dx*dx);
    single DEz = -1.0f*( (S_lef[1] - 2.0f*S_cen[1] + S_rig[1]) / (dx*dx) - (S_top[0] - 2.0f*S_cen[0] + S_bot[0]) / (dy*dy) );

    /* Assign the derivatives */
    Sprime[n*Ns+0] = DHx;
    Sprime[n*Ns+1] = DHy;
    Sprime[n*Ns+2] = DEz;
}


/*! Evaluates the slopes from the vector differential equation.
 *      This is the major function which the user has to modify to model their
 *      own problem.
 *  \param Sprime Return pointer for computed dS/dt vector field
 *  \param sp Simulation parameters
 *  \param S State vector */
void differentiateStateVectorField_d( single *Sprime, single *S, simParam sp )
{
    int Nxy = sp.Ny * sp.Nx;    // size of array

    /* set up and allocate device memory */
    cutilSafeCall( cudaSetDevice(sp.useGPUnum) );
    single *Sprime_d = NULL;
    single *S_d = NULL;
    single *bcMtop_d = NULL;
    single *bcMbot_d = NULL;
    single *bcMrig_d = NULL;
    single *bcMlef_d = NULL;
    cutilSafeCall( cudaMalloc( (void**)&Sprime_d,   sp.Ns*Nxy*sizeof(single) ) );
    cutilSafeCall( cudaMalloc( (void**)&S_d,        sp.Ns*Nxy*sizeof(single) ) );
    cutilSafeCall( cudaMalloc( (void**)&bcMtop_d,   sp.Ns*sp.Nx*sizeof(single) ) );
    cutilSafeCall( cudaMalloc( (void**)&bcMbot_d,   sp.Ns*sp.Nx*sizeof(single) ) );
    cutilSafeCall( cudaMalloc( (void**)&bcMrig_d,   sp.Ns*sp.Ny*sizeof(single) ) );
    cutilSafeCall( cudaMalloc( (void**)&bcMlef_d,   sp.Ns*sp.Ny*sizeof(single) ) );
    assert(Sprime_d != NULL && S_d != NULL && bcMtop_d != NULL && bcMbot_d != NULL && bcMrig_d != NULL && bcMlef_d != NULL);

    /* copy State Variables to device global memory */
    cutilSafeCall( cudaMemset( Sprime_d, 0, sp.Ns*Nxy*sizeof(single) ) );
    cutilSafeCall( cudaMemcpy( S_d, S, sp.Ns*Nxy*sizeof(single), cudaMemcpyHostToDevice ) );
    cutilSafeCall( cudaMemcpy( bcMtop_d, sp.bcMtop, sp.Ns*sp.Nx*sizeof(single), cudaMemcpyHostToDevice ) );
    cutilSafeCall( cudaMemcpy( bcMbot_d, sp.bcMbot, sp.Ns*sp.Nx*sizeof(single), cudaMemcpyHostToDevice ) );
    cutilSafeCall( cudaMemcpy( bcMrig_d, sp.bcMrig, sp.Ns*sp.Ny*sizeof(single), cudaMemcpyHostToDevice ) );
    cutilSafeCall( cudaMemcpy( bcMlef_d, sp.bcMlef, sp.Ns*sp.Ny*sizeof(single), cudaMemcpyHostToDevice ) );

    /* set up kernel parameters */
    #ifdef __DEVICE_EMULATION__
        #define DIM 64
    #else
        #define DIM 512
    #endif
    dim3 grid = ceil(Nxy / (single)DIM);
    dim3 threads(DIM, 1, 1);
    assert(threads.x <= DIM);    // max_threads_per_block

    /* launch the kernel */
    kernel <<<grid, threads>>> (Sprime_d, S_d, bcMtop_d, bcMbot_d, bcMrig_d, bcMlef_d, sp.Nx, sp.Ny, sp.Ns, sp.dx, sp.dy );

    /* copy Sprime_d (result of the kernel) back to host main memory */
    cutilSafeCall( cudaMemcpy( Sprime, Sprime_d, sp.Ns*Nxy*sizeof(single), cudaMemcpyDeviceToHost ) );

    /* clean-up */
    // must deallocate all the reserved memory, otherwise huge performance penalty!!!
    cutilSafeCall( cudaFree(Sprime_d) );
    cutilSafeCall( cudaFree(S_d) );
    cutilSafeCall( cudaFree(bcMtop_d) );
    cutilSafeCall( cudaFree(bcMbot_d) );
    cutilSafeCall( cudaFree(bcMrig_d) );
    cutilSafeCall( cudaFree(bcMlef_d) );
}

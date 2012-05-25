/* ==========================================================================
 * modelingEquations.c
 *
 * This is the C source file that is included in odeStepComp.c and needs to be
 * edited by the user to model different problems.
 *
 * Project:     Massive ODE Solver
 * Author:      Aamir Ahmed Khan (akhan3@nd.edu)
 * Copyright:   2012, University of Notre Dame
 * =========================================================================*/

#include "globals.h"

/*! Evaluates the slopes from the vector differential equation.
 *      This is the major function which the user has to modify to model their
 *      own problem.
 *  \param Sprime Return pointer for computed dS/dt vector field
 *  \param sp Simulation parameters
 *  \param S State vector */
void differentiateStateVectorField( single *Sprime, single *S, simParam sp )
{
// Execute the differential equations
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( int ix = 0; ix < sp.Nx; ++ix ) {
        for( int iy = 0; iy < sp.Ny; ++iy ) {
            /* Pick neighbours */
            single *S_top, *S_bot, *S_rig, *S_lef;  // pointers to neighbors
            // this function must be called to take care of boundary conditions and pick proper neighbors
            pickNeighbors(&S_top, &S_bot, &S_rig, &S_lef, S, sp, ix, iy );

            /* Execute the differential equations */
            // Variable order : S[0] - Hx        S[1] - Hy       S[2] - Ez
            // Calculate Laplacian
            single *S_cen = &S[iy*sp.Ns + ix*sp.Ny*sp.Ns];   // pointer to variable
            single DHx = -1.0*(S_top[2] - 2.0*S_cen[2] + S_bot[2]) / (sp.dy*sp.dy);
            single DHy = +1.0*(S_lef[2] - 2.0*S_cen[2] + S_rig[2]) / (sp.dx*sp.dx);
            single DEz = -1.0*( (S_lef[1] - 2.0*S_cen[1] + S_rig[1]) / (sp.dx*sp.dx) - (S_top[0] - 2.0*S_cen[0] + S_bot[0]) / (sp.dy*sp.dy) );

            /* Assign the derivatives */
            Sprime[iy*sp.Ns+ix*sp.Ny*sp.Ns+0] = DHx;
            Sprime[iy*sp.Ns+ix*sp.Ny*sp.Ns+1] = DHy;
            Sprime[iy*sp.Ns+ix*sp.Ny*sp.Ns+2] = DEz;
        }   // iy
    }   // ix
}




/*! Compensate the S vector after updating in the ODE stepper
 *      This is the smaller function which the user has to modify to model their own problem.
 *      Also note that this function operates point-wise at a particular location.
 *  \param S Return pointer for the state vector at that location
 *  \param Sprev pointer for state vector at that location at previous time instant.
 *                  Can be used to freeze the state at that location.
 *  \param sp Simulation parameters
 *  \param ixy single index number indicating the location on xy-plane */
void compensateStateVector( single *S, single *Sprev, const simParam sp, int ixy ) {
    // restore norm of S vector
    // int preserveNorm = (int)sp.P[14];
    // if(preserveNorm) {
        // single Ms = sp.P[0];
        // single magnitude = sqrt( S[0]*S[0] + S[1]*S[1] + S[2]*S[2] );
        // S[0] *= Ms / magnitude;
        // S[1] *= Ms / magnitude;
        // S[2] *= Ms / magnitude;
    // }
    // // check for frozen locations
    // int frozen = (int)sp.Pxy[ixy*sp.Npxy + 0];
    // if(frozen) {
        // S[0] = Sprev[0];
        // S[1] = Sprev[1];
        // S[2] = Sprev[2];
    // }
}

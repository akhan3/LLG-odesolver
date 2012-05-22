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


/*! Evaluates the slopes from the vector differential equation.
 *      This is the major function which the user has to modify to model their
 *      own problem.
 *  \param Sprime Return pointer for computed dS/dt vector field
 *  \param sp Simulation parameters
 *  \param S State vector */
void differentiateStateVectorField( single *Sprime, single *S, simParam sp )
{
    // /* TODO: Maybe some optimization room in this function */
    // const single mu0 = 4*M_PI*1e-7;   // Vacuum permeability in SI units [N/A^2]
    // int Nxy = sp.Ny * sp.Nx;    // size of array

    // /* Allocate memory for as many auxiliary H-fields as necessary
     // * but don't forget to deallocate at the end of this function! */
    // single *H = (single*)calloc( sp.Ns*Nxy, sizeof(single) );

    // /* Effective field computation */
    // {
        // /* Start with the external field */
        // #ifdef _OPENMP
        // #pragma omp parallel for
        // #endif
        // for( int ixy = 0; ixy < Nxy; ++ixy ) {
            // int i = ixy*sp.Ns;
            // single HextX = sp.Pxy[ixy*sp.Npxy + 2];
            // single HextY = sp.Pxy[ixy*sp.Npxy + 3];
            // single HextZ = sp.Pxy[ixy*sp.Npxy + 4];
            // H[i+0] = HextX;
            // H[i+1] = HextY;
            // H[i+2] = HextZ;
        // }

        // /* Add demagnetization field (interaction only with itself) */
        // single demagX = sp.P[11];
        // single demagY = sp.P[12];
        // single demagZ = sp.P[13];
        // #ifdef _OPENMP
        // #pragma omp parallel for
        // #endif
        // for( int ixy = 0; ixy < Nxy; ++ixy ) {
            // int i = ixy*sp.Ns;
            // H[i+0] += -demagX * S[i+0];
            // H[i+1] += -demagY * S[i+1];
            // H[i+2] += -demagZ * S[i+2];
        // }

        // /* Add anisotropy field (interaction only with itself) */
        // single Ms = sp.P[0];
        // single Kanis = sp.P[4];
        // single anisX = sp.P[5];
        // single anisY = sp.P[6];
        // single anisZ = sp.P[7];
        // const single anisConstant = 2.0 * Kanis / (mu0 * Ms*Ms);
        // #ifdef _OPENMP
        // #pragma omp parallel for
        // #endif
        // for( int ixy = 0; ixy < Nxy; ++ixy ) {
            // int i = ixy*sp.Ns;
            // single Mdotk = S[i+0]*anisX + S[i+1]*anisY + S[i+2]*anisZ;
            // H[i+0] += anisConstant * Mdotk * anisX;
            // H[i+1] += anisConstant * Mdotk * anisY;
            // H[i+2] += anisConstant * Mdotk * anisZ;
        // }

        // /* Add exchange and coupling fields from nearest neighbours
         // * Iterate over all the dots in column-major-order
         // * since the data is coming from MATLAB
         // * This step is most time consuming (bottleneck) */
        // single Aexch = sp.P[3];
        // single couplX = sp.P[8];
        // single couplY = sp.P[9];
        // single couplZ = sp.P[10];
        // const single exchangeConstant = 2.0 * Aexch / (mu0 * Ms*Ms);
        // #ifdef _OPENMP
        // #pragma omp parallel for
        // #endif
        // for( int ix = 0; ix < sp.Nx; ++ix ) {
            // for( int iy = 0; iy < sp.Ny; ++iy ) {
                // /* Pick neighbours */
                // single *S_top, *S_bot, *S_rig, *S_lef;  // pointers to neighbors
                // // this function must be called to take care of boundary conditions and pick proper neighbors
                // pickNeighbors(&S_top, &S_bot, &S_rig, &S_lef, S, sp, ix, iy );
                // /* Calculate Laplacian */
                // single *S_cen = &S[iy*sp.Ns + ix*sp.Ny*sp.Ns];   // pointer to itself
                // single L_Mx = (S_rig[0] - 2*S_cen[0] + S_lef[0]) / (sp.dx*sp.dx) +
                              // (S_top[0] - 2*S_cen[0] + S_bot[0]) / (sp.dy*sp.dy);
                // single L_My = (S_rig[1] - 2*S_cen[1] + S_lef[1]) / (sp.dx*sp.dx) +
                              // (S_top[1] - 2*S_cen[1] + S_bot[1]) / (sp.dy*sp.dy);
                // single L_Mz = (S_rig[2] - 2*S_cen[2] + S_lef[2]) / (sp.dx*sp.dx) +
                              // (S_top[2] - 2*S_cen[2] + S_bot[2]) / (sp.dy*sp.dy);
                // /* Add exchange field */
                // H[iy*sp.Ns+ix*sp.Ny*sp.Ns+0] += exchangeConstant * L_Mx;
                // H[iy*sp.Ns+ix*sp.Ny*sp.Ns+1] += exchangeConstant * L_My;
                // H[iy*sp.Ns+ix*sp.Ny*sp.Ns+2] += exchangeConstant * L_Mz;
                // /* Add nearest neighbor coupling field */
                // H[iy*sp.Ns+ix*sp.Ny*sp.Ns+0] += couplX * (S_top[0] + S_bot[0] + S_rig[0] + S_lef[0]);
                // H[iy*sp.Ns+ix*sp.Ny*sp.Ns+1] += couplY * (S_top[1] + S_bot[1] + S_rig[1] + S_lef[1]);
                // H[iy*sp.Ns+ix*sp.Ny*sp.Ns+2] += couplZ * (S_top[2] + S_bot[2] + S_rig[2] + S_lef[2]);
            // }   // iy
        // }   // ix
    // }   // Field computation ended

    // /* Execute the differential equations */
    // #ifdef _OPENMP
    // #pragma omp parallel for
    // #endif
    // for( int ixy = 0; ixy < Nxy; ++ixy ) {
        // int i = ixy*sp.Ns;
        // single Ms = sp.P[0];
        // single gamma = sp.P[1];
        // single alpha = sp.Pxy[ixy*sp.Npxy + 1];
        // Sprime[i+0] = -gamma * (S[i+1]*H[i+2]-S[i+2]*H[i+1]) - (alpha*gamma/Ms) * ( S[i+0]*(S[i+1]*H[i+1]+S[i+2]*H[i+2]) - H[i+0]*(S[i+1]*S[i+1]+S[i+2]*S[i+2]) );
        // Sprime[i+1] = -gamma * (S[i+2]*H[i+0]-S[i+0]*H[i+2]) - (alpha*gamma/Ms) * ( S[i+1]*(S[i+2]*H[i+2]+S[i+0]*H[i+0]) - H[i+1]*(S[i+2]*S[i+2]+S[i+0]*S[i+0]) );
        // Sprime[i+2] = -gamma * (S[i+0]*H[i+1]-S[i+1]*H[i+0]) - (alpha*gamma/Ms) * ( S[i+2]*(S[i+0]*H[i+0]+S[i+1]*H[i+1]) - H[i+2]*(S[i+0]*S[i+0]+S[i+1]*S[i+1]) );
    // }

    // /* clean-up */
    // // must deallocate all the memory reserved for auxiliary fields, otherwise huge performance penalty!!!
    // free(H);


////////////////////////////////////////////////////////////////////////////////

// Execute the differential equations
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( int ix = 0; ix < sp.Nx; ++ix ) {
        for( int iy = 0; iy < sp.Ny; ++iy ) {
            // Pick neighbours
            single *S_top, *S_bot, *S_rig, *S_lef;  // pointers to neighbors
            // this function must be called to take care of boundary conditions and pick proper neighbors
            pickNeighbors(&S_top, &S_bot, &S_rig, &S_lef, S, sp, ix, iy );

            // Variable order : S[0] - Hx        S[1] - Hy       S[2] - Ez

            // Calculate Laplacian
            single *S_cen = &S[iy*sp.Ns + ix*sp.Ny*sp.Ns];   // pointer to variable
            single DHx = -1.0*(S_top[2] - 2.0*S_cen[2] + S_bot[2]) / (sp.dy*sp.dy);
            single DHy = +1.0*(S_lef[2] - 2.0*S_cen[2] + S_rig[2]) / (sp.dx*sp.dx);
            single DEz = -1.0*( (S_lef[1] - 2.0*S_cen[1] + S_rig[1]) / (sp.dx*sp.dx) - (S_top[0] - 2.0*S_cen[0] + S_bot[0]) / (sp.dy*sp.dy) );

            // And now it actually calls and solves the differential equation
            Sprime[iy*sp.Ns+ix*sp.Ny*sp.Ns+0] = DHx;
            Sprime[iy*sp.Ns+ix*sp.Ny*sp.Ns+1] = DHy;
            Sprime[iy*sp.Ns+ix*sp.Ny*sp.Ns+2] = DEz;
        }   // iy
    }   // ix

////////////////////////////////////////////////////////////////////////////////

}




/*! Compensate the S vector after updating in the ODE stepper
 *      This is the smaller function which the user has to modify to model their own problem.
 *      Also note that this function operates point-wise at a particular location.
 *  \param S Return pointer for the state vector at that location
 *  \param Sprev pointer for state vector at that location at previous time instant.
 *                  Can be used to freeze the state at that location.
 *  \param sp Simulation parameters
 *  \param ixy single index number indicating the location on xy-plane */
void compensateStateVector(single *S, single *Sprev, const simParam sp, int ixy ) {
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

/* TBD: Do we need to change the operation of compensateStateVector function
 *      from point-wise to iterate over all?
 *      This will result in some performance penalty!!! */

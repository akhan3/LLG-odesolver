clc
clear

%===============================================================================
%% Add mexDir to current path and call the function to compile MEX file
%===============================================================================
    mexDir = 'src';
    addpath(mexDir);    % add to path
    mexSetup(mexDir);   % compile the MEX file


%===============================================================================
%% Simulation parameters
%===============================================================================
    % Stuff all the parameters in a simParam structure
    % Some book-keeping information
    sp.simName = 'Untitled Simulation';
    sp.startTimeStamp = clock;  % record the start wall-clock-time
    sp.finishedWithSuccess = 0; % initially mark as unsuccessful
    % simulation time
    sp.ti = 0;      % initial time [s]
    sp.tf = .5e-9;  % final time [s]
    sp.dt = 1e-13;  % time step [s]
    sp.t = [sp.ti:sp.dt:sp.tf]; % time array [s]
    sp.Nt = length(sp.t);       % number of time points
    % number and size of the dots
        % Nx*Ny MUST NOT exceed (5000)^2
        % Total GPU memory = 2.8177e+09 Bytes
        % Each dot requires 3*4 Bytes
        % 2.8177e+09/12/10 to accommodate S(t),S(t+1),Hext,... on the GPU
    sp.Ny = 60;     % #rows of dots in the plane
    sp.Nx = 80;     % #columns of dots in the plane
    sp.dy = 100e-9; % y-length of a dot [m]
    sp.dx = 100e-9; % x-width of a dot [m]
    sp.Ns = 3;    % numer of state variables at each mesh point (size of state vector)

    % ==========================================================================
    % material parameters
    % ==========================================================================
    % Paramters that are uniform over the mesh
    sp.P = single(zeros(1,1));  % initialize to be of singe precision
    sp.P( 1) = 8.6e5;       % Saturation magnetization [A/m]
    sp.P( 2) = 2.21e5;      % gamma. Gyromagnetic Ratio [1/(A/m/s)]
    sp.P( 3) = 0.05;        % Damping factor [dim-less]
    sp.P( 4) = 1.3e-11;     % Exchange constant [J/A]
    sp.P( 5) = 1e5;         % Anisotropy constant [J/m^3]
    sp.P( 6) = 0;   % x-component of unit vector defining anisotropy axes [dim-less]
    sp.P( 7) = 0;   % y-component of unit vector defining anisotropy axes [dim-less]
    sp.P( 8) = 1;   % z-component of unit vector defining anisotropy axes [dim-less]
    % Coupling coefficient [dim-less]
        % It is defined to be negative here and used as it is in field calculation
    sp.P( 9) = -.2;    % x-coupling coefficient [dim-less]
    sp.P(10) = -.2;    % y-coupling coefficient [dim-less]
    sp.P(11) = -.2;    % z-coupling coefficient [dim-less]
    % Demag factor [dim-less]
        % It is defined to be positive here and used as negative in field calculation
    sp.P(12) = .4;      % x-demag factor [dim-less]
    sp.P(13) = .4;      % y-demag factor [dim-less]
    sp.P(14) = .2;      % z-demag factor [dim-less]
    sp.P(15) = 1;       % preserveNorm. If 1, all S vectors will be
                        %   re-normalized to Ms after the ODE-step is taken,

    % Paramters that vary over the mesh. P(y,x)
    sp.Pxy = single(zeros(1, sp.Ny,sp.Nx));
    sp.Pxy(1, :,:) = 0 * ones(sp.Ny,sp.Nx);      % frozenMask. S will be kept frozen where this variable is 1
    %sp.Pxy(1, 5:20,10:45) = 1;
    sp.Pxy(2, :,:) = .05 * ones(sp.Ny,sp.Nx);    % Position dependent damping factor
    % External field. Time dependence will be defined in time marching loop
    sp.Pxy(3, :,:) = 0;  % x-component of Hext [A/m]
    sp.Pxy(4, :,:) = 0;  % y-component of Hext [A/m]
    sp.Pxy(5, :,:) = 0;  % z-component of Hext [A/m]
    sp.Np = length(sp.P);         % numer of spatially uniform paramters
    sp.Npxy = size(sp.Pxy, 1);    % numer of spatially varying paramters
    % ==========================================================================

    % ODE Solver selection
    sp.useRK4 = 0;  % if 1, RK4-ODE-solver will be used, otherwise Euler's
    sp.useGPU = 0;  % if 1, GPU will be used


%===============================================================================
%% Allocate bulk data for state vectors
    % Dimensions of S must be treated as follows
    % S(i,y,x,t)    =>  i = 1:Ns for the index of variable in state vector
    %                   y = 1:Ny for y-index (row number) in plane
    %                   x = 1:Nx for x-index (column number) in plane
    %                   t = 1:Nt for time index
    % It is better to define them as single precision,
    % otherwise they will be converted by the validateSimParam routine anyway.
%===============================================================================
    S = single(zeros(sp.Ns,sp.Ny,sp.Nx,sp.Nt)); % can be huge in size


%%===============================================================================
%%% Define External field's time and spatial dependence
    %% Hext can either be defined completely beforehand here before entering
    %% time-marching loop. It can, of course, also be modified in the
    %% time-marching loop
%%===============================================================================
    %% This defines a 10GHz disc shaped oscillator near the center of the matrix
    %r = round(sp.Ny/6);
    %cx = round(sp.Nx/2);
    %cy = round(sp.Ny/3);
    %f = 10e9;
    %amplitude = 5;
    %sinwave = sin(2*pi*f*sp.t);
    %for y = 1:sp.Ny
        %for x = 1:sp.Nx
            %if (y-cy)^2 + (x-cx)^2 < r^2
                %for it = 1:sp.Nt
                    %Hext(3,y,x,it) = amplitude*Ms * sinwave(it);
                %end
            %end
        %end
    %end


%===============================================================================
%% intial condtion for S
%===============================================================================
    random = 1;
    if random
        rng(2012);  % seed the random generator first to reproduce results
        theta = pi .* rand(sp.Ny, sp.Nx);
        phi = 2*pi .* rand(sp.Ny, sp.Nx);
    else
        theta =  45 * pi/180 .* ones(sp.Ny, sp.Nx);
        phi   =  45 * pi/180 .* ones(sp.Ny, sp.Nx);
        %phi = 2*pi .* rand(sp.Ny, sp.Nx);
    end
    Ms = sp.P(1);
    initCond.S = zeros(sp.Ns,sp.Ny,sp.Nx);
    initCond.S(1,:,:) = Ms .* sin(theta) .* cos(phi);
    initCond.S(2,:,:) = Ms .* sin(theta) .* sin(phi);
    initCond.S(3,:,:) = Ms .* cos(theta);
    % assign initCond.S to S(r,t=1) - DON'T FORGET!
    S(:,:,:,1) = initCond.S;


%===============================================================================
%% Boundary conditions for S
    % S_top defines the state of S-vectors just above the row of dots
    %   that is represented by maximum y value, i.e. (y=Ny).
    %   Please note that this is in accordance with xy-Cartesian axes
    %   as opposed to usual Matrix indices where smallest row coordinate
    %   indicates the top-most row. Please take care in defining the
    %   S_top and S_bot boundary conditions in accordance with this Cartesian
    %   notation.
    %
    % Also, the boundary condition vectors are not part of the
    %   simulation domain. They will only interact with the top-most,
    %   bottom-most, right-most and left-most S-vectors in the simulation
    %   domain.
    %
    % Boundary conditions defined only once and outside of time-marching loop
    %   are essentially Dirichlet boundary conditions. To define Neumann and
    %   other types of boundary conditions, modify the sp.boundCond structure
    %   at each iteration within the time-marcing loop.
%===============================================================================
    % Embed boundary conditions as boundCond member in the simParam structure
    % newly defined boundary conditions
    sp.boundCond.S_top = zeros(sp.Ns,sp.Nx);    % +y top
    sp.boundCond.S_bot = zeros(sp.Ns,sp.Nx);    % -y bottom
    sp.boundCond.S_rig = zeros(sp.Ns,sp.Ny);    % +x right
    sp.boundCond.S_lef = zeros(sp.Ns,sp.Ny);    % -x left
    % set some boundary conditions. They can be changed in time-marching
    sp.boundCond.S_lef(3,:) = -5*Ms;    % set left boundary to all down -z
    sp.boundCond.S_rig(3,:) = +5*Ms;    % set right boundary to all up +z
    sp.boundCond.S_top(3,:) = +5*Ms;    % set top boundary to all up +z
    sp.boundCond.S_bot(3,:) = -5*Ms;    % set bottom boundary to all down -z


%===============================================================================
%% Validate the parameters - VERY IMPORTANT!
%===============================================================================
fprintf('INFO: Starting simulation...\n');
fprintf('INFO: Verifying simulation parameters...\n');
[success,sp,S] = validateSimParam(sp,S);
if success == 0
    fprintf('ERROR: Simulation cannot start due to bad parameters.\n');
    return;
else
    fprintf('INFO: All simulation parameters have been verified.\n');
    sp
end


%===============================================================================
%% Time marching
%===============================================================================
fprintf('INFO: Time marching for %d points...\n', sp.Nt);
S(:,:,:,1) = initCond.S;  % assign initCond.S to S(t=1)

for k = 1:length(sp.t)-1   % solve ODE for all time points but last
    fprintf('INFO: %.1f%% t(%d)=%gs: ', 100*single(k)/single(sp.Nt), k, sp.t(k));
    tic;    % code instrumentation
    S(:,:,:,k+1) = odeStepComp(S(:,:,:,k), sp);
    timeTaken = toc;
    fprintf('ODE step executed in %.2f ms runtime\n', timeTaken*1000);

    % Hext time dependence is defined here
    %===============================================================================
    % This defines a 10GHz disc shaped oscillator near the center of the matrix
    r = round(sp.Ny/6);     % radius of disc
    cx = round(sp.Nx/2);    % x-center of disc
    cy = round(sp.Ny/3);    % y-center of disc
    [X Y] = meshgrid(1:sp.Nx, 1:sp.Ny);
    dicsMask = [(Y-cy).^2 + (X-cx).^2 < r^2];
    f = 10e9;
    amplitude = 5;
    sinwave = sin(2*pi*f*sp.t);
    % drive the z-component of Hext only
    sp.Pxy(5,:,:) = dicsMask .* amplitude*Ms * sinwave(k);

    % Flip the boundary conditions for second half of the simulation
    %%===============================================================================
    if k >= int32(sp.Nt/2)
        sp.boundCond.S_lef(3,:) = +5*Ms;
        sp.boundCond.S_rig(3,:) = -5*Ms;
        sp.boundCond.S_top(3,:) = -5*Ms;
        sp.boundCond.S_bot(3,:) = +5*Ms;
    end

    %break
end

%===============================================================================
%% More book-keeping information in simParam structure
%===============================================================================
sp.stopTimeStamp = clock;   % record the stop wall-clock-time
sp.runTime = datevec(datenum(sp.stopTimeStamp) - datenum(sp.startTimeStamp));
sp.finishedWithSuccess = 1; % now mark as a successful finish
fprintf('INFO: simParam structure\n');
sp
fprintf('INFO: Simulation Finished!\n');


%===============================================================================
%% Post processing (optional): save the data and visualize the result
%===============================================================================
%animateDots

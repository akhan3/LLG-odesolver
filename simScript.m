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
        % 2.8177e+09/12/10 to accommodate M(t),M(t+1)Hext,... on the GPU
    sp.Ny = 100;     % #rows of dots in the plane
    sp.Nx = 100;     % #columns of dots in the plane
    sp.dy = 100e-9;     % height of a dot [m]
    sp.dx = 100e-9;     % width of a dot [m]
    % material parameters
    sp.Ms = 8.6e5;          % Saturation Magnetization [A/m]
    sp.gamma = 2.21e5;      % Gyromagnetic Ratio [1/(A/m/s)]
    sp.alpha = 0.05;        % Damping factor [dim-less]
    sp.Aexch = 1.3e-11;     % Exchange constant [J/A]
    sp.Kanis = 1e5;         % Anisotropy constant [J/m^3]
    sp.anisVec = [0 0 1];   % Unit vector defining anisotropy axes [dim-less]
    sp.couplVec = [-.2 -.2 -.2];    % Coupling coefficient [dim-less]
        % It is defined to be negative here and used as it is in field calculation
    sp.demagVec = [.4 .4 .2];       % Demag factor [dim-less]
        % It is defined to be positive here and used as negative in field calculation
    % ODE Solver selection
    sp.useGPU = 0;  % if 1, GPU will be used
    sp.useRK4 = 0;  % if 1, RK4-ODE-solver will be used, otherwise Euler's
    sp.preserveNorm = 1;    % if 1, after the ODE-step is taken,
                            %   all M vectors will be re-normalized to Ms


%===============================================================================
%% Allocate bulk data
    % Dimensions of M and Hext must be treated as follows
    % M(v,r,c,t)    =>  v = 1,2,3 for x,y,z components of the vector
    %                   r = row coordinate
    %                   c = column coordinate in plane
    %                   t = time coordinate
    % It is better to define them as single precision,
    % otherwise they will be converted by the validateSimParam routine anyway.
%===============================================================================
    M = single(zeros(3,sp.Ny,sp.Nx,sp.Nt)); % can be huge in size
    Hext = single(zeros(size(M)));          % same size as of M


%===============================================================================
%% Define External field's time and spatial dependence
    % Hext can either be defined completely beforehand here before entering
    % time-marching loop. It can, of course, also be modified in the
    % time-marching loop
%===============================================================================
    % This defines a 10GHz disc shaped oscillator near the center of the matrix
    r = round(sp.Ny/6);
    cx = round(sp.Nx/2);
    cy = round(sp.Ny/3);
    f = 10e9;
    amplitude = 5;
    sinwave = sin(2*pi*f*sp.t);
    for y = 1:sp.Ny
        for x = 1:sp.Nx
            if (y-cy)^2 + (x-cx)^2 < r^2
                for it = 1:sp.Nt
                    Hext(3,y,x,it) = amplitude*sp.Ms * sinwave(it);
                end
            end
        end
    end
    %for i = 1:sp.Nt
        %Hext(3, round(sp.Ny/3):round(sp.Ny/2), ...
                %round(sp.Nx/3):round(sp.Nx/2), i) = 10*sp.Ms * sinwave(i);
    %end


%===============================================================================
%% intial condtion for M
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
    initCond.M = zeros(3,sp.Ny,sp.Nx);
    initCond.M(1,:,:) = sp.Ms .* sin(theta) .* cos(phi);
    initCond.M(2,:,:) = sp.Ms .* sin(theta) .* sin(phi);
    initCond.M(3,:,:) = sp.Ms .* cos(theta);
    % assign initCond.M to M(r,t=1) - DON'T FORGET!
    M(:,:,:,1) = initCond.M;


%===============================================================================
%% Boundary conditions for M
    % Mtop defines the state of M-vectors just above the row of dots
    %   that is represented by maximum y value, i.e. (y=Ny).
    %   Please note that this is in accordance with xy-Cartesian axes
    %   as opposed to usual Matrix indices where smallest row coordinate
    %   indicates the top-most row. Please take care in defining the
    %   Mtop and Mbot boundary conditions in accordance with this Cartesian
    %   notation.
    %
    % Also, the boundary condition vectors are not part of the
    %   simulation domain. They will only interact with the top-most,
    %   bottom-most, right-most and left-most M-vectors in the simulation
    %   domain.
    %
    % Boundary conditions defined only once and outside of time-marching loop
    %   are essentially Dirichlet boundary conditions. To define Neumann and
    %   other types of boundary conditions, modify the sp.boundCond structure
    %   at each iteration within the time-marcing loop.
%===============================================================================
    % Embed boundary conditions as boundCond member in the simParam structure
    sp.boundCond.Mtop = zeros(3,sp.Nx);    % +y top
    sp.boundCond.Mbot = zeros(3,sp.Nx);    % -y bottom
    sp.boundCond.Mrig = zeros(3,sp.Ny);    % +x right
    sp.boundCond.Mlef = zeros(3,sp.Ny);    % -x left
    % set some boundary conditions. They can be changed in time-marching
    sp.boundCond.Mlef(3,:,:) = -5*sp.Ms;    % set left boundary to all down -z
    sp.boundCond.Mrig(3,:,:) = +5*sp.Ms;    % set right boundary to all up +z
    sp.boundCond.Mtop(3,:,:) = +5*sp.Ms;    % set top boundary to all up +z
    sp.boundCond.Mbot(3,:,:) = -5*sp.Ms;    % set bottom boundary to all down -z


%===============================================================================
%% Validate the parameters - VERY IMPORTANT!
%===============================================================================
fprintf('INFO: Starting simulation...\n');
fprintf('INFO: Verifying simulation parameters...\n');
[success,sp,M,Hext] = validateSimParam(sp,M,Hext);
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
M(:,:,:,1) = initCond.M;  % assign initCond.M to M(t=1)
% H = zeros(3,sp.Ny,sp.Nx); % could pre-allocate Hfield to save time
for k = 1:length(sp.t)-1   % solve ODE for all time points but last
    fprintf('INFO: %.1f%% t(%d)=%gs: ', ...
                100*single(k)/single(sp.Nt), k, sp.t(k));
    tic;    % code instrumentation
    M(:,:,:,k+1) = odeStepComp(sp, M(:,:,:,k), Hext(:,:,:,k));
    timeTaken = toc;
    fprintf('ODE step executed in %.2f ms runtime\n', timeTaken*1000);

    % Flip the boundary conditions for second half of the simulation
    if k >= int32(sp.Nt/2)
        sp.boundCond.Mlef(3,:,:) = +5*sp.Ms;
        sp.boundCond.Mrig(3,:,:) = -5*sp.Ms;
        sp.boundCond.Mtop(3,:,:) = -5*sp.Ms;
        sp.boundCond.Mbot(3,:,:) = +5*sp.Ms;
    end
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

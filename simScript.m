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
    sp.ti = 0;  % initial time
    sp.tf = 1e-9;   % final time
    sp.dt = 1e-13;  % time step
    sp.t = [sp.ti:sp.dt:sp.tf]; % time array
    sp.Nt = length(sp.t)+1;   % number of time points
    % number of the dots
        % Nx*Ny MUST NOT exceed (5000)^2
        % Total GPU memory = 2.8177e+09 Bytes
        % Each dot requires 3*4 Bytes
        % 2.8177e+09/12/10 to accommodate M(t),M(t+1)Hext,... on the GPU
    sp.Ny = 100;     % #rows of dots in the plane
    sp.Nx = 100;     % #columns of dots in the plane
    % material parametrs
        % TODOL: spatially varying
    sp.Ms = 8.6e5;
    sp.gamma = 2.21e5;
    sp.alpha = 0.05;
    sp.cCoupl = [-.2 -.2 -.2];  % cCoupl is defined to be nagative here
                                %   and used as it is in field calculation
    sp.cDemag = [.4 .4 .2];     % cDemag is defined to be positive here
                                %   and used as negative in field calculation
    % ODE Solver selection
    sp.useRK4 = 0;  % if 1, RK4-ODE-solver will be used, otherwise Eueler's
    sp.useGPU = 0;  % if 1, GPU will be used


%===============================================================================
%% Allocate bulk data
%===============================================================================
    % Dimensions of M and Hext must be treated as follows
    % M(v,r,c,t)    =>  v = 1,2,3 for x,y,z components of the vector
    %                   r = row coordinate
    %                   c = column coordinate in plane
    %                   t = time coordinate
    % It is better to define them as single precision,
    % otherwise they will be converted by the validateSimParam routine anyway.
    M = single(zeros(3,sp.Ny,sp.Nx,sp.Nt)); % can be huge in size
    Hext = single(zeros(size(M)));          % same size as of M


%===============================================================================
%% Define External field's time and spatial dependence
%===============================================================================
    Hext(1, round(sp.Ny/3):round(sp.Ny/2), ...
            round(sp.Nx/3):round(sp.Nx/2), :) = 5*sp.Ms;


%===============================================================================
%% intial condtion for M
%===============================================================================
    random = 0;
    if random
        theta = pi .* rand(sp.Ny, sp.Nx);
        phi = 2*pi .* rand(sp.Ny, sp.Nx);
    else
        theta = 40 * pi/180 .* ones(sp.Ny, sp.Nx);
        phi   = 2*pi .* rand(sp.Ny, sp.Nx);
    end
    initCond.M = zeros(3,sp.Ny,sp.Nx);
    initCond.M(1,:,:) = sp.Ms .* sin(theta) .* cos(phi);
    initCond.M(2,:,:) = sp.Ms .* sin(theta) .* sin(phi);
    initCond.M(3,:,:) = sp.Ms .* cos(theta);
    % assign initCond.M to M(r,t=1) - DON'T FORGET!
    M(:,:,:,1) = initCond.M;


%===============================================================================
%% Boundary condtions for M
%===============================================================================
    sp.boundCond.Mtop = zeros(3,sp.Nx);    % +x top
    sp.boundCond.Mbot = zeros(3,sp.Nx);    % -x bottom
    sp.boundCond.Mrig = zeros(3,sp.Ny);    % +y right
    sp.boundCond.Mlef = zeros(3,sp.Ny);    % -y left


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
animateDots

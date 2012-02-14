% clc
clear

addpath src

%% Compile MEX file first
disp 'Compiling...'
cd src
mex -g odeStepComp.c
cd ..


disp 'Running simulation...'

%% Simulation parameters
% Build a simParam structure
sp.simName = 'Untitled Simulation';
sp.tf = 1e-9;
sp.dt = 1e-13;
sp.Ny = 50;   %  Nx*Ny MUST NOT exceed certain value to keep in GPU memory
sp.Nx = 30;
sp.Ms = 8.6e5;
sp.gamma = 2.21e5;
sp.alpha = 0.05;
sp.cCoupl = [-.2 -.2 -.2];
sp.cDemag = [.4 .4 .2];
sp.useRK4 = 0;   % if false, Euler-ODE-solver will be used
sp.useGPU = 0;          % if true, GPU will be used

%% allocate data
t = [0:sp.dt:sp.tf];
sp.Nt = length(t);
M = zeros(3,sp.Ny,sp.Nx,sp.Nt);
Hext = zeros(size(M));

%% External field
Hext(1, sp.Ny/3:sp.Ny/2, sp.Nx/3:sp.Nx/2, :) = 5*sp.Ms;



%% initial condition for M
% start from random
random = 0;
if random
    theta = pi .* rand(sp.Ny, sp.Nx);
    phi = 2*pi .* rand(sp.Ny, sp.Nx);
else
    theta = 40 * pi/180 .* ones(sp.Ny, sp.Nx);
    phi   = 2*pi .* rand(sp.Ny, sp.Nx);
end
ic.M(1,:,:) = sp.Ms .* sin(theta) .* cos(phi);
ic.M(2,:,:) = sp.Ms .* sin(theta) .* sin(phi);
ic.M(3,:,:) = sp.Ms .* cos(theta);


%% Boundary conditions for M
bc.Mtop = zeros(3,sp.Nx);    % +x top
bc.Mbot = zeros(3,sp.Nx);    % -x bottom
bc.Mrig = zeros(3,sp.Ny);    % +y right
bc.Mlef = zeros(3,sp.Ny);    % -y left
            % merge all bc in one large matrix
            % bcM = zeros(3,sp.Ny+2,sp.Nx+2);
            % bcM(:,1,2:end-1) = bc.Mtop;
            % bcM(:,end,2:end-1) = bc.Mbot;
            % bcM(:,2:end-1,1) = bc.Mrig;

% add boundary conditions to the simParam structure
sp.bc = bc;


%% validate the parameters
[success,sp,M,Hext] = validateSimParam(sp,M,Hext);
if success == 0
    fprintf('ERROR: Simulation cannot start due to bad parameters.\n');
    return;
else
    fprintf('INFO: All simulation parameters are verified.\n');
end

% return

%% Time marching
fprintf('Time marching for %d points... ', sp.Nt);
M(:,:,:,1) = ic.M;
% H = zeros(3,sp.Ny,sp.Nx); % could pre-allocate Hfield to save time
for it = 1:length(t)-1      % all but last
    % fprintf('%d\n', it);
    % fprintf('t = %g s\n', t(it));
    M(:,:,:,it+1) = odeStepComp(sp, M(:,:,:,it), Hext(:,:,:,it));
    % timeTaken = toc;
    % fprintf('ODE step executed in %g runtime seconds\n', timeTaken);
end
fprintf('done\n');

%% Animate the result
% animateDots


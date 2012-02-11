clc
clear

%% Messages
disp 'Running simulation...'
mex odeStepComp.c

%% parameters
sp.tf = 5e-13;
sp.dt = 1e-13;
sp.Ny = 4;   %  Nx*Ny MUST NOT exceed (3,000)^2 = 9,000,000
sp.Nx = 5;
sp.Ms = 8.6e5;
sp.gamma = 2.21e5;
sp.alpha = 0.05;
sp.cCoupl = -0.2;
sp.cDemag = diag([.4 .4 .2]);
sp.odeSolver_rk4 = 0;   % if false, Euler-ODE-solver will be used
sp.useGPU = 0;          % if true, GPU will be used

%% allocate data
t = [0:sp.dt:sp.tf];
sp.Nt = length(t);
M = zeros(3,sp.Ny,sp.Nx,sp.Nt);
Hext = zeros(size(M));
% if sp.useGPU
%     M = gpuArray(M);
%     Hext = gpuArray(Hext);
% end

%% initial condition for M
% start from random
theta = pi.*rand(sp.Ny, sp.Nx);
phi = 2*pi.*rand(sp.Ny, sp.Nx);
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

%% Time marching
M(:,:,:,1) = ic.M;
% H = zeros(3,sp.Ny,sp.Nx); % could pre-allocate Hfield to save time
for it = 1:length(t)-1      % all but last
    A = odeStep(sp, bc, M(:,:,:,it), Hext(:,:,:,it));
    % M(:,:,:,it+1) = odeStep(sp, bc, M(:,:,:,it), Hext(:,:,:,it));
    % timeTaken = toc;
    % fprintf('ODE step executed in %g runtime seconds\n', timeTaken);
end

%% Animate
% visualize_dots

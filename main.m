% clc
clear

%% parameters
sp.tf = 5e-13;
sp.dt = 5e-13;
sp.Nx = 2000;
sp.Ny = 2000;
sp.Ms = 8.6e5;
sp.gamma = 2.21e5;
sp.alpha = 0.05;
sp.c = -0.2;
sp.N = diag([.4 .4 .2]);
sp.odeSolver_rk4 = false;    % if false, Euler-ODE-solver will be used
sp.useGPU = 1;          % if true, GPU will be used

%% allocate data
t = [0:sp.dt:sp.tf];
sp.Nt = length(t);
M = zeros(3,sp.Ny,sp.Nx,sp.Nt);
Hext = zeros(size(M));
if sp.useGPU
    M = gpuArray(M);
    Hext = gpuArray(Hext);
end

%% initial condition for M
% start from random
theta = pi.*rand(sp.Ny, sp.Nx);
phi = 2*pi.*rand(sp.Ny, sp.Nx);
ic.M(1,:,:) = sp.Ms * sin(theta)*cos(phi);
ic.M(2,:,:) = sp.Ms * sin(theta)*sin(phi);
ic.M(3,:,:) = sp.Ms * cos(theta);

%% Boundary conditions for M
bc.Mtop = zeros(3,sp.Nx);    % +x top
bc.Mbot = zeros(3,sp.Nx);    % -x bottom
bc.Mrig = zeros(3,sp.Ny);    % +y right
bc.Mlef = zeros(3,sp.Ny);    % -y left

%% Time marching
% H = zeros(3,sp.Ny,sp.Nx);   % could pre-allocate Hfield to save time
%     H = zeros(3,sp.Ny,sp.Nx);   % but then have to reset H before each
%     ODE step

M(:,:,:,1) = ic.M;
for it = 1:length(t)-1      % all but last
    tic;
    M(:,:,:,it+1) = odeStep(sp, bc, M(:,:,:,it), Hext(:,:,:,it));
    timeTaken = toc;
    fprintf('ODE step executed in %g runtime seconds\n', timeTaken);
end

%% Animate
% visualize_dots


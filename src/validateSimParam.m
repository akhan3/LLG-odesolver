function [success,sp,S] = validateSimParam(sp,S)

    success = 1;

%===============================================================================
%% Default sp (simParam structure)
%===============================================================================
    % Book-keeping information
    spDefault.simName = 'Untitled Simulation';
    spDefault.startTimeStamp = clock;
    spDefault.stopTimeStamp = inf * clock;
    spDefault.runTime = inf * clock;
    spDefault.finishedWithSuccess = int32(0);
    % Simulation time
    spDefault.ti = single(0);
    spDefault.tf = single(1e-11);
    spDefault.dt = single(1e-13);
    spDefault.t = single([spDefault.ti:spDefault.dt:spDefault.tf]);
    spDefault.Nt = int32(length(spDefault.t));
    % number and size of the dots
    spDefault.Ny = int32(20);
    spDefault.Nx = int32(20);
    spDefault.dy = single(100e-9);
    spDefault.dx = single(100e-9);
    % material parametrs
    spDefault.Ns = int32(1);                        % numer of state variables at each mesh point
    spDefault.P = single(zeros(1,1));               % Paramters that are uniform over the mesh
    spDefault.Pxy = single(zeros(1, sp.Ny,sp.Nx));  % Paramters that vary over the mesh
    spDefault.Np = length(sp.P);         % numer of spatially uniform paramters
    spDefault.Npxy = size(sp.Pxy, 1);    % numer of spatially varying paramters
    % ODE Solver selection
    spDefault.useRK4 = int32(0);
    spDefault.useGPU = int32(0);
    spDefault.useGPUnum = int32(0);
    % Boundary condition defaults to ZERO
    spDefault.boundCond.S_top = single(zeros(spDefault.Ns,sp.Nx));    % +x top
    spDefault.boundCond.S_bot = single(zeros(spDefault.Ns,sp.Nx));    % -x bottom
    spDefault.boundCond.S_rig = single(zeros(spDefault.Ns,sp.Ny));    % +y right
    spDefault.boundCond.S_lef = single(zeros(spDefault.Ns,sp.Ny));    % -y left


%===============================================================================
%% Valiade sp (simParam structure)
%===============================================================================
    fnames = fieldnames(spDefault);
    for i = 1:length(fnames)
        if ~fieldExists(sp, fnames(i))
            fprintf('INFO: sp.%s field does not exist, so creating with default value\n', char(fnames(i)));
            cmd = sprintf('sp.%s = spDefault.%s;', char(fnames(i)), char(fnames(i)));
            eval(cmd);
        end
    end

    % sanity check for sp.t and sp.Nt
    if (sum(sp.t ~= [sp.ti:sp.dt:sp.tf]) ~= 0)  % if any element in time array is bad
        fprintf('INFO: sp.t is badly formed, so fixed it!\n');
        sp.t = single([sp.ti:sp.dt:sp.tf]);
    end
    if (sp.Nt ~= length(sp.t))
        fprintf('INFO: sp.Nt does not match the length of time array, so fixed it!\n');
        sp.Nt = int32(length(sp.t));
    end

    % sanity check for boundary conditions
    if( size(sp.boundCond.S_top, 2) ~= sp.Nx ) || ...
            ( size(sp.boundCond.S_bot, 2) ~= sp.Nx ) || ...
            ( size(sp.boundCond.S_rig, 2) ~= sp.Ny ) || ...
            ( size(sp.boundCond.S_lef, 2) ~= sp.Ny )
        fprintf('ERROR: Size of sp.boundCond does not match that of the number of dots!\n');
        success = 0;
        return;
    end


%===============================================================================
%% convert everything to single precision or 32-bit integer
%===============================================================================
    sp.ti = single(sp.ti);
    sp.tf = single(sp.tf);
    sp.dt = single(sp.dt);
    sp.t = single(sp.t);
    sp.Nt = int32(sp.Nt);
    sp.Ny = int32(sp.Ny);
    sp.Nx = int32(sp.Nx);
    sp.dy = single(sp.dy);
    sp.dx = single(sp.dx);
    sp.Ns = int32(sp.Ns);
    sp.Np = int32(sp.Np);
    sp.Npxy = int32(sp.Npxy);
    sp.P = single(sp.P);
    sp.Pxy = single(sp.Pxy);
    sp.useRK4 = int32(sp.useRK4);
    sp.useGPU = int32(sp.useGPU);
    sp.useGPUnum = int32(sp.useGPUnum);
    sp.boundCond.S_top = single(sp.boundCond.S_top);
    sp.boundCond.S_bot = single(sp.boundCond.S_bot);
    sp.boundCond.S_rig = single(sp.boundCond.S_rig);
    sp.boundCond.S_lef = single(sp.boundCond.S_lef);
    % S(r,t) is a bulky array. convert only if not already single
    if strcmp(class(S), 'single') == 0
        fprintf('INFO: converting S to single precision... ');
        S = single(S);
        fprintf('done\n');
    end
end


%===============================================================================
%% check if a field exists in a structure
%===============================================================================
function present = fieldExists(struct, fname)
    present = any(strcmp(fname,fieldnames(struct)));
end

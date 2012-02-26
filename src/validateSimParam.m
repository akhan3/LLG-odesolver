function [success,sp,M,Hext] = validateSimParam(sp,M,Hext)

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
    spDefault.Ms = single(8.6e5);
    spDefault.gamma = single(2.21e5);
    spDefault.alpha = single(0.05);
    spDefault.Aexch = single(1.3e-11);
    spDefault.cCoupl = single([-.2 -.2 -.2]);
    spDefault.cDemag = single([.4 .4 .2]);
    % ODE Solver selection
    spDefault.useGPU = int32(0);
    spDefault.useRK4 = int32(0);
    spDefault.preserveNorm = int32(1);
    % Boundary condition defaults to ZERO
    spDefault.boundCond.Mtop = single(zeros(3,sp.Nx));    % +x top
    spDefault.boundCond.Mbot = single(zeros(3,sp.Nx));    % -x bottom
    spDefault.boundCond.Mrig = single(zeros(3,sp.Ny));    % +y right
    spDefault.boundCond.Mlef = single(zeros(3,sp.Ny));    % -y left


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

    % sanity check for sp.Aexch, sp.cCoupl and sp.cDemag
    if (length(sp.Aexch) ~= 1)
        fprintf('ERROR: sp.Aexch must be a scalar!\n');
        success = 0;
        return;
    elseif (length(sp.cCoupl) ~= 3)
        fprintf('ERROR: sp.cCoupl must contain only 3 elements!\n');
        success = 0;
        return;
    elseif (length(sp.cDemag) ~= 3)
        fprintf('ERROR: sp.cDemag must contain only 3 elements!\n');
        success = 0;
        return;
    end

    % sanity check for boundary conditions
    if( length(sp.boundCond.Mtop) ~= sp.Nx ) || ...
            ( length(sp.boundCond.Mbot) ~= sp.Nx ) || ...
            ( length(sp.boundCond.Mrig) ~= sp.Ny ) || ...
            ( length(sp.boundCond.Mlef) ~= sp.Ny )
        fprintf('ERROR: Size of sp.boundCond does not match that of the number of dots!\n');
        success = 0;
        return;
    end

    %if ~fieldExists(sp,'Ms')
        %sp.Ms = spDefault.Ms;
    %elseif sp.Ms == 0
        %fprintf('Ms cannot have a value of %g [A/m]\n', sp.Ms);
        %success = 0;
    %end


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
    sp.Ms = single(sp.Ms);
    sp.gamma = single(sp.gamma);
    sp.alpha = single(sp.alpha);
    sp.Aexch = single(sp.Aexch);
    sp.cCoupl = single(sp.cCoupl);
    sp.cDemag = single(sp.cDemag);
    sp.useRK4 = int32(sp.useRK4);
    sp.useGPU = int32(sp.useGPU);
    sp.boundCond.Mtop = single(sp.boundCond.Mtop);
    sp.boundCond.Mbot = single(sp.boundCond.Mbot);
    sp.boundCond.Mrig = single(sp.boundCond.Mrig);
    sp.boundCond.Mlef = single(sp.boundCond.Mlef);
    % M(r,t) and Hext(r,t) bulky arrays
    if strcmp(class(M), 'single') == 0
        fprintf('INFO: converting M to single precision... ');
        M = single(M);
        fprintf('done\n');
    end
    if strcmp(class(Hext), 'single') == 0
        fprintf('INFO: converting Hext to single precision... ');
        Hext = single(Hext);
        fprintf('done\n');
    end
end


%===============================================================================
%% check if a field exists in a structure
%===============================================================================
function present = fieldExists(struct, fname)
    present = any(strcmp(fname,fieldnames(struct)));
end

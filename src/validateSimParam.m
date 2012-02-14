function [success,sp,M,Hext] = validateSimParam(sp,M,Hext)

    success = 1;

    %% Default sp (simParam structure)
    spDefault.simName = 'Untitled Simulation';
    spDefault.tf = single(1e-11);
    spDefault.dt = single(1e-13);
    spDefault.Nt = int32(floor(spDefault.tf/spDefault.dt) + 1);
    spDefault.Ny = int32(20);
    spDefault.Nx = int32(20);
    spDefault.Ms = single(8.6e5);
    spDefault.gamma = single(2.21e5);
    spDefault.alpha = single(0.05);
    spDefault.cCoupl = single([-.2 -.2 -.2]);
    spDefault.cDemag = single([.4 .4 .2]);
    spDefault.useRK4 = int32(0);
    spDefault.useGPU = int32(0);


    %% valiade sp (simParam structure)
    if ~fieldExists(sp,'simName')
        sp.simName = spDefault.simName;
    end

    if ~fieldExists(sp,'tf')
        sp.tf = spDefault.tf;
    end

    if ~fieldExists(sp,'dt')
        sp.dt = spDefault.dt;
    end

    if ~fieldExists(sp,'Nt') || ~(sp.Nt == floor(sp.tf/sp.dt) + 1)
        sp.Nt = int32(floor(sp.tf/sp.dt) + 1);
    end

    if ~fieldExists(sp,'Nx')
        sp.Nx = spDefault.Nx;
    end

    if ~fieldExists(sp,'Ny')
        sp.Ny = spDefault.Ny;
    end

    if ~fieldExists(sp,'Ms')
        sp.Ms = spDefault.Ms;
    elseif sp.Ms == 0
        fprintf('Ms cannot have a value of %g [A/m]\n', sp.Ms);
        success = 0;
    end

    if ~fieldExists(sp,'gamma')
        sp.gamma = spDefault.gamm;
    elseif sp.gamma == 0
        fprintf('gamma cannot have a value of %g\n', sp.gamma);
        success = 0;
    end

    if ~fieldExists(sp,'alpha')
        sp.alpha = spDefault.alpha;
    end

    if ~fieldExists(sp,'cCoupl')
        sp.cCoupl = spDefault.cCoupl;
    end

    if ~fieldExists(sp,'cDemag')
        sp.cDemag = spDefault.cDemag;
    end

    if ~fieldExists(sp,'useRK4')
        sp.useRK4 = spDefault.useRK4;
    end

    if ~fieldExists(sp,'useGPU')
        sp.useGPU = spDefault.useGPU;
    end

    %% convert everything to single precision or 32-bit integer
    sp.tf = single(sp.tf);
    sp.dt = single(sp.dt);
    sp.Nt = int32(sp.Nt);
    sp.Nx = int32(sp.Nx);
    sp.Ny = int32(sp.Ny);
    sp.Ms = single(sp.Ms);
    sp.gamma = single(sp.gamma);
    sp.alpha = single(sp.alpha);
    sp.cCoupl = single(sp.cCoupl);
    sp.cDemag = single(sp.cDemag);
    sp.useRK4 = int32(sp.useRK4);
    sp.useGPU = int32(sp.useGPU);
    % boundary conditions
    sp.bc.Mtop = single(sp.bc.Mtop);
    sp.bc.Mbot = single(sp.bc.Mbot);
    sp.bc.Mrig = single(sp.bc.Mrig);
    sp.bc.Mlef = single(sp.bc.Mlef);
    % M(r,t) and Hext(r,t) bulky arrays
    M = single(M);
    Hext = single(Hext);
end



%% check if a field exists in a structure
function present = fieldExists(struct, fname)
    present = any(strcmp(fname,fieldnames(struct)));
end

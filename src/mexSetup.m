function mexSetup(mexDir)

%% Directory setup and MEX file compilation
SRC_DIR = fullfile(mexDir);
MEX_OPT = fullfile(SRC_DIR, 'mexopts.sh');
MEX_CFILE = fullfile(SRC_DIR, 'odeStepComp.c');
addpath(SRC_DIR);

compileCmd = ['mex -f ',MEX_OPT,' -outdir ',SRC_DIR,' ',MEX_CFILE];
fprintf('INFO: Compiling MEX file...\n%s\n', compileCmd);
eval(compileCmd);
fprintf('INFO: MEX filed compiled successfully.\n');

end

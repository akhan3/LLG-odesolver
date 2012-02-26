function mexSetup(mexDir)

%% Directory setup and MEX file compilation
SRC_DIR = fullfile(mexDir);
MEX_CFILE = fullfile(SRC_DIR, 'odeStepComp.c');
addpath(SRC_DIR);
MEX_OPT = 'CFLAGS="\$CFLAGS -std=c99 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"';

compileCmd = ['mex ', MEX_OPT, ' -outdir ',SRC_DIR,' ',MEX_CFILE];
fprintf('INFO: Compiling MEX file...\n%s\n', compileCmd);
eval(compileCmd);
fprintf('INFO: MEX file compiled successfully.\n');

end

function [status] = mexSetup(mexDir)

%% Source files
SRC_DIR = fullfile(mexDir);
MEX_FILE = fullfile(SRC_DIR, 'odeStepComp.cpp');
CUDA_FILE = fullfile(SRC_DIR, 'kernel.cu');
CUDA_OBJFILE = fullfile(SRC_DIR, 'kernel.o');

%% Includes and libraries
CUDA_INC = sprintf('-I%s/include -I%s/C/common/inc', getenv('CUDA_ROOT'), getenv('CUDA_SDK'));
CUDA_LIB = sprintf('-L%s/lib64 -lcuda -lcudart  -L%s/C/lib/ -lcutil_x86_64', getenv('CUDA_ROOT'), getenv('CUDA_SDK'));
MAT_INC = sprintf('-I%s/extern/include', matlabroot);

% Convert cu file to cpp
cudaCmd = sprintf('nvcc -c -g %s --compiler-options "-fPIC -fopenmp -fno-strict-aliasing" %s -o %s', CUDA_INC, CUDA_FILE, CUDA_OBJFILE);
fprintf('INFO: Converting CUDA file...\n%s\n\n', cudaCmd);
status = system(cudaCmd);

if(status ~= 0)
    fprintf('ERROR: Cannot convert CUDA file %s due to errors\n', CUDA_FILE);
    return
end

% Compile and link with mex
mexCmd = sprintf('mex -g %s %s CFLAGS="\\$CFLAGS -fopenmp -fPIC" CXXFLAGS="\\$CXXFLAGS -fopenmp -fPIC" LDFLAGS="\\$LDFLAGS -fopenmp -fPIC" -outdir %s %s', MEX_FILE, CUDA_OBJFILE, SRC_DIR, CUDA_LIB);
fprintf('INFO: Compiling MEX file...\n%s\n', mexCmd);
eval(mexCmd);
fprintf('INFO: MEX file compiled successfully.\n');


end

%% SETUP EXTERNAL LIBRARIES
% AWFUL
try
  au_sparse(int32([1 2 3 4]), int32([1 3 5 7]), randn(4,1));
catch
  disp('Setting up AWFUL...');
  cd External/awful/matlab
  addpath(genpath(pwd));
  au_mexall
  cd ../../..
end

% Manopt
if exist('manopt_version', 'file') ~= 2
  disp('Setting up Manopt...');
  cd 'External/manopt';
  importmanopt;
  cd ../..
end

%% ADD PATHS
disp('Adding paths...')
cd Experiment
addpath(genpath(pwd));
cd ..

cd Algorithms
addpath(genpath(pwd));
cd JE
rmpath([pwd, '/ce_src/build']);
rmpath([pwd, '/ce_src/LRMC']);
rmpath([pwd, '/ce_src']);
cd ..

%% COMPILE MEX FILES
% JE
cd JE
disp('Compiling mex files for JE');
mex je_form_hessian.cpp
mex je_add_hessian_layer.cpp
cd ..

% NB
cd NB
disp('Compiling mex files for NB');
mex -largeArrayDims xt_nb_spmaskmult.c
mex -largeArrayDims xt_nb_setsparseentries.c
mex -largeArrayDims -lmwlapack xt_nb_cholsolvecell.c
mex -largeArrayDims -lmwlapack -lmwblas xt_nb_spbuildmatrix.c
mex -largeArrayDims -lmwlapack -lmwblas xt_nb_buildAchol.c
cd ..

% PG
cd PG
disp('Compiling mex files for PG');
mex xt_pg_kronmex.c
cd ..

cd ..

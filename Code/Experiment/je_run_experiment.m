function je_run_experiment(dataset, r, varargin)
% JE_RUN_EXPERIMENT

%% OPTIONS
% Extend varargin to include tolerance, maximum number of iterations and
% display.
warning('off', 'MATLAB:nearlySingularMatrix'); % turn off warning for better visual
varargin = [
  {
  'tol=1e-10', ...
  'max_iter=300', ...
  'display=0', ...
  'nu=0.1', ...
  'max_als=10', ...
  'output_filename=0', ...
  }, ...
  varargin, ...
  ];

% Convert varargin to opts.
opts = au_opts(...
  'sample_start=1;sample_end=0;num_samples=100', ...
  'overwrite=0', ...
  varargin{:});

% Make the first letter of the dataset name uppercase.
dataset = regexprep(lower(dataset),'(\<[a-z])','${upper($1)}');
if ~opts.output_filename,
  file = ['Results/', dataset, '/', lower(dataset), '_r', int2str(r)];
else
  file = ['Results/', dataset, '/', lower(opts.output_filename), '_r', int2str(r)];
end

% If a previous version of the file exists, load it.
clearvars -global lrmc_results
if exist([file, '.mat'], 'file') == 2, load(file);
end

% Set lrmc_results as a global variable.
global lrmc_results

%% LIST OF ALGORITHMS
if ~opts.sample_end
  opts.sample_end = opts.sample_start + opts.num_samples - 1;
end

data = load(['Datasets/' dataset '/' lower(dataset)]); % Load the measurement matrix.
if ~issparse(data.M), data.M(data.W==0) = 100; % So we get an alert if some algo accesses zeros in M.
end

% Write the measurement matrix and the weight matrix in binary.
sample_id = lower(dataset);
filename = ['Temp/', sample_id, '_r', num2str(r)];

if ~issparse(data.M)
  je_write_binary_matrix([filename, '_M.bin'], data.M);
  je_write_binary_matrix([filename, '_W.bin'], data.W);
end

algorithms = {
  % Latest custom-coded algorithms
  'ALS' @(sample_num) je_alternation(data.M, data.W, r, sample_num, nan, varargin{:}, 'retraction=0')
  'PF' @(sample_num) je_alternation(data.M, data.W, r, sample_num, nan, varargin{:}, 'retraction=1')
  'DW' @(sample_num) je_wiberg(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', 'retraction=0')
  'DRW1' @(sample_num) je_wiberg(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', 'alg_type=1', 'constraint=0')
  'DRW1_UWOFF' @(sample_num) je_wiberg(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', 'alg_type=1', 'constraint=0', 'search_unique_weights=0')
  'DRW1P' @(sample_num) je_wiberg(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', 'alg_type=1')
  'DRW1P_UWOFF' @(sample_num) je_wiberg(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', 'alg_type=1', 'search_unique_weights=0')
  'DRW2' @(sample_num) je_wiberg(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', 'constraint=0')
  'DRW2_UWOFF' @(sample_num) je_wiberg(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', 'constraint=0', 'search_unique_weights=0')
  'DRW2P' @(sample_num) je_wiberg(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0')
  'DRW2P_UWOFF' @(sample_num) je_wiberg(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', 'search_unique_weights=0')
  ...
  % Latest Ceres-solver LMs
  'CE_LM' @(sample_num) je_ceres_lrmc_wrapper(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', 'max_als=0', 'write_binary_mw=0', ['sample_id=', sample_id], 'retraction=0');
  'CE_LMI' @(sample_num) je_ceres_lrmc_wrapper(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', 'max_als=0', 'use_inner_iters=1', 'write_binary_mw=0', ['sample_id=', sample_id], 'retraction=0');
  'CE_ALM' @(sample_num) je_ceres_lrmc_wrapper(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', sprintf('max_als=%d', opts.max_als), 'write_binary_mw=0', ['sample_id=', sample_id], 'retraction=0');
  'CE_ALMI' @(sample_num) je_ceres_lrmc_wrapper(data.M, data.W, r, sample_num, nan, varargin{:}, 'nu=0.0', sprintf('max_als=%d', opts.max_als), 'use_inner_iters=1', 'write_binary_mw=0', ['sample_id=', sample_id], 'retraction=0');
  'CE_ARULM' @(sample_num) je_ceres_lrmc_wrapper(data.M, data.W, r, sample_num, nan, varargin{:}, sprintf('nu=%d', opts.nu), sprintf('max_als=%d', opts.max_als), 'reg_unreg=1', 'write_binary_mw=0', ['sample_id=', sample_id], 'retraction=0');
  'CE_ARULMI' @(sample_num) je_ceres_lrmc_wrapper(data.M, data.W, r, sample_num, nan, varargin{:}, sprintf('nu=%d', opts.nu), sprintf('max_als=%d', opts.max_als), 'reg_unreg=1', 'use_inner_iters=1', 'write_binary_mw=0', ['sample_id=', sample_id], 'retraction=0');
  ...
  % External algorithms
  'PG_CSF' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=PG_CSF');
  'TO_DW' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=TO_DW');
  'CH_LM_S' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=CH_LM_S');
  'CH_LM_S_GN' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=CH_LM_S_GN');
  'CH_LM_S_RW2' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=CH_LM_S_RW2');
  'CH_LM_M' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=CH_LM_M');
  'CH_LM_M_GN' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=CH_LM_M_GN');
  'CH_LM_M_RW2' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=CH_LM_M_RW2');
  'CO_LM_S' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=CO_LM_S');
  'CO_LM_M' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=CO_LM_M');
  'CO_LM_M_GN' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=CO_LM_M_GN');
  'NB_RTRMC' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=NB_RTRMC');
  'RC_ALM' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=RC_ALM');
  'RC_RCALM' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=RC_RCALM', sprintf('r_init=%d', r + 2));
  'DB_BALM' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=DB_BALM');
  'DO_BALM' @(sample_num) je_external_wrapper(data.M, data.W, r, sample_num, varargin{:}, 'alg=DO_BALM');
  };

% Results now stored in struct array, so that e.g. time for LM_EU
% are in lrmc_results.(dataset)(sample).LM_EU.time
% To get an array of results from one algorithm, use
% lm_eu = [lrmc_results.(dataset).LM_EU];
% [lm_eu.runtime; lm_eu.cost; lm_eu.iters]'

%% EXPERIMENT
for sample_num = opts.sample_start : opts.sample_end
  
  fprintf('Experiments running on [%s - rank %d] sample [%d] \n', dataset, r, sample_num);
  lrmc_results.(dataset)(sample_num).sample_num = sample_num;
  lrmc_results.(dataset)(sample_num).rank = r;
  
  for k = 1 : length(algorithms)
    
    % If the input files are accidentally deleted, write them again.
    if ~issparse(data.M) && (exist([file, '.mat'], 'file') ~= 2)
      je_write_binary_matrix([filename, '_M.bin'], data.M);
      je_write_binary_matrix([filename, '_W.bin'], data.W);
    end
    
    alg = algorithms{k,1};
    alg_caller = algorithms{k,2};
    if isfield(opts, (alg))
      if opts.(alg)
        if ~opts.overwrite && ...
            isfield(lrmc_results.(dataset)(sample_num), alg) && ...
            isfield(lrmc_results.(dataset)(sample_num).(alg), 'cost')
          fprintf('Algorithm [%s] already run on [%s - rank %d] sample [%d]: %.6f\n', ...
            alg, dataset, r, sample_num, lrmc_results.(dataset)(sample_num).(alg).cost);
        else
          fprintf('Algorithm [%s] running on [%s - rank %d] sample [%d]', alg, dataset, r, sample_num);
          time_start = tic;
          [U, V, err, iter] = alg_caller(sample_num);
          runtime = toc(time_start);
          
          lrmc_results.(dataset)(sample_num).(alg).iters = iter;
          lrmc_results.(dataset)(sample_num).(alg).cost = err;
          lrmc_results.(dataset)(sample_num).(alg).runtime = runtime;
          lrmc_results.(dataset)(sample_num).(alg).U = U;
          lrmc_results.(dataset)(sample_num).(alg).V = V;
          
          fprintf(': %.6f, %03d iters, %.3f iters per sec\n', ...
            lrmc_results.(dataset)(sample_num).(alg).cost, ...
            lrmc_results.(dataset)(sample_num).(alg).iters, ...
            lrmc_results.(dataset)(sample_num).(alg).iters / ...
            lrmc_results.(dataset)(sample_num).(alg).runtime)
          
          if ~exist('Results', 'dir'), mkdir('Results');
          end
          
          if ~exist(['Results/', dataset], 'dir'), mkdir(['Results/', dataset]);
          end
          
          save(file, 'lrmc_results');
        end
      end
    end
  end
end

if ~issparse(data.M)
  delete([filename, '_M.bin']);
  delete([filename, '_W.bin']);
end

clearvars -global lrmc_results

end

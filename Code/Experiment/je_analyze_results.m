function je_analyze_results(dataset, r, varargin)
% JE_ANALYZE_RESULTS

% These flags say which methods to run
%% OPTIONS
opts = au_opts(...
  'sample_start=1;sample_end=0;num_samples=100', ...
  'tol=1e-6', ...
  'csv=0', ...
  'display=1', ...
  'output_filename=0', ...
  'folder_suffix=0', ...
  'display_reconstructions=0', ...
  'save_reconstructions=0', ...
  varargin{:});

alg_names = {
  % Latest custom-coded algorithms
  'ALS', 'New'         % Alternation (Undamped RW3)
  'DW', 'New'          % Damped Wiberg (Damped RW2 + projection constraint)
  'DRW1_UWOFF', '(NO UW)'
  'DRW1', '(UW)'
  'DRW1P_UWOFF', '(NO UW)'
  'DRW1P', '(UW)'   % Damped RW1 with projection constraint (Damped RW1 + projection constraint + retraction)
  'DRW2_UWOFF', '(NO UW)'
  'DRW2', '(UW)'
  'DRW2P_UWOFF', '(NO UW)'
  'DRW2P', '(UW)'   % Damped RW2 with projection constraint (Damped RW2 + projection constraint + retraction)
  % External algorithms
  'PG_CSF', 'Orig.'      % Gotardo's CSF (Damped RW2)
  'TO_DW', 'Orig.'       % Okatani's Damped Wiberg (Damped RW2 + projection constraint)
  'CH_LM_S', 'Mod.'      % Chen's LM_S (Damped full Newton)
  'CH_LM_S_GN', 'Mod.'   % Chen's LM_S modified to use GN-matrix (Damped RW1)
  'CH_LM_S_RW2', 'Mod.'  % Chen's LM_S modified to use approximate GN-matrix (Damped RW2)
  'CH_LM_M', 'Mod.'      % Chen's LM_M (Reduced damped full Newton)
  'CH_LM_M_GN', 'Mod.'   % Chen's LM_M_GN (Reduced damped RW1)
  'CH_LM_M_RW2', 'Mod.'  % Chen's LM_M modified to use approximate GN-matrix (Reduced damped RW2)
  'CO_LM_S', 'Orig.'     % Chen's original implementation of LM_S
  'CO_LM_M', 'Orig.'     % Chen's original implementation of LM_M
  'CO_LM_M_GN', 'Orig.'  % Chen's original implementation of LM_M_GN
  'NB_RTRMC', 'Orig.'    % Boumal's original implementation of RTRMC
  'RC_ALM', 'Unopt.'     % Cabral's ALM (New but unoptimized implementation)
  'RC_RCALM', 'Unopt.'   % Cabral's ALM with rank continuation (New but unoptimized implementation)
  'DB_BALM', 'Unopt.'    % Del Bue's BALM (New but unoptimized implementation)
  'DO_BALM', 'Al. Orig.' % Del Bue's "almost-original" implementation of BALM with few tweaks
  % Latest Ceres-solver LMs
  'CE_LM', 'Ceres'       % CERES LM
  'CE_LMI', 'Ceres'      % CERES LM INNER
  'CE_ALM', 'Ceres'      % 10 ALS + CE_LM
  'CE_ALMI', 'Ceres'     % 10 ALS + CE_LM_I
  'CE_ARULM', 'Ceres'    % 10 reg. ALS + reg. CE_LM + unreg. CE_LM
  'CE_ARULMI', 'Ceres'   % 10 reg. ALS + reg. CE_LMI + unreg. CE_LMI
  };

% Make the first letter of the dataset name uppercase.
dataset = regexprep(lower(dataset),'(\<[a-z])','${upper($1)}');

% Add folder suffix
file = 'Results';
if opts.folder_suffix ~= 0, file = [file, '_', opts.folder_suffix];
end

if ~opts.output_filename,
  file = [file, '/', dataset, '/', lower(dataset), '_r', int2str(r)];
else
  file = [file, '/', dataset, '/', lower(opts.output_filename), '_r', int2str(r)];
end

% If a previous version of the file exists, load it.
if exist([file, '.mat'], 'file') == 2, load(file);
else
  error('The specified results file does not exist.');
end

best_min = nan;
N = size(alg_names, 1);

% Create a list of used algorithms and search for the best optimum.
j = 1;
alg_list.idx = zeros(N, 1);
for i=1:N
  alg_name = alg_names{i,1};
  if isfield(opts, (alg_name))
    if opts.(alg_name)
      alg_list.idx(j) = i;
      j = j + 1;
      
      % If the specified algorithm does not exist, continue with the rest
      % of the for loop.
      if ~isfield([lrmf_results.(dataset)], (alg_name)), continue
      end
      
      if opts.sample_end == 0,
        alg = [lrmf_results.(dataset)(opts.sample_start:end).(alg_name)];
      else
        alg = [lrmf_results.(dataset)(opts.sample_start:opts.sample_end).(alg_name)];
      end
      if (isnan(best_min)) || (best_min > min([alg.cost])), best_min = min([alg.cost]);
      end
    end
  end
end
alg_list.idx = alg_list.idx(alg_list.idx ~= 0);
N = length(alg_list.idx);

% Initialize arrays list.algorithms = algorithms{alg_idx,1};

alg_list.runs.all = nan(N,1);
alg_list.runs.conv = nan(N,1);
alg_list.iters.all.mean = nan(N,1);
alg_list.iters.all.median = nan(N,1);
alg_list.iters.success.mean = nan(N,1);
alg_list.iters.success.median = nan(N,1);
alg_list.iters.fail.mean = nan(N,1);
alg_list.iters.fail.median = nan(N,1);
alg_list.runtime.all.mean = nan(N,1);
alg_list.runtime.all.median = nan(N,1);
alg_list.runtime.success.mean = nan(N,1);
alg_list.runtime.success.median = nan(N,1);
alg_list.runtime.fail.mean = nan(N,1);
alg_list.runtime.fail.median = nan(N,1);
alg_list.ips.all.mean = nan(N,1);
alg_list.ips.all.median = nan(N,1);
alg_list.RUSSO.conv = nan(N,1);
alg_list.RUSSO.all = nan(N,1);
alg_list.RUSSO.optima = cell(N, 1);
alg_list.RUSSO.runtime = cell(N,1);
alg_list.triple_RUSSO.conv = nan(N, 1);
alg_list.triple_RUSSO.all = nan(N, 1);
alg_list.triple_RUSSO.optima = cell(N, 1);
alg_list.triple_RUSSO.runtime = cell(N,1);
alg_list.TSS.mean = nan(N,1);
alg_list.TSS.std = nan(N,1);
alg_list.TSS.mean_est = nan(N,1);
alg_list.TSS.std_est = nan(N,1);
alg_list.TRUSSO.mean = nan(N, 1);
alg_list.TRUSSO.std = nan(N, 1);

for i=1:N
  alg_name = alg_names{alg_list.idx(i),1};
  
  % If the specified algorithm exists, do the calculations. Otherwise,
  % continue with the rest of the for loop.
  if opts.(alg_name) && isfield([lrmf_results.(dataset)], (alg_name))
    if opts.sample_end == 0
      alg = [lrmf_results.(dataset)(opts.sample_start:end).(alg_name)];
    else
      alg = [lrmf_results.(dataset)(opts.sample_start:opts.sample_end).(alg_name)];
    end
    alg_conv = [alg.cost] - best_min < opts.tol  *best_min;
    alg_all_iters = [alg.iters];
    alg_success_iters = alg_all_iters(alg_conv);    alg_fail_iters = alg_all_iters(~alg_conv);
    alg_all_runtime = [alg.runtime];
    alg_success_runtime = alg_all_runtime(alg_conv);
    alg_fail_runtime = alg_all_runtime(~alg_conv);
    
    alg_list.runs.all(i) = sum([alg.cost]>0);
    alg_list.runs.conv(i) = sum(alg_conv);
    alg_list.iters.all.mean(i) = mean(alg_all_iters);
    alg_list.iters.all.median(i) = median(alg_all_iters);
    alg_list.iters.success.mean(i) = mean(alg_success_iters);
    alg_list.iters.success.median(i) = median(alg_success_iters);
    alg_list.iters.fail.mean(i) = mean(alg_fail_iters);
    alg_list.iters.fail.median(i) = median(alg_fail_iters);
    alg_list.runtime.all.mean(i) = mean(alg_all_runtime);
    alg_list.runtime.all.median(i) = median(alg_all_runtime);
    alg_list.runtime.success.mean(i) = mean(alg_success_runtime);
    alg_list.runtime.success.median(i) = median(alg_success_runtime);
    alg_list.runtime.fail.mean(i) = mean(alg_fail_runtime);
    alg_list.runtime.fail.median(i) = median(alg_fail_runtime);
    alg_list.ips.all.mean(i) = mean(alg_all_iters ./ alg_all_runtime);
    alg_list.ips.all.median(i) = median(alg_all_iters ./ alg_all_runtime);
    
    p = alg_list.runs.conv(i) / alg_list.runs.all(i);
    
    % If no success run is found, estimate the lower bound of time to 2
    % successes.
    if p == 0,
      p = 1 / (alg_list.runs.all(i) + 1);
      alg_list.runtime.success.mean(i) = 0;
    end
    
    % RUSSO-X
    alg_list.RUSSO.optima{i} = nan(alg_list.runs.all(i) - 1, 1);
    alg_list.RUSSO.runtime{i} = nan(alg_list.runs.all(i) - 1, 1);
    if alg_list.runs.conv(i) >= 2,
      % Randperm RUSSO
      if alg_list.runs.all(i) < 10,
        % If the number of runs < 10, compute randperm runs up to 10.
        max_runs = min(10, factorial(alg_list.runs.all(i)));
      else
        max_runs = 100;
      end
      for j = 1 : max_runs
        % RUSSO optimum
        [alg_list.RUSSO.optima{i}(j), alg_list.RUSSO.runtime{i}(j)] = je_compute_russo([alg.cost], [alg.runtime], j, sprintf('tol=%e', opts.tol));
        
        % Triple-RUSSO optimum
        num_extra_russos = 2;
        extra_russo_opts = nan(1, num_extra_russos);
        extra_russo_times = nan(1, num_extra_russos);
        for k = 1 : num_extra_russos
          [extra_russo_opts(k), extra_russo_times(k)] = je_compute_russo([alg.cost], [alg.runtime], j + k * max_runs, sprintf('tol=%e', opts.tol));
        end
        alg_list.triple_RUSSO.optima{i}(j) = min([alg_list.RUSSO.optima{i}(j), extra_russo_opts]);
        alg_list.triple_RUSSO.runtime{i}(j) = sum([alg_list.RUSSO.runtime{i}(j), extra_russo_times]);
      end
      
      % Mean time for RUSSO
      alg_list.TRUSSO.mean(i) = mean(alg_list.RUSSO.runtime{i});
      alg_list.TRUSSO.std(i) = std(alg_list.RUSSO.runtime{i});
      % No replacement RUSSO
      %  for j = 1 : alg_list.runs.all(i) - 1
      %    alg_list.RUSSO.optima{i}(j) = je_compute_russo([alg(j:end).cost]);
      %  end
      % Compute how many have converged to the global minimum.
      alg_list.RUSSO.conv(i) = sum( ...
        abs(alg_list.RUSSO.optima{i} - best_min) < ...
        opts.tol * best_min);
      alg_list.RUSSO.all(i) = sum(~isnan(alg_list.RUSSO.optima{i}));
      alg_list.triple_RUSSO.conv(i) = sum( ...
        abs(alg_list.triple_RUSSO.optima{i} - best_min) < ...
        opts.tol * best_min);
      alg_list.triple_RUSSO.all(i) = sum(~isnan(alg_list.triple_RUSSO.optima{i}));
    end
    
    % TSS
    succ = find(([alg.cost] - best_min) < opts.tol * best_min);
    if length(succ) < 2,
      alg_list.TSS.mean(i) = sum([alg.runtime]);
    else
      k = succ(length(succ) - 1) - 1;
      TSS = zeros(k, 1);
      
      for ii = 1 : k
        next_2s = succ(succ >= ii);
        next_2s = next_2s(1:2);
        TSS(ii) = sum([alg(ii:next_2s(end)).runtime]);
      end
      alg_list.TSS.mean(i) = mean(TSS);
      alg_list.TSS.std(i) = std(TSS);
    end
    
    % TSS estimate using geometric distribution
    alg_list.TSS.mean_est(i) = 2 * alg_list.runtime.success.mean(i);
    % If fail mean is not nan, then add the mean fail time bit.
    if ~isnan(alg_list.runtime.fail.mean(i))
      alg_list.TSS.mean_est(i) = alg_list.TSS.mean_est(i)...
        + 2 * (1 - p) / p * alg_list.runtime.fail.mean(i);
    end
    
    alg_list.TSS.std_est(i) = 0;
    % If fail is not nan, then compute the variance.
    if ~isnan(alg_list.runtime.fail.mean(i))
      alg_list.TSS.std_est(i) = alg_list.runtime.fail.mean(i) * ...
        sqrt(2 * (1 - p)) / p;
    end
    %     MTSSest = alg_list.TSS.mean_est(i);
    %     STSSest = alg_list.TSS.std_est(i);
  end
end

% a = zeros(2, 1);
% a(1) = alg_list.runtime.t2s.mean(i);
% a(2) = alg_list.runtime.t2s.std(i);

if opts.display
  fprintf(['\n--- [ ', dataset,' - rank ', num2str(r), ' ] --- \n']);
end

T = table(...
  alg_names(alg_list.idx, 1), ...
  alg_names(alg_list.idx, 2), ...
  [ alg_list.runs.conv alg_list.runs.all], ...
  alg_list.runs.conv ./ alg_list.runs.all * 100, ...
  alg_list.runtime.all.median, ...
  alg_list.TSS.mean, ...
  alg_list.TSS.std, ...
  [ alg_list.RUSSO.conv alg_list.RUSSO.all], ...
  alg_list.RUSSO.conv ./ alg_list.RUSSO.all * 100, ...
  alg_list.TRUSSO.mean, ...
  alg_list.TRUSSO.std, ...
  [ alg_list.triple_RUSSO.conv alg_list.triple_RUSSO.all], ...
  alg_list.triple_RUSSO.conv ./ alg_list.triple_RUSSO.all * 100, ...
  alg_list.iters.success.median, ...
  alg_list.iters.fail.median, ...
  alg_list.iters.all.median, ...
  alg_list.runtime.success.median, ...
  alg_list.runtime.fail.median, ...
  alg_list.ips.all.median, ...
  alg_list.iters.success.mean, ...
  alg_list.iters.fail.mean, ...
  alg_list.iters.all.mean, ...
  alg_list.runtime.success.mean, ...
  alg_list.runtime.fail.mean, ...
  alg_list.runtime.all.mean, ...
  alg_list.ips.all.mean, ...
  'VariableNames',{'ALGORITHM', 'NOTE', 'CONV_RUNS', 'CONV_RATE', 'TIME_MED', 'M_TSS', 'STD_TSS', 'RUSSO_CONV', 'RUSSO_RATE', 'M_TRUSSO', 'STD_TRUSSO', 'TRIRUSSO_CONV', 'TRIRUSSO_RATE', 'SI_MED', 'FI_MED', 'ITERS_MED', 'ST_MED', 'FT_MED', 'IPS_MED', 'SI_AVR', 'FI_AVR', 'ITERS_AVR', 'ST_AVR', 'FT_AVR', 'TIME_AVR', 'IPS_AVR'});

if opts.display, disp(T(:, [1 2 3 4 5 8 9 10 11 12 13]));
end

if opts.csv
  writetable(T, [file, '.csv']);
end

if opts.display
  fprintf('Best minimum: \t\t%.6f\n', best_min);
  fprintf('Point of convergence: \t%.6f\n\n', best_min + opts.tol * best_min);
end

end

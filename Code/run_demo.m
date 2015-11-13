% Run the setup file prior to performing any experiment.
run setup

% Perform experiment on one algorithm with default settings ...
% [ max_iter = 300, tol = 1e-10, sample_start = 1, sample_end = 100 ]
% We will limit sample_end to 10.
je_run_experiment('dino_trimmed', 4, 'DRW2P=1', 'sample_end=10');

% Perform experiment on multiple algorithms with specific settings
% [ overwrite = 1 will re-run experiment on the specified sample range. ]
je_run_experiment('dino_trimmed', 4, 'DRW2P=1', 'DRW1P=1', 'ALS=1', ...
  'TO_DW=1', 'CH_LM_S_RW2=1', 'NB_RTRMC=1', ... % 'PG_CSF=1'
  'sample_start=1', 'sample_end=3', 'max_iter=300', 'overwrite=1');

% Analyze algorithm performance
% [ csv = 1 will write the csv file to
% Results/<Dataset_name>/<Dataset_name_along with rank>.csv. ]
je_analyze_results('dino_trimmed', 4, 'DRW2P=1', 'DRW1P=1', 'ALS=1', ...
  'TO_DW=1', 'CH_LM_S_RW2=1', 'NB_RTRMC=1', ... % 'PG_CSF=1'
  'csv=1');

%% LIST OF ALGORITHMS
% % Latest custom-coded algorithms
% 'ALS', 'New'         % Alternation (Undamped RW3)
% 'DW', 'New'          % Damped Wiberg (Damped RW2 + projection constraint)
% 'DRW1_UWOFF', '(NO UW)'
% 'DRW1', '(UW)'
% 'DRW1P_UWOFF', '(NO UW)'
% 'DRW1P', '(UW)'   % Damped RW1 with projection constraint (Damped RW1 + projection constraint + retraction)
% 'DRW2_UWOFF', '(NO UW)'
% 'DRW2', '(UW)'
% 'DRW2P_UWOFF', '(NO UW)'
% 'DRW2P', '(UW)'   % Damped RW2 with projection constraint (Damped RW2 + projection constraint + retraction)
% % External algorithms
% 'PG_CSF', 'Orig.'      % Gotardo's CSF (Damped RW2)
% 'TO_DW', 'Orig.'       % Okatani's Damped Wiberg (Damped RW2 + projection constraint)
% 'CH_LM_S', 'Mod.'     % Chen's LM_S (Damped full Newton)
% 'CH_LM_S_GN', 'Mod.'  % Chen's LM_S modified to use GN-matrix (Damped RW1)
% 'CH_LM_S_RW2', 'Mod.' % Chen's LM_S modified to use approximate GN-matrix (Damped RW2)
% 'CH_LM_M', 'Mod.'     % Chen's LM_M (Reduced damped full Newton)
% 'CH_LM_M_GN', 'Mod.'  % Chen's LM_M_GN (Reduced damped RW1)
% 'CH_LM_M_RW2', 'Mod.' % Chen's LM_M modified to use approximate GN-matrix (Reduced damped RW2)
% 'CO_LM_S', 'Orig.'     % Chen's original implementation of LM_S
% 'CO_LM_M', 'Orig.'     % Chen's original implementation of LM_M
% 'CO_LM_M_GN', 'Orig.'  % Chen's original implementation of LM_M_GN
% 'NB_RTRMC', 'Orig.'    % Boumal's original implementation of RTRMC
% 'RC_ALM', 'Unopt.'     % Cabral's ALM (New but not necessarily optimized)
% 'RC_RCALM', 'Unopt.'   % Cabral's ALM with rank continuation (New but not necessarily optimized)
% 'DB_BALM', 'Unopt.'    % Del Bue's BALM (New but not necessarily optimized)
% % Latest Ceres-solver LMs
% 'CE_LM', 'Ceres'       % CERES LM
% 'CE_LMI', 'Ceres'      % CERES LM INNER
% 'CE_ALM', 'Ceres'      % 10 ALS + CE_LM
% 'CE_ALMI', 'Ceres'     % 10 ALS + CE_LM_I
% 'CE_ARULM', 'Ceres'    % 10 reg. ALS + reg. CE_LM + unreg. CE_LM
% 'CE_ARULMI', 'Ceres'   % 10 reg. ALS + reg. CE_LMI + unreg. CE_LMI

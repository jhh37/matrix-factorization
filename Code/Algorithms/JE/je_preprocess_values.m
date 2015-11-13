function [opts, store] = je_preprocess_values(M, W, r, opts)
% JE_PREPROCESS_VALUES
%   Preprocess variables that are required for LRMC-solving algorithms.
%   The processed variables are as follows:
%   - m, n, r
%   - wid
%   - nnz
%   - wnz

%% CUSTOM FUNCTIONS
vec = @(X) X(:);

%% DATASET DIMENSION
% Obtain the size of the dataset.
[store.dim.m, store.dim.n] = size(W);
store.dim.mr = store.dim.m * r;
store.dim.mn = store.dim.m * store.dim.n;

%% WEIGHT MATRIX TYPE AND RANGES
% Vectorize the weight and the measurement values.
[w_vec, ~, ws] = find(W(:));
store.dim.nnz = size(w_vec, 1);

% Store the ranges column-wise (and row-wise if necessary).
store.wc.width = store.dim.n;
store.wc.rng = cumsum([1 sum(W ~= 0, 1)]');
if ~opts.init_v_only
  store.wr.width = store.dim.m;
  store.wr.rng = cumsum([1; sum(W ~= 0, 2)]);
end
% Set the indicator weight flag to 1.
if sum(ws ~= 1) == 0, opts.use_indicator_weight = 1;
end

%% WEIGHTED MEASUREMENT
% Compute m-tilde (Wt * vec(M)).
store.wc.mt = M(w_vec);
if issparse(store.wc.mt), store.wc.mt = full(store.wc.mt);
end
if opts.use_indicator_weight, store.wc.mt = store.wc.mt .* ws;
end

%% ROW-WISE WEIGHT INDICES
if ~opts.init_v_only
  % Obtain the column-to-row indices and save m-tilde in row-major format as well.
  mnt_vec = 1 : store.dim.mn;
  mnt_vec = ceil(mnt_vec / store.dim.m) + store.dim.n * (mod((mnt_vec - 1), store.dim.m));
  [~, store.crid] = sort(mnt_vec(w_vec));
  store.wr.mt = store.wc.mt(store.crid);
  
  if opts.search_unique_weights
    % Obtain the corresponding id of each weight row to avoid repetitions in forming Ut_Q.
    % May need a better method for Netflix. (8~9 sec for AIT8050)
    [W2, ~, store.wr.wid] = unique(W, 'rows');
    
    % Form a vector of weights in ascending order of columns.
    store.wr.unique.count = size(W2, 1); % Store the number of unique weight columns.
    [w_vec, ~, ws] = find(W2');
    store.wr.unique.rng = cumsum([1; sum(W2 ~= 0, 2)]);
  else
    store.wr.wid = 1 : store.dim.m;
    
    % Form a vector of weights in ascending order of columns.
    store.wr.unique.count = store.dim.m; % Store the number of unique weight columns.
    [w_vec, ~, ws] = find(W');
    store.wr.unique.rng = store.wr.rng;
  end
  
  % Make a list of nz for each weight column.
  store.wr.unique.nz = cell(store.wr.unique.count, 1);
  store.wr.unique.nnz = zeros(store.wr.unique.count, 1);
  for j = 1 : store.wr.unique.count
    store.wr.unique.nz{j} = w_vec(store.wr.unique.rng(j) : store.wr.unique.rng(j+1) - 1);
    store.wr.unique.nnz(j) = size(store.wr.unique.nz{j}, 1);
  end
  
  % Also store the individual weight values if a skyline-type weight matrix is used.
  if ~opts.use_indicator_weight
    store.wr.unique.wnz = cell(store.wr.unique.count, 1);
    for j = 1 : store.wr.unique.count
      store.wr.unique.wnz{j} = ws(store.wr.unique.rng(j) : store.wr.unique.rng(j+1) - 1);
    end
  else
  end
end

%% COLUMN-WISE WEIGHT INDICES
% Obtain the corresponding id of each weight column to avoid repetitions in forming Ut_Q.
% May need a better method for Netflix. (8~9 sec for AIT8050)
if opts.search_unique_weights
  [W2, ~, store.wc.wid] = unique(W', 'rows');
  
  % Form a vector of weights in ascending order of columns.
  store.wc.unique.count = size(W2, 1); % Store the number of unique weight columns.
  [w_vec, ~, ws] = find(W2');
  store.wc.unique.rng = cumsum([1; sum(W2 ~= 0, 2)]);
  
else
  store.wc.wid = 1 : store.dim.n;
  
  % Form a vector of weights in ascending order of columns.
  store.wc.unique.count = store.dim.n; % Store the number of unique weight columns.
  [w_vec, ~, ws] = find(W);
  store.wc.unique.rng = store.wc.rng;
end

% Store the assigned columns for each unique weight column.
% This is required for layer-wise computation.
if opts.use_layerwise_computation
  store.wc.unique.ac = cell(store.wc.unique.count, 1);
  for j = 1 : store.wc.unique.count
    store.wc.unique.ac{j} = find(store.wc.wid == j);
  end
end

% Make a list of nz for each weight column.
store.wc.unique.nz = cell(store.wc.unique.count, 1);
for j = 1 : store.wc.unique.count
  store.wc.unique.nz{j} = w_vec(store.wc.unique.rng(j) : store.wc.unique.rng(j+1) - 1);
end
if ~opts.init_v_only,
  store.wc.unique.nnz = zeros(store.wr.unique.count, 1);
  for j = 1 : store.wc.unique.count
    store.wc.unique.nnz(j) = size(store.wc.unique.nz{j}, 1);
  end
end

% Also store the individual weight values if a skyline-type weight matrix is used.
if ~opts.use_indicator_weight
  store.wc.unique.wnz = cell(store.wc.unique.count, 1);
  for j = 1 : store.wc.unique.count
    store.wc.unique.wnz{j} = ws(store.wc.unique.rng(j) : store.wc.unique.rng(j+1) - 1);
  end
else
end

%% Vt indices
% Permutated Vt version 2
if ~opts.use_layerwise_computation
  % If not using the layer-wise computation option store the entire ranges for
  % generating sparse matrices.
  store.vt.irng = bsxfun(@plus, (vertcat(store.wc.unique.nz{store.wc.wid}) - 1) * r, 1:r)';
  store.vt.jrng = bsxfun(@times, 1 : store.dim.nnz, ones(1,r)');
  for i = 1 : store.wc.width
    store.vt.v(store.wc.rng(i) : (store.wc.rng(i + 1) - 1)) = i;
  end
else
  % Otherwise, store in blocks.
  store.vt.irng = cell(store.wc.unique.count, 1);
  store.vt.jrng = cell(store.wc.unique.count, 1);
  for j = 1 : store.wc.unique.count
    store.vt.irng{j} = [vec(bsxfun(@plus, (store.wc.unique.nz{j} - 1) * r, 1:r)'); store.dim.mr];
    store.vt.jrng{j} = [vec(bsxfun(@times, 1 : store.wc.unique.nnz(j), ones(1,r)')); store.wc.unique.nnz(j)];
  end
end
end


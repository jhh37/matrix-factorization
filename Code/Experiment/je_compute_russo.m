function [russo_optimum, russo_time] = je_compute_russo(cost_seq, time_seq, rand_seq, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

opts = au_opts( ...
  'base=1348032014', ...
  'width=10', ...
  'tol=1e-6', ...
  varargin{:});

if nargin > 2,
  if numel(rand_seq) == 1,
    % If sample number given, generate a random permutation sequence.
    rng(opts.base + rand_seq * opts.width, 'twister');
    rand_seq = randperm(numel(cost_seq));
  end
  if numel(cost_seq) == numel(rand_seq),
    % If the actual random permutation sequence is given, use it to sort.
    cost_seq = cost_seq(rand_seq);
    if nargout > 1
      % If we need to compute meantime, then 
      time_seq = time_seq(rand_seq);
    end
  end
end

so_far_min = cost_seq(1);
for n = 2 : numel(cost_seq)
  cost_diff = cost_seq(n) - so_far_min;
  if abs(cost_diff) < opts.tol * so_far_min,
    % If the cost difference is less than the function tolerance value,
    % stop and output russo_opt.
    russo_optimum = so_far_min;
    if nargout > 1, russo_time = sum(time_seq(1 : n));
    end
    return
  elseif cost_diff < 0,
    % Else if the cost difference is less than 0, this means a better
    % optimum is found. Update so_far_min.
    so_far_min = cost_seq(n);
  end
end

% If fail, time = total time taken, the russo optimum is nan.
russo_time = sum(time_seq);
russo_optimum = nan;

end


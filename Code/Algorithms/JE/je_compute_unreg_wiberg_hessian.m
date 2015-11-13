function hess = je_compute_unreg_wiberg_hessian(A, AtQ, AtR, B, R, r, opts, store)
% JE_COMPUTE_UNREG_WIBERG_HESSIAN
%   Compute Hessian of unregularized Wiberg-varian algorithm.

%% HESSIAN FORMING
% Sum over the layers.
% Hess = zeros(store.dim.mr);

%% ALGORITHM III GAUSS-NEWTON MATRIX (Ruhe and Wedin '80)
% Blockwise computation - same as ALS.
% tic
VtT_Vt = cell(store.wr.unique.count, 1);
if opts.use_indicator_weight
  for j = 1 : store.wr.unique.count
    VtT_Vt{j} = B(store.wr.unique.nz{j}, :);
    VtT_Vt{j} = sparse(VtT_Vt{j}' * VtT_Vt{j});
  end
else
  for j = 1 : store.wr.unique.count
    VtT_Vt{j} = bsxfun(@times, store.wr.unique.wnz{j}, B(store.wr.unique.nz{j}, :));
    VtT_Vt{j} = sparse(VtT_Vt{j}' * VtT_Vt{j});
  end
end

hess = blkdiag(VtT_Vt{store.wr.wid});
% toc

%% ALGORITHM II GAUSS-NEWTON MATRIX (Ruhe and Wedin '80)
% Slow for large dataset.
% tic
if opts.alg_type < 3
  % Add n Hessian layers of size mr x mr.
  hess = full(hess);
  % tic
  if opts.use_indicator_weight
    for j = 1 : store.wc.unique.count
      B_j = B(store.wc.unique.ac{j}, :);
      je_add_hessian_layer(hess, AtQ{j} * AtQ{j}', - (B_j' * B_j), int32(store.wc.unique.nz{j} - 1));
    end
  else
    for j = 1 : store.wc.unique.count
      B_j = bsxfun(@times, store.wc.unique.wnz{j}, B_j);
      je_add_hessian_layer(hess, AtQ{j} * AtQ{j}', - (B_j' * B_j), int32(store.wc.unique.nz{j} - 1));
    end
  end
  % toc
  
  %% ALGORITHM I GAUSS-NEWTON MATRIX
  % tic
  if opts.alg_type < 2
    % Add n Hessian layers of size mr x mr.
    if opts.use_indicator_weight
      for j = 1 : store.wc.unique.count
        R_j = R(store.wc.unique.nz{j}, store.wc.unique.ac{j});
        AtRinv_j = eye(r) / AtR{j};
        je_add_hessian_layer(hess, R_j * R_j', AtRinv_j * AtRinv_j', int32(store.wc.unique.nz{j} - 1));
      end
    else
      for j = 1 : store.wc.unique.count
        R_j = bsxfun(@times, store.wc.unique.wnz{j}, R_j);
        AtRinv_j = eye(r) / AtR{j};
        je_add_hessian_layer(hess, R_j * R_j', AtRinv_j * AtRinv_j', int32(store.wc.unique.nz{j} - 1));
      end
    end
  end
end
% toc

% Form the full matrix from its triangular components. If using the
% projection constraint, add it to Hessian.
if strcmpi(opts.constraint, 'projection')
  je_form_hessian(hess, r, A * A');
else je_form_hessian(hess, r);
end
% hess = hess + K * kron(speye(r), A * A') * K';
%   else
%     % If not using layer-wise computation, use sparse matrix %
%     multiplication. This is the fastest option if you have enough memory.
%     AtQ2 = cell(store.wc.unique.count, 1); for j = 1 :
%     store.wc.unique.count
%       AtQ2{j} = sparse(AtQ{j} * AtQ{j}');
%     end % Compute the GN matrix. AtQ22 = blkdiag(AtQ2{store.wc.wid}); if
%     opts.use_indicator_weight
%       BtT = au_sparse(int32(store.vt.irng), int32(store.vt.jrng),
%       B(store.vt.v, :)'); hess = hess - BtT * AtQ22 * BtT';
%     end
end
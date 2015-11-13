function xt_pg_check_buildkron()

try
  xt_pg_kronmex( eye(2), eye(3) );       % simply attempt to run mex-file
catch
  disp('ERROR running "kronmex" mex-file. Will now compile source code:')
  mex xt_pg_kronmex.c
end

end
function data = doForSeveralSimulations(fun, basedir, data, addargs, recursive, reorderfun)
%for each simulation, run a task to process it

if ~exist('recursive', 'var')
  recursive = false;
end

if ~exist('reorderfun', 'var')
  reorderfun = [];
end

if iscell(basedir)
  %batch case
  if not(isempty(reorderfun))
    basedir = basedir(reorderfun(basedir));
  end
  for k=1:numel(basedir)
    data = doForSeveralSimulations(fun, basedir{k}, data, addargs, recursive);
  end
else

  %recursive case
  if recursive
    
    ds = dir(basedir);
    
    ds = ds([ds.isdir]);
    
    nm = {ds.name};
    
    notdot = cellfun(@(x)any(x~='.'), nm);
    
    nm = nm(notdot)';
    
    if not(isempty(reorderfun))
      nm = nm(reorderfun(nm));
    end
    
    for k=1:numel(nm)
      data = doForSeveralSimulations(fun, [basedir filesep nm{k}], data, addargs, recursive);
    end
    
  end
  
  %check if this is a simulation directory
  fls = cellfunc(@(x)[basedir filesep x], {'poblacion.txt' 'arbol.txt' 'timestats.txt' ['..' filesep 'estado.mat'] ['..' filesep 'readableEstado.txt']})';
  theyare = cellfun(@(x)exist(x, 'file'), fls);
%   fprintf('Mira 1: %s\n', any2str(fls));
%   fprintf('Mira 2: %s\n', any2str(theyare));
  if all(theyare)
    fprintf('Processing simulation found in %s\n', basedir);
    data = fun(basedir, data, addargs{:});
  else
    fprintf('     No simulation found in %s\n', basedir);
  end
  
end

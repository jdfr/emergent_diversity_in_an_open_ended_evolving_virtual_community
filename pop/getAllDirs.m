function dirs = getAllDirs(basedir, show)

if nargin<2
  show = false;
end

d = dir(basedir);
istree = any(strcmp('arbol.txt', {d.name}));
if istree
  dirs = {basedir};
else
  ds = d([d.isdir]);
  dsn = {d.name}';
  dsn = dsn(cellfun(@(x)x(1)~='.',dsn));
  if isempty(dsn)
    dirs = {};
  else
    dirs = cellfun(@(x)getAllDirs([basedir filesep x], false), dsn, 'uniformoutput', false);
    dirs = vertcat(dirs{:});
  end
end
dirs = dirs(:);

if show
  fprintf('dirs = { ...\n');
  for k=1:numel(dirs)
    fprintf('  %s ; ...\n', any2str(dirs{k}));
  end
  fprintf('  };\n');
end

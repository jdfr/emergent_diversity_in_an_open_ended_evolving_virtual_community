function dr = findDirWW(basedir, with, without)

d = dir(basedir);
dn = {d.name}';
[w1 w2] = cellfun(@(x)deal(strcmp(x, with), strcmp(x, without)), dn);

if any(w1) && not(any(w2))
  %fprintf('FOUND: %s\n', basedir);
  dr1 = {basedir};
else
  dr1 = {};
end

dn = dn([d.isdir]');
dn = dn(cellfun(@(x)x(1)~='.', dn));
%fprintf('%s: Searcing in %d\n', basedir, numel(dn));
if isempty(dn)
  dr = dr1;
else
  dr = cell(size(dn));
  for k=1:numel(dn)
    dr{k} = findDirWW([basedir filesep dn{k}], with, without);
  end
  dr = vertcat(dr1, dr{:});
end
%fprintf('RES: %s\n', any2str(dr));



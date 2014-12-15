function [varargout] = generateSimulationCompressions(basedir, poph, maxg, level)

name = [basedir filesep sprintf('compression_%d', level) '.mat'];
if exist(name, 'file')
  fprintf('ALREADY CALCULATED!!!!');
  return
end

if isempty(poph)
  if exist([basedir filesep 'poph.mat'], 'file')
    fprintf('LOADIN MAT...\n');
    poph = load([basedir filesep 'poph.mat']);
    poph = poph.poph;
  else
    fprintf('LOADIN SIM...\n');
    poph = loadSimulation(basedir);
  end
  fprintf('LOADIN POPH DONE\n');
end

%WE RELY ON ThESE TWO FACTS TO EASE COMPUTATIONS:
if not(all(poph.idx==poph.tree.idxD(1:numel(poph.idx))))
  error('INDEXES INCOSISTENT FOR %s!!!!', basedir);
end
if not(all(poph.idx(poph.tree.parents(poph.tree.parents>0))==poph.tree.idxA(poph.tree.parents>0)))
  error('PARENT INDEXES INCOSISTENT FOR %s!!!!', basedir);
end


compression = nan(size(poph.generation));

ming = 0;

if isempty(maxg)
  maxg = max(poph.generation);
end

parents = poph.tree.parents;
changes = poph.tree.change;
minG = poph.minGenome;


for g=ming:maxg
  thisg = find(poph.generation==g);
  ng = numel(thisg);
  %include individuals without parents (actually, this only the individual
  %in the first generation)
  notToCalculate = (parents(thisg)~=0);
  %include individuals which have not changed
  notToCalculate(notToCalculate) = cellfun(@(x)x(1)=='=', changes(thisg(notToCalculate)));
  %include individuals whose minGenomes are different
  notToCalculate(notToCalculate) = cellfun(@isequal, minG(thisg(notToCalculate)), minG(thisg(notToCalculate)));
  compression(thisg(notToCalculate)) = compression(parents(thisg(notToCalculate)));
  inds = not(notToCalculate);
  thisg = thisg(inds);
  fprintf('GENERATION %04d/%04d %04d/%04d inds\n', g, maxg, numel(thisg), ng);
  if not(isempty(thisg))
    inds = find(inds);
    nm = sprintf('P%03g', g);
    P = load([basedir filesep nm]);
    P = P.(nm);
    for k=1:numel(inds)
      compression(thisg(k)) = getCompression(P.raster{inds(k)}, P.dimensions(inds(k),:), level);
    end
  end
  clear raster;
end

if nargout>0
  varargout{1} = compression;
end

compression = struct('compression', {compression}, 'level', {level}); %#ok<NASGU>

save(name, 'compression');



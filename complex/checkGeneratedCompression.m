function checkGeneratedCompression(basedir, poph, compression, level)
%to check the result of generateSimulationCompressions

if isempty(compression)
  name = [basedir filesep sprintf('compression_%d', level) '.mat'];
  fprintf('LOADING %s\n', name);
  compression = load(name);
  compression = compression.compression;
  level = compression.level;
  compression = compression.compression;
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

ming = 0;
maxg = poph.generation(numel(compression));

for g=ming:maxg
  thisg = find(poph.generation==g);
  fprintf('GENERATION %d (%d)\n', g, numel(thisg));
  inds = poph.individual(thisg);
  nm = sprintf('P%03g', g);
  P = load([basedir filesep nm]);
  P = P.(nm);
  for k=1:numel(thisg)
    c = getCompression(P.raster{inds(k)}, P.dimensions(inds(k),:), level);
    if c~=compression(thisg(k))
      error('For G=%d,I=%d(idx=%d), recorded=%s, calculated=%s', g, inds(k), thisg(k), mat2str(compression(thisg(k))), mat2str(c));
    end
  end
end
fprintf('NO ERRORS!!!!!!!!!!!!!!!!!!!!!!!');

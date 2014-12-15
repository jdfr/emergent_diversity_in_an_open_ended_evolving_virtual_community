function makeMaxHistograms(basedir)
%gráficas para comparar entre sí varias simulaciones (situadas en un mismo
%superdirectorio común)

filestats = [basedir filesep 'statsend.mat'];
if exist(filestats, 'file')
  load(filestats, 'stats');
else
  stats = getAllDirs(basedir);
  save(filestats, 'stats');
end

maxgens                      = cellfun(@(x)x.gen,             stats);
maxHeightInLastGeneration    = cellfun(@(x)x.max.height,      stats);
maxWidthInLastGeneration     = cellfun(@(x)x.max.width,       stats);
maxNumGInLastGeneration      = cellfun(@(x)x.max.numG,        stats);
maxLongGInLastGeneration     = cellfun(@(x)x.max.longGenome,  stats);
maxBranchesInLastGeneration  = cellfun(@(x)x.max.nbranches,   stats);

meanHeightInLastGeneration    = cellfun(@(x)x.mean.height,      stats);
meanWidthInLastGeneration     = cellfun(@(x)x.mean.width,       stats);
meanNumGInLastGeneration      = cellfun(@(x)x.mean.numG,        stats);
meanLongGInLastGeneration     = cellfun(@(x)x.mean.longGenome,  stats);
meanBranchesInLastGeneration  = cellfun(@(x)x.mean.nbranches,   stats);

figargs  = {'Visible', 'off', 'Units', 'pixels', 'PaperUnits', 'points', 'PaperPosition', [0 0 1000 1000]};
lineargs = {'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'k'};

fig = figure(figargs{:});
hist(maxgens, numel(maxgens));
grid on; title('max generation histogram');
drawnow; saveas(fig,[basedir filesep 'maxgen_hist.png'],'png'); close(fig);

plotCloud(maxgens, maxHeightInLastGeneration,    'max',  'height',     basedir, figargs, lineargs);
plotCloud(maxgens, maxWidthInLastGeneration,     'max',  'width',      basedir, figargs, lineargs);
plotCloud(maxgens, maxNumGInLastGeneration,      'max',  'countG',     basedir, figargs, lineargs);
plotCloud(maxgens, maxLongGInLastGeneration,     'max',  'longGenome', basedir, figargs, lineargs);
plotCloud(maxgens, maxBranchesInLastGeneration,  'max',  'nbranches',  basedir, figargs, lineargs);

plotCloud(maxgens, meanHeightInLastGeneration,   'mean', 'height',     basedir, figargs, lineargs);
plotCloud(maxgens, meanWidthInLastGeneration,    'mean', 'width',      basedir, figargs, lineargs);
plotCloud(maxgens, meanNumGInLastGeneration,     'mean', 'countG',     basedir, figargs, lineargs);
plotCloud(maxgens, meanLongGInLastGeneration,    'mean', 'longGenome', basedir, figargs, lineargs);
plotCloud(maxgens, meanBranchesInLastGeneration, 'mean', 'nbranches',  basedir, figargs, lineargs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotCloud(maxgens, measure, aggregatename, fieldname, basedir, figargs, lineargs)
fig = figure(figargs{:});
line(maxgens, measure, lineargs{:});
grid on;
xlabel('last generation');
ylabel([aggregatename ' ' fieldname ' in last generation']);
title(['Aggregate point cloud: last generation / ' aggregatename ' ' fieldname ' in last generation']);
drawnow;
saveas(fig, [basedir filesep 'gen_' aggregatename '_' fieldname 'InMaxGen.png'], 'png');
close(fig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stats = getAllDirs(basedir)
fprintf('Probing dir <%s>...\n', basedir);
dirr        = dir(basedir);
subdirs     = arrayfun(@(x) x.isdir && (x.name(1)~='.'), dirr);
stats       = cell(0,1);
if     exist([basedir filesep 'arbol.txt'], 'file') && ...
       exist([basedir filesep 'poblacion.txt'], 'file') && ...
       exist([basedir filesep '..' filesep 'estado.mat'], 'file') 
  fprintf('Compiling statistics for dir <%s>...\n', basedir);
  try
    stats = {getStatistics(basedir)};
    stats{1}.basedir = basedir;
    fprintf('DONE! Statistics compiled for dir <%s>...\n', basedir);
  catch ME
    fprintf('HOW SAD!!!! An Error has been generated:\n');
    fprintf(showError(ME));
  end
end
if any(subdirs)
  subdirs = dirr(subdirs);
  for k=1:numel(subdirs)
    substats = getAllDirs([basedir filesep subdirs(k).name]);
    if ~isempty(substats)
      stats = [stats; substats]; %#ok<AGROW>
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function st = getStatistics(basedir)
t1=cputime;
poph = loadSimulation(basedir);
t2=cputime;
[poph.longGenome, poph.numG] = cellfun(@(x)deal(numel(x), sum(x=='G')), poph.genome);
t3=cputime;
fprintf('TIME 1: %s, TIME 2: %s\n', mat2str(t2-t1), mat2str(t3-t2));
maxgen    = max(poph.generation);
lastGen   = find(poph.generation==maxgen);

funcs       = {'max', 'mean', 'median'};
funcsargs   = reshape([funcs; cell(size(funcs))], 1, []);
hands       = cellfun(@str2func, funcs, 'uniformoutput', false);
fields      = {'height', 'width', 'nbranches', 'nleafs', 'nleafsOnTop', 'numG', 'longGenome'};
fieldsargs  = reshape([fields; cell(size(fields))], 1, []);

allMax = struct(fieldsargs{:});
for k=1:numel(fields)
  allMax.(fields{k}) = max(poph.(fields{k}));
end

st = struct('gen', maxgen, funcsargs{:}, 'allMax', allMax);
for m=1:numel(funcs)
  for k=1:numel(fields)
    st.(funcs{m}).(fields{k}) = hands{m}(poph.(fields{k})(lastGen)); %#ok<FNDSB>
  end
end




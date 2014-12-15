function makeAllGraphics(arg1, toFile, fieldForBest)
%gráficas de estadísticas a lo largo de una simulación

if ~exist('toFile', 'var')
  toFile = true;
end

if ~exist('fieldForBest', 'var')
  fieldForBest = 'nbranches';%'fitness';
%         'fitness',       (datos{4}), ...
%         'xpos',          (datos{5}), ...
%         'height',        (datos{6}), ...
%         'width',         (datos{7}), ...
%         'nbranches',     (datos{8}), ...
%         'nleafs',        (datos{9}), ...
%         'nPixelsInSoil', (datos{10}), ...
%         'nleafsOnTop',   (datos{11}), ...
end

%load data if necessary
if ischar(arg1)
  basedir = arg1;
  %load simulation data
  poph = loadSimulation(basedir);
elseif isstruct(arg1)
  poph = arg1;
  basedir = '.';
elseif iscell(arg1)
  poph    = arg1{1};
  basedir = arg1{2};
else
  error('unrecognized first argument!!');
end

poph.fitness(poph.fitness<0) = -1;

mingen = min(poph.generation);
maxgen = max(poph.generation);
ngens  = maxgen-mingen+1; %#ok<NASGU>

params = putDefaultParameters(poph.ACIParams);

%compute best (at the end) individual's observable statistics
lastGen                      = find(poph.generation==max(poph.generation));
[nevemind, idxbest]          = max(poph.(fieldForBest)(lastGen)); %#ok<ASGLU>
idxbest                      = lastGen(idxbest);
individualBest               = poph.individual(idxbest);
fprintf('Printed individual for best <%s>: generation=%d, individual=%d\n', fieldForBest, poph.generation(lastGen(1)), individualBest);
[absBest, relBest, genoBest] = compileIndividualHistory(poph, params, poph.generation(lastGen(1)), 1, individualBest, [min(poph.generation) poph.generation(lastGen(1))]); %#ok<ASGLU>

%compute observable variables' aggregate statistics
[agreg agregenome agregdisrup popsize] = agregateObservablesByGeneration(poph, params);

%%%%%%%%%%%%%%%%%%%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gs = (mingen:maxgen)';

if toFile
  figargs = {'Visible', 'off', 'Units', 'pixels', 'PaperUnits', 'points', 'PaperPosition', [0 0 1000 1000]};
else
  figargs = {};
end

ejes01 = [mingen maxgen 0 1.01];

%enviroment's width
envwidth = zeros(ngens, 1);
g=mingen;
for k=1:ngens
  thisGen = poph.xpos(poph.generation==g);
  zzz = max(thisGen)-min(thisGen)+1;
  envwidth(k) = zzz;
  g=g+1;
end
clear thisGen;
toPlot = {...
          {'environment''s width',   'environment''s width along generations', ...
            {envwidth         'b'}...
          }};
plotStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, 'envWidth');

toPlot = {...
          {'fitness',        'best (at last generation) individual''s fitness history', ...
            {agreg.fitness.median,   'g', 'population median'; ...
             absBest.fitness         'b', 'best individual''s'}...
          } ...
          {'disruption',     'best (at last generation) individual''s disruption    history', ...
            {agregdisrup.disruption.median, 'g', 'population median'; ...
             absBest.disruption             'b', 'best individual''s'} ...
          ejes01} ...
          {@()plotOperators(gs, genoBest.changes)} ...
          {'CEILING fitness',          'best (at last generation) individual''s CEILING fitness (maximum possible fitness without competition) history', ...
            {agreg.ceilFitness.median,     'g', 'population median'; ...
             absBest.ceilFitness,          'b', 'best individual''s'}, ...
          } ...
          {'RATIO fitness',         'best (at last generation) individual''s fitness RATIO (actual fitness divided by ceiling fitness) history', ...
            {agreg.ratioFitness.median,    'g', 'population median'; ...
             absBest.ratioFitness,         'b', 'best individual''s'}...
          ejes01} ...
          };
plotStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, 'bestRecordingFitness');

%population size
toPlot = {...
          {'population size',   'population size along generations', ...
            {popsize         'b', 'best individual''s'}...
          }};
plotStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, 'popSize');

%fitness
toPlot = {{'max'     'b', 'max fitness'} ...
          {'mean',   'g', 'mean fitness'} ...
          {'median', 'r', 'median fitness'} ...
          {'std',    'c', 'fitness'' std'} ...
          };
plotAggregateStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, agreg, 'fitness', 'Fitness statistics by generation');
toPlot = {{'max'     'b', 'max CEILING fitness'} ...
          {'mean',   'g', 'mean CEILING fitness'} ...
          {'median', 'r', 'median CEILING fitness'} ...
          {'std',    'c', 'CEILING fitness'' std'} ...
          };
plotAggregateStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, agreg, 'ceilFitness', 'CEILING fitness (maximum possible fitness without competition) statistics by generation');
toPlot = {{'max'     'b', 'max fitness RATIO'} ...
          {'mean',   'g', 'mean fitness RATIO'} ...
          {'median', 'r', 'median fitness RATIO'} ...
          {'std',    'c', 'fitness RATIO''s std'} ...
          {'min',    'm', 'min fitness RATIO'} ...
          };
plotAggregateStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, agreg, 'ratioFitness', 'Fitness RATIO (actual fitness divided by ceiling fitness) statistics by generation', ejes01);
  
%width
toPlot = {{'max'     'b', 'max width'} ...
          {'mean',   'g', 'mean width'} ...
          {'median', 'r', 'median width'} ...
          {'std',    'c', 'width''s std'} ...
          {'min',    'm', 'min width'} ...
          };
plotAggregateStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, agreg, 'width', 'Width statistics by generation');

%coverage
toPlot = {{'max'     'b', 'max coverage'} ...
          {'mean',   'g', 'mean coverage'} ...
          {'median', 'r', 'median coverage'} ...
          {'std',    'c', 'coverage''s std'} ...
          {'min',    'm', 'min coverage'} ...
          };
plotAggregateStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, agreg, 'coverage', 'Coverage (potential light catch divided by width) statistics by generation', ejes01);

%height
toPlot = {{'max'     'b', 'max height'} ...
          {'mean',   'g', 'mean height'} ...
          {'median', 'r', 'median height'} ...
          {'std',    'c', 'height''s std'} ...
          };
plotAggregateStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, agreg, 'height', 'Height statistics by generation');

%genome length
toPlot = {{'max'     'b', 'max genome length'} ...
          {'mean',   'g', 'mean genome length'} ...
          {'median', 'r', 'median genome length'} ...
          {'std',    'c', 'genome length''s std'} ...
          {'min',    'm', 'min genome length'} ...
          };
plotAggregateStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, agregenome, 'long', 'Genome legth statistics by generation');

%disruption
toPlot = {{'mean',   'b', 'mean disruption'} ...
          {'median', 'g', 'median disruption'} ...
          {'std',    'r', 'disruption''s std'} ...
          };
plotAggregateStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, agregdisrup, 'disruption', 'Disruption statistics by generation');

%last best individual's recordings
toPlot = {...
          {'fitness',        'best (at last generation) individual''s fitness history', ...
            {agreg.fitness.median,   'g', 'population median'; ...
             absBest.fitness         'b', 'best individual''s'}...
          } ...
          {'disruption',     'best (at last generation) individual''s disruption    history', ...
            {agregdisrup.disruption.median, 'g', 'population median'; ...
             absBest.disruption             'b', 'best individual''s'} ...
          ejes01} ...
          {'width',          'best (at last generation) individual''s width   history', ...
            {agreg.width.median,     'g', 'population median'; ...
             absBest.width,          'b', 'best individual''s'}, ...
          } ...
          {'height',         'best (at last generation) individual''s height  history', ...
            {agreg.height.median,    'g', 'population median'; ...
             absBest.height,         'b', 'best individual''s'}...
          } ...
          {'coverage',         'best (at last generation) individual''s coverage  history', ...
            {agreg.coverage.median,    'g', 'population median'; ...
             absBest.coverage,         'b', 'best individual''s'}...
          ejes01} ...
          };
plotStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, 'bestRecordingPhenotype');
toPlot = {...
          {'fitness',        sprintf('best (at last generation, in <%s>) individual''s fitness history, gen=%d, ind=%d', fieldForBest, poph.generation(lastGen(1)), individualBest), ...
            {agreg.fitness.median,   'g', 'population median'; ...
             absBest.fitness         'b', 'best individual''s'}...
          } ...
          {'disruption',     'best (at last generation) individual''s disruption    history', ...
            {agregdisrup.disruption.median, 'g', 'population median'; ...
             absBest.disruption             'b', 'best individual''s'} ...
          ejes01} ...
          {'genome length',  'best (at last generation) individual''s genome length history', ...
            {agregenome.long.median, 'g', 'population median'; ...
             genoBest.long           'b', 'best individual''s'} ...
          } ...
          {'number of G symbols',  'best (at last generation) individual''s G count history', ...
            {agregenome.numG.median, 'g', 'population median'; ...
             genoBest.numG           'b', 'best individual''s'} ...
          } ...
          {'tree depth',  'best (at last generation) individual''s G tree depth of genome history', ...
            {agregenome.treeDepth.median, 'g', 'population median'; ...
             genoBest.treeDepth           'b', 'best individual''s'} ...
          } ...
          };
plotStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, 'bestRecordingGenome');

plotHistogramByGeneration(basedir, toFile, figargs, 'fitnessHistogram', poph, mingen, maxgen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotAggregateStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, structure, field, titulo, ejes)
fig        = figure(figargs{:});
hold on;
for k=1:numel(toPlot);
  [xs ys] = stairs(gs, structure.(field).(toPlot{k}{1}));
  line(xs, ys, 'Color', toPlot{k}{2});
  clear xs ys;
end
xlim([gs(1) gs(end)]);
legendargs = cellfun(@(x)x{end}, toPlot, 'uniformoutput', false);
legend(legendargs, 'Location', 'BestOutside');
title(titulo);
xlabel('generation');
ylabel(field);
if exist('ejes', 'var')
  axis(ejes);
end
grid on;
if toFile
  drawnow;
  saveas(fig,[basedir filesep 'graph_' field 'ByGen.png'],'png');
  close(fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotStatisticsByGeneration(basedir, toFile, figargs, toPlot, gs, fileName)
fig        = figure(figargs{:});
for k=1:numel(toPlot)
  subplot(numel(toPlot),1,k);
  hold on;
  if isa(toPlot{k}{1}, 'function_handle')
    toPlot{k}{1}();
  else
    for m=1:size(toPlot{k}{3}, 1)
      [xs ys] = stairs(gs, toPlot{k}{3}{m,1});
      line(xs, ys, 'Color', toPlot{k}{3}{m,2});
      clear xs ys;
    end
    xlim([gs(1) gs(end)]);
    if size(toPlot{k}{3}, 2)>=3
      legend(toPlot{k}{3}(:,3), 'Location', 'BestOutside');
    end
    xlabel('generations');
    ylabel(toPlot{k}{1});
    title(toPlot{k}{2});
    grid on;
    if size(toPlot{k},2)>=4
      axis(toPlot{k}{4});
    end
  end
end
if toFile
  drawnow;
  saveas(fig,[basedir filesep 'graph_' fileName 'ByGen.png'],'png');
  close(fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotHistogramByGeneration(basedir, toFile, figargs, fileName, poph, mingen, maxgen)
fitnesses        = poph.fitness;
generations      = poph.generation;
fig = figure(figargs{:});
hold on;
colors = cool(maxgen-mingen+1);
kk=1;
fitnessBinEdges    = [-inf 0 linspace(0, max(fitnesses), 20)]';
fitnessBinEdges(3) = eps(1);
fitnessBinCenters  = [-1; 0; (fitnessBinEdges(3:end-1)+diff(fitnessBinEdges(3:end))/2)];
for k=mingen:maxgen
  thisGen = generations==k;
  n = reshape(histc(fitnesses(thisGen), fitnessBinEdges), [], 1);
  n = [n(1:end-2); (n(end-1)+n(end))];
  [xs ys] = stairs(fitnessBinCenters,n);
  line(xs, ys, 'Color', colors(kk,:));
  clear xs ys;
  kk=kk+1;
end
grid on;
xlabel('fitness');
ylabel('counts');
legend(arrayfun(@(x)sprintf('gen. %d', x), mingen:maxgen, 'uniformoutput', false), 'FontSize', 5);%, 'Location', 'Best');
title('fitness histogram by generation');
if toFile
  drawnow;
  saveas(fig,[basedir filesep 'graph_' fileName 'ByGen.png'],'png');
  close(fig);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [agreg agregenome agregdisrup popsize]= agregateObservablesByGeneration(poph, params)
%compile and show statistics by generation

mingen = min(poph.generation);
maxgen = max(poph.generation);

ngens       = maxgen-mingen+1;
popsize     = zeros(ngens,1);
mea         = {struct('max', zeros(ngens,1), 'min', zeros(ngens, 1), 'std', zeros(ngens,1), 'mean', zeros(ngens,1), 'median', zeros(ngens,1))};
agreg       = struct('fitness', mea, 'height', mea, 'width', mea, 'nbranches', mea, 'nleafs', mea, 'nleafsOnTop', mea);
agregenome  = struct('long', mea, 'numG', mea, 'treeDepth', mea);
agregdisrup = struct('disruption', mea, 'distOrigen', mea);
disrupSinCeros = [true false];
vars        = fieldnames(agreg);
arevars     = cellfun(@(x)ismember(x, fieldnames(poph)), vars);
disrups     = fieldnames(agregdisrup);
aredisrups  = cellfun(@(x)ismember(x, fieldnames(poph.tree)), disrups);
meas        = fieldnames(mea{1});
nchar       = 0;
g           = mingen;

agreg.coverage     = mea{1};
agreg.ceilFitness  = mea{1};
agreg.ratioFitness = mea{1};
thereAreLeafs  = isfield(poph, 'nleafsOnTop');

for k=1:ngens
  thisGen                    = poph.generation==g;
  thisGen(thisGen)           = poph.nPixelsInSoil(thisGen)==0;
  popsize(k)                 = sum(thisGen);
  if popsize(k)>0
    for n=1:numel(vars)
      if arevars(n)
        var                    = vars{n};
        thisGenVar             = poph.(var)(thisGen);
        for m=1:numel(meas)
          mea                  = meas{m};
          agreg.(var).(mea)(k) = feval(mea, thisGenVar);
        end
      end
    end
    if thereAreLeafs
      coverage = poph.nleafsOnTop(thisGen)./poph.width(thisGen);
      coverage(isnan(coverage)) = 0;
      ceilFitness  = calculateCeilFitness(params, poph.height(thisGen), poph.nleafsOnTop(thisGen), poph.nbranches(thisGen), poph.nleafs(thisGen), poph.nPixelsInSoil(thisGen));
      ratioFitness = poph.fitness(thisGen)./ceilFitness;
      ratioFitness(isnan(ratioFitness)) = 0;
      for m=1:numel(meas)
        mea                         = meas{m};
        agreg.coverage.(mea)(k)     = feval(mea, coverage);
        agreg.ceilFitness.(mea)(k)  = feval(mea, ceilFitness);
        agreg.ratioFitness.(mea)(k) = feval(mea, ratioFitness);
      end
    end
    [longGenomes, numGsGenomes, treeDepths] = cellfun(@(x)deal(numel(x), sum(x=='G'), safeMax(cumsum((x=='[')-(x==']')))), poph.genome(thisGen));
    for m=1:numel(meas)
      mea                           = meas{m};
      agregenome.long.(mea)(k)      = feval(mea, longGenomes);
      agregenome.numG.(mea)(k)      = feval(mea, numGsGenomes);
      agregenome.treeDepth.(mea)(k) = feval(mea, treeDepths);
    end
  end
  thisGenTree                     = poph.tree.gD==g;
  if any(thisGenTree)
    for n = 1:numel(disrups)
      if aredisrups(n)
        disrup                         = disrups{n};
        thisGenDisrup                  = poph.tree.(disrup)(thisGenTree);
        if disrupSinCeros(n)
          thisGenDisrup                = thisGenDisrup(thisGenDisrup~=0);
        end
        if numel(thisGenDisrup)>0
          for m=1:numel(meas)
            mea                          = meas{m};
            agregdisrup.(disrup).(mea)(k)   = feval(mea, thisGenDisrup);
          end
        end
      end
    end
  end
  if nchar>0
    fprintf(repmat('\b', 1, nchar));
  end
  str   = sprintf('Processed generation %04d', g);
  fprintf(str);
  nchar = numel(str);
  g=g+1;
end
if nchar>0
  fprintf(repmat('\b', 1, nchar));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ceilFitness = calculateCeilFitness(params, height, nleafsOnTop, nbranches, nleafs, negY)
haveNotTooManyNeg  = negY<=params.negYThreshold;
if (~params.tooTallArePlaced) && (~isempty(params.heightT))
  areNotTooTall    = height<=params.heightT;
else
  areNotTooTall    = true(size(haveNotTooManyNeg));
end
haveSomeFitness    = haveNotTooManyNeg & areNotTooTall;
ceilFitness        = zeros(size(haveNotTooManyNeg));
ceilFitness(haveSomeFitness) = nleafsOnTop(haveSomeFitness);
switch params.divideFitness
  case 'byLeafs'
    % dividiendo por el numero de hojas del arbol
    do       = nleafs>0; % |G|= 10
    ceilFitness(do)  = (ceilFitness(do)-negY(do))./((nleafs(do)*params.factorG).^params.alphaG);
    ceilFitness(~do) = 0;
  case 'byBranches'
    % dividiendo por el numero de ramas distinguibles del arbol *|G|
    do               = nbranches>0; % |G|= 10
    ceilFitness(do)  = (ceilFitness(do)-negY(do))./((nbranches(do)*params.factorG).^params.alphaG);
    ceilFitness(~do) = 0;
end
ceilFitness(ceilFitness<0) = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mx = safeMax(x)
if isempty(x)
  mx = 0;
else
  mx = max(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [absolute, relative, genomes] = compileIndividualHistory(poph, params, generation, rangeid, individual, range)
%print individual statistics (fitness, etc.) throughout evolutionary
%history

gens          = max(min(poph.generation), range(1)):min(max(poph.generation), range(2));
absolute.fitness       = repmat(nan, size(gens));
absolute.height        = repmat(nan, size(gens));
absolute.width         = repmat(nan, size(gens));
absolute.nbranches     = repmat(nan, size(gens));
absolute.nleafs        = repmat(nan, size(gens));
absolute.nleafsOnTop   = repmat(nan, size(gens));
relative.fitness       = repmat(nan, size(gens));
relative.height        = repmat(nan, size(gens));
relative.width         = repmat(nan, size(gens));
relative.nbranches     = repmat(nan, size(gens));
relative.nleafs        = repmat(nan, size(gens));
relative.nleafsOnTop   = repmat(nan, size(gens));
genomes.long           = repmat(nan, size(gens));
genomes.numG           = repmat(nan, size(gens));
genomes.treeDepth      = repmat(nan, size(gens));
genomes.changes        = repmat(false, 7,numel(gens));
changesStr             = cell(1, numel(gens));
%absolute values whose relative ones have not been calculated (because
%either it is impossible, or too complicated, or I'm tired)
absolute.nPixelsInSoil = repmat(nan, size(gens));
absolute.disruption    = repmat(nan, size(gens));
absolute.distOrigen    = repmat(nan, size(gens));
absolute.coverage      = repmat(nan, size(gens));
absolute.ceilFitness   = repmat(nan, size(gens));
absolute.ratioFitness  = repmat(nan, size(gens));

names = fieldnames(absolute);

arenames = cellfun(@(x)ismember(x, fieldnames(poph)), names);

indiv = poph.triplet2idx(generation, rangeid, individual);

nchar = 0;

fillDisruption = isfield(poph.tree, 'disruption');
fillDistOrigen = isfield(poph.tree, 'distOrigen');
thereAreLeafs  = isfield(poph, 'nleafsOnTop');

while ~isempty(indiv)
  idx = find(poph.idx==indiv);
  
  if isempty(idx)
    break;
  else
    gen                = find(gens==poph.generation(idx));
    
    for k=1:numel(names)
      if arenames(k)
        name = names{k};
        absolute.(name)(gen)       = poph.(name)(idx);
        thisGen                    = poph.generation==gens(gen);
        mn                         = min(poph.(name)(thisGen));
        range                      = max(poph.(name)(thisGen))-mn;
        if range==0
          relative.(name)(gen)     = 0.5;
        else
          relative.(name)(gen)     = (absolute.(name)(gen)-mn)/range;
        end
      end
    end
    genomes.long(gen)            = numel(poph.genome{idx});
    genomes.numG(gen)            = sum(poph.genome{idx}=='G');
    genomes.treeDepth(gen)       = max(cumsum((poph.genome{idx}=='[')-(poph.genome{idx}==']')));
    inTree                       = find(poph.tree.idxD==indiv);
    if fillDisruption
      absolute.disruption(gen)   = poph.tree.disruption(inTree);
    end
    if fillDistOrigen
      absolute.distOrigen(gen)   = poph.tree.distOrigen(inTree);
    end
    genomes.changes(:,gen) = ismember('ARLCDIT', poph.tree.change{inTree})';
    changesStr{gen}        = poph.tree.change{inTree};
    indiv                        = poph.tree.idxA(inTree);
  end
  if nchar>0
    fprintf(repmat('\b', 1, nchar));
  end
  str   = sprintf('Processed generation %04d', gens(gen));
  fprintf(str);
  nchar = numel(str);
end
if nchar>0
  fprintf(repmat('\b', 1, nchar));
end
if thereAreLeafs
  absolute.coverage     = absolute.nleafsOnTop./absolute.width;
  absolute.coverage(isnan(absolute.coverage)) = 0;
  absolute.ceilFitness  = calculateCeilFitness(params, absolute.height, absolute.nleafsOnTop, absolute.nbranches, absolute.nleafs, absolute.nPixelsInSoil);
  absolute.ratioFitness = absolute.fitness./absolute.ceilFitness;
  absolute.ratioFitness(isnan(absolute.ratioFitness)) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotOperators(gs, changes)

%     1 2   op_alteracion        A simbolo  posicion
%     2 3   op_dupaleatoria      R segmento posicion
%     3 4   op_dupnivel          L segmento posicion
%     4 5   op_dupsecuencia      C segmento posicion
%     5 6   op_eliminacion       D posicion
%     6 7   op_insercion         I simbolo  posicion
%     7 8   op_transferencia     T segmento posicion

[filas columnas] = find(changes);

columnas = reshape(gs(columnas), [], 1);

           %1 2 3 4 5 6 7
newfilas = [7 2 3 4 5 6 1];
newfilas = reshape(newfilas(filas), [], 1);

axis([gs(1) gs(end) 0 7.1]);
line(columnas, newfilas, 'LineStyle', 'none', 'Marker', 'd', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');%, 'MarkerSize', 1);

xlabel('generations')
ylabel('operators')
set(gca,'YTick',0:7);
set(gca,'YTickLabel','None|op_transferencia|op_dupaleatoria|op_dupnivel|op_dupsecuencia|op_eliminacion|op_insercion|op_alteracion');
grid on;
subplot(5,1,2); position1 = get(gca, 'Position');
subplot(5,1,3); position2 = get(gca, 'Position');
pause(0.01);
set(gca, 'Position', [position2(1:2) position1(3) position2(4)]);


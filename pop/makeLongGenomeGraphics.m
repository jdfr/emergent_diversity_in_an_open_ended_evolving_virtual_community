function makeLongGenomeGraphics(arg1, gen, toFile, basedir)
%gráficas de scatterplot para una determinada generación: correlación entre
%longitud de genoma y otras variables

%load data if necessary
if ischar(arg1)
  basedir = arg1;
  %load simulation data
  poph = loadSimulation(basedir);
elseif isstruct(arg1)
  poph = arg1;
elseif iscell(arg1)
  poph = arg1{1};
  basedir = arg1{2};
else
  error('unrecognized first argument!!');
end

if ~exist('basedir', 'var')
  basedir = '.';
end

poph.fitness(poph.fitness<0) = -1;

if isinf(gen)
  ThisGen = find(poph.generation==max(poph.generation));
else
  ThisGen = find(poph.generation==gen);
end


longGenomeThisGen     = cellfun(@numel, poph.genome(ThisGen));
numGThisGen           = cellfun(@(x)sum(x=='G'), poph.genome(ThisGen));
widthThisGen          = poph.width(ThisGen);
heightThisGen         = poph.height(ThisGen);
nbranchesThisGen      = poph.nbranches(ThisGen);
fitnessThisGen        = poph.fitness(ThisGen);
nleafsThisGen         = poph.nleafs(ThisGen);
nleafsOnTopThisGen    = poph.nleafsOnTop(ThisGen);
xposThisGen           = poph.xpos(ThisGen);

if toFile
  figArgs = {'Visible', 'off', 'Units', 'pixels', 'PaperUnits', 'points', 'PaperPosition', [0 0 1000 1000]};
else
  figArgs = {};
end

lineArgs = {'LineStyle', 'none', 'Marker', '.'};

plotScatter(longGenomeThisGen, numGThisGen,           'numG',            'scatterplot: genome length vs. amount of G symbols',      toFile, basedir, figArgs, lineArgs, gen);
plotScatter(longGenomeThisGen, nbranchesThisGen,      'numBranches',     'scatterplot: genome length vs. amount of branches',       toFile, basedir, figArgs, lineArgs, gen);
plotScatter(longGenomeThisGen, nleafsThisGen,         'numLeafs',        'scatterplot: genome length vs. amount of leafs',          toFile, basedir, figArgs, lineArgs, gen);
plotScatter(longGenomeThisGen, nleafsOnTopThisGen,    'numLeafsOnTop',   'scatterplot: genome length vs. amount of leafs on top',   toFile, basedir, figArgs, lineArgs, gen);
plotScatter(longGenomeThisGen, heightThisGen,         'height',          'scatterplot: genome length vs. tree height',              toFile, basedir, figArgs, lineArgs, gen);
plotScatter(longGenomeThisGen, widthThisGen,          'width',           'scatterplot: genome length vs. tree width',               toFile, basedir, figArgs, lineArgs, gen);
plotScatter(longGenomeThisGen, xposThisGen,           'position',        'scatterplot: genome length vs. tree placement',           toFile, basedir, figArgs, lineArgs, gen);
plotScatter(longGenomeThisGen, fitnessThisGen,        'fitness',         'scatterplot: genome length vs. tree fitness',             toFile, basedir, figArgs, lineArgs, gen);

fig = figure(figArgs{:});
hist(longGenomeThisGen, 100);
title('histogram: genome length');
grid on;
if toFile
  saveas(fig, [basedir filesep 'genomeLengthHISTOGRAM.png'], 'png');
  close(fig);
end

function plotScatter(longGenomeThisGen, measure, ylab, titulo, toFile, basedir, figArgs, lineArgs, gen)
fig = figure(figArgs{:});
line(longGenomeThisGen, measure, lineArgs{:});
grid on;
xlabel('genome length');
ylabel(ylab);
title(titulo);
if toFile
  if isinf(gen)
    lastgen = '_AtLastGen';
  else
    lastgen = ['_AtGen' mat2str(gen)];
  end
  saveas(fig, [basedir filesep 'genomeLengthVS' ylab lastgen '.png'], 'png');
  close(fig);
end
